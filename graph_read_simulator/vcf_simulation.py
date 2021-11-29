import logging
import numpy as np
import random
import obgraph
from obgraph.variants import VcfVariants

default_vcf_header = \
    """##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##contig=<ID=1>
"""


def make_simulated_individual_genotypes_from_genotype_matrix(genotype_matrix, consistency_rate=0.96):
    # follows an individual in matrix with some consistency
    n_individuals, n_variants = genotype_matrix.shape
    genotypes = []
    individual = None
    for variant in range(n_variants):
        # find new random individual at first variant and with a prob elsewhere
        if variant == 0 or random.random() < consistency_rate:
            individual = np.random.randint(0, n_individuals)

        genotypes.append(genotype_matrix[individual, variant])

    genotypes = VcfSimulator.numeric_genotypes_to_literal(genotypes)

    return genotypes


class Variant:
    def __init__(self, chromosome, position, ref_sequence, alt_sequence, genotypes=None):
        self._chromosome = chromosome
        self._position = position
        self._ref_sequence = ref_sequence
        self._alt_sequence = alt_sequence
        self._genotypes = genotypes

    def set_genotypes(self, genotypes):
        self._genotypes = genotypes

    def allele_frequency(self):
        frequency = 0
        for genotype in self._genotypes:
            frequency += genotype.count("1")
        return frequency / (2*len(self._genotypes))

    def to_vcf_line(self):
        genotypes = "\t".join(self._genotypes)
        return "%d\t%d\t.\t%s\t%s\t.\tPASS\tAF=%.2f\tGT\t%s\n" % \
               (self._chromosome, self._position, self._ref_sequence, self._alt_sequence, self.allele_frequency(), genotypes)

class Variants:
    def __init__(self, variants):
        self._variants = variants

    def _sort_variants(self):
        logging.info("Sorting variants")
        self._variants = sorted(self._variants, key=lambda v: v._position)

    def __getitem__(self, item):
        return self._variants[item]

    def to_vcf_file(self, file_name, use_header_from_file=None, skip_variants_with_N=True):
        self._sort_variants()
        header_lines = ""
        if use_header_from_file is not None:
            logging.info("Using header from %s" % use_header_from_file)
            f = open(use_header_from_file)
            header_lines = ""
            for line in f:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        break
                    header_lines += line
        else:
            header_lines = default_vcf_header

        # add chrom line
        individual_names = "\t".join(["SAMPLE" + str(d) for d in range(len(self._variants[0]._genotypes))])
        header_lines += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % individual_names

        f = open(file_name, "w")
        f.write(header_lines)
        for variant in self._variants:
            if "N" in variant._ref_sequence:
                continue
            f.writelines([variant.to_vcf_line()])
        f.close()


class GenotypeMatrixSimulator:
    def __init__(self, n_individuals, n_variants, correlation_rate=0.95, mutation_rate=0.01):
        assert n_individuals
        self._n_individuals = n_individuals
        self._n_variants = n_variants
        self._correlation_rate = correlation_rate
        self._mutation_rate = mutation_rate
        self.matrix = None
        self.simulate_genotypes_realistically()

    def _random_sparse_genotype_matrix(self, n_individuals, nonzero_prob=0.15):
        matrix = np.zeros((n_individuals, self._n_variants)) + 1
        for individual in range(n_individuals):
            for variant in range(self._n_variants):
                if np.random.random() < nonzero_prob:
                    if np.random.random() < 0.33:
                        matrix[individual, variant] = 2
                    else:
                        matrix[individual, variant] = 3

        return matrix

    def simulate_single_individual_not_in_matrix(self):
        return self._genotype_matrix_from_exisisting(self.matrix, 1)[0]

    def subsample(self, n_individuals):
        rows = np.arange(0, self.matrix.shape[0])
        random_rows = np.random.choice(rows, n_individuals, replace=False)
        return self.matrix[random_rows, :]

    def _genotype_matrix_from_exisisting(self, existing_matrix, n_new_individuals):
        # follow existing individuals with some prob
        new_matrix = np.zeros((n_new_individuals, self._n_variants))
        n_existing_individuals = existing_matrix.shape[0]
        logging.info("Genotyping %d individuals from existing matrix with %d individuals" % (n_new_individuals, n_existing_individuals))
        for i in range(n_new_individuals):
            following_individual = np.random.randint(0, n_existing_individuals)
            for variant in range(self._n_variants):
                if np.random.random() < self._mutation_rate:
                    # random mutation
                    genotype = np.random.randint(1, 4)
                else:
                    # follow existing individual
                    genotype = existing_matrix[following_individual, variant]

                # some possibility of changing individual
                if np.random.random() > self._correlation_rate:
                    following_individual = np.random.randint(0, n_existing_individuals)

                new_matrix[i, variant] = genotype
        return new_matrix

    def simulate_genotypes_realistically(self):
        # simulate first a few indidividuals, then simulate more from them, and so on
        n_individuals = [max(2, self._n_individuals // 20), max(2, self._n_individuals // 15), self._n_individuals // 3]
        n_individuals.append(self._n_individuals-sum(n_individuals))
        logging.info("Will simulate in batches of individuals: %s" % n_individuals)

        n_tot = sum(n_individuals)
        matrix = np.zeros((n_tot, self._n_variants)) + 1
        n_individuals_done = 0
        for i, n in enumerate(n_individuals):
            if i == 0:
                matrix[0:n, :] = self._random_sparse_genotype_matrix(n)
            else:
                # make new random individuals based on the previous
                logging.info("Making %d new individuals from exissting %d" % (n, n_individuals_done))
                matrix[n_individuals_done:n_individuals_done + n, :] = self._genotype_matrix_from_exisisting(
                    matrix[0:n_individuals_done, :], n)

            n_individuals_done += n

        self.matrix = matrix


class VcfSimulator:
    def __init__(self, reference_sequence, genotype_matrix):
        self._ref = reference_sequence
        self._genome_size = len(self._ref)
        self._genotype_matrix = genotype_matrix
        self._n_individuals, self._n_variants = self._genotype_matrix.shape
        logging.info("Simulating population vcf with %d individuals and %d variants" % (self._n_individuals, self._n_variants))
        self._variants = []


    @staticmethod
    def random_nucleotides(length=1, except_nucleotide=None):
        nucleotides = ["A", "C", "T", "G"]
        if except_nucleotide is not None and except_nucleotide.upper() in nucleotides:
            nucleotides.remove(except_nucleotide.upper())
        return ''.join(np.random.choice(nucleotides, length))

    def _simulate_variant_at_position(self, position):
        if np.random.random() < 0.7:
            # SNP
            ref_sequence = self._ref[position]
            alt_sequence = VcfSimulator.random_nucleotides(1, ref_sequence)
        else:
            # 50/50 insertion/deletion
            if np.random.randint(0, 2) == 0:
                size = np.random.randint(2, 6)
            else:
                size = np.random.randint(2, 40)

            if np.random.random() < 0.5:
                ref_sequence = ''.join(self._ref[position:position+size])
                alt_sequence = self._ref[position]
            else:
                # insertion
                ref_sequence = self._ref[position]
                alt_sequence = ref_sequence + self.random_nucleotides(size, ref_sequence)

            assert ref_sequence.upper() != alt_sequence.upper()

        variant = Variant(1, position+1, ref_sequence, alt_sequence)
        self._variants.append(variant)

    @staticmethod
    def numeric_genotypes_to_literal(genotypes):
        map = {1: "0|0", 2: "1|1", 3: "0|1"}
        return [map[int(genotype)] for genotype in genotypes]


    def simulate(self):
        variant_positions = np.random.choice(np.arange(10, self._genome_size-10), self._n_variants, replace=False)  # choice to avoid duplicate positions
        assert len(set(variant_positions)) == self._n_variants, "Got only %d variants" % len(set(variant_positions))
        for pos in variant_positions:
            self._simulate_variant_at_position(pos)

        for variant in range(self._n_variants):
            self._variants[variant].set_genotypes(VcfSimulator.numeric_genotypes_to_literal(self._genotype_matrix[:,variant]))

        return Variants(self._variants)










