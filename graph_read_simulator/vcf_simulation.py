import logging
import numpy as np
import random

class Variant:
    def __init__(self, chromosome, position, ref_sequence, alt_sequence, genotypes=None):
        self._chromosome = chromosome
        self._position = position
        self._ref_sequence = ref_sequence
        self._alt_sequence = alt_sequence
        self._genotypes = genotypes

    def set_genotypes(self, genotypes):
        self._genotypes = genotypes

    def to_vcf_line(self):
        genotypes = "\t".join(self._genotypes)
        return "%d\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n" % \
               (self._chromosome, self._position, self._ref_sequence, self._alt_sequence, genotypes)

class Variants:
    def __init__(self, variants):
        self._variants = variants

    def to_vcf_file(self, file_name, use_header_from_file=None):
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

        # add chrom line
        individual_names = "\t".join(["SAMPLE" + str(d) for d in range(len(self._variants[0]._genotypes))])
        header_lines += "#CRHOM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % individual_names

        f = open(file_name, "w")
        f.write(header_lines)
        for variant in self._variants:
            f.writelines([variant.to_vcf_line()])
        f.close()


class VcfSimulator:
    def __init__(self, reference_sequence, n_variants, n_individuals):
        self._ref = reference_sequence
        self._n_individuals = n_individuals
        self._n_variants = n_variants
        self._genome_size = len(self._ref)
        self._variants = []

    @staticmethod
    def _random_nucleotides(cls, length=1, except_nucleotide=None):
        nucleotides = ["A", "C", "T", "G"]
        if except_nucleotide is not None:
            nucleotides.remove(except_nucleotide.upper())
        return ''.join(np.random.choice(nucleotides, length))

    def _simulate_variant_at_position(self, position):
        if np.random.random() < 0.9:
            # SNP
            ref_sequence = self._ref[position]
            alt_sequence = VcfSimulator._random_nucleotides(ref_sequence, 1)
        else:
            # insertion
            ref_sequence = self._ref[position]
            alt_sequence = ref_sequence + self._random_nucleotides(ref_sequence, np.random.randint(1, 5))

        variant = Variant(1, position+1, ref_sequence, alt_sequence)
        self._variants.append(variant)

    def _numeric_genotypes_to_literal(self, genotypes):
        map = {1: "0|0", 2: "1|1", 3: "0|1"}
        return [map[genotype] for genotype in genotypes]

    def _simulate_genotypes(self, correlation_between_variants=0.8):
        matrix = np.zeros((self._n_individuals, self._n_variants))
        for individual in range(self._n_individuals):
            # at first variant, give random genotype
            # at next variant, there is a certain prob that individual will have previous genotype % 3
            for variant in range(self._n_variants):
                if variant == 0:
                    genotype = np.random.randint(1, 4)
                else:
                    if np.random.random() < correlation_between_variants:
                        genotype = 1 + (matrix[individual, variant - 1] % 3)
                    else:
                        genotype = np.random.randint(1, 4)
                matrix[individual, variant] = genotype

        for variant in range(self._n_variants):
            self._variants[variant].set_genotypes(self._numeric_genotypes_to_literal(matrix[:,variant]))

    def simulate(self):
        variant_positions = np.random.randint(0, self._genome_size, self._n_variants)
        for pos in variant_positions:
            self._simulate_variant_at_position(pos)

        self._simulate_genotypes()
        return Variants(self._variants)










