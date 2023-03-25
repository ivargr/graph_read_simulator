import logging
logging.basicConfig(level=logging.INFO)
import argparse
import dataclasses
import sys
from .util import make_haplotype_fasta, vg_path_to_obg_interval, make_haplotype_paths
from .simulation import simulate_reads
from .id_assignment import assign_ids
from .diploid_reference_builder import DiploidReferenceBuilder
import random
from pyfaidx import Fasta
from .vcf_simulation import VcfSimulator
import numpy as np
from .vcf_simulation import Variants
from .simulation import MultiChromosomeCoordinateMap
from shared_memory_wrapper import from_file, to_file
import bionumpy as bnp


def liftover(args):
    coordinate_map = from_file(args.coordinate_map)
    assert isinstance(coordinate_map, MultiChromosomeCoordinateMap)
    reverse = args.reverse

    file = bnp.open(args.input, buffer_type=bnp.Bed6Buffer)
    out_file = bnp.open(args.output, "w", buffer_type=bnp.Bed6Buffer)

    for chunk in file.read_chunks():
        chromosomes = chunk.chromosome.tolist()
        starts = chunk.start
        stops = chunk.stop

        new_starts = np.array(
            [coordinate_map.convert(chromosome, start, reverse=reverse) for chromosome, start in zip(chromosomes, starts)]
        )
        new_stops = new_starts + (stops-starts)

        chunk = dataclasses.replace(chunk, start=new_starts, stop=new_stops)
        out_file.write(chunk)


def assign_ids_wrapper(args):
    assign_ids(args.truth_file_name, args.fasta_file_name)


def simulate_reads_new_wrapper(args):
    chromosome = args.chr_haplotype.split()[0]
    haplotype = args.chr_haplotype.split()[1]
    random_seed = random.randint(0, 2**32-1)

    simulate_reads(chromosome, haplotype, args.coverage,
                   args.read_length, args.snv_prob, args.deletion_prob, args.insertion_prob,
                   random_seed, data_base_name=args.data_base_name)


def simulate_reads_wrapper(args):
    chromosome = args.chr_haplotype.split()[0]
    haplotype = args.chr_haplotype.split()[1]
    random_seed = int(haplotype)

    haplotype_interval_file_name = "haplotype_" + chromosome + "__" + haplotype + ".intervalcollection.indexed"
    haplotype_fasta_file_name = "haplotype_" + chromosome + "__" + haplotype + ".fasta"
    haplotype_reference_interval_file_name = "individual_reference_path_" + chromosome + ".intervalcollection.indexed"

    simulate_reads(chromosome, haplotype_fasta_file_name, haplotype_interval_file_name,
                   haplotype_reference_interval_file_name, args.coverage,
                   args.read_length, args.snv_prob, args.deletion_prob, args.insertion_prob,
                   random_seed)


def vg_path_to_obg_interval_wrapper(args):
    vg_path_to_obg_interval(args.path_file_name, args.out_file_name)


def make_haplotype_paths_wrapper(args):
    make_haplotype_paths(args.graph_file_name, args.linear_ref_path_file_name, args.haplotype0_file_name,
                         args.haplotype1_file_name, args.out_base_name, args.chromosome)


def prepare_simulation(args):
    builder = DiploidReferenceBuilder(args.reference, args.vcf, args.chromosome, args.haplotype, args.base_output_name)
    builder.build()


def simulate_from_linear_ref(args):
    from .simulation import _simulate, OneToOneCoordinateMap
    ref = Fasta(args.reference_fasta)
    coordinate
    _simulate(ref, args.chromosome, OneToOneCoordinateMap(), args.coverage, args.read_length,
              include_reverse_complements=False, skip_reads_with_missing_bases=False)



def sample_from_linear_ref(args):
    ref = Fasta(args.reference_fasta)
    spacing = args.spacing
    read_length = args.read_length
    chromosomes = list(ref.keys())
    logging.info("Will sample from chromosomes %s" % chromosomes)

    for chromosome in chromosomes:
        ref_sequence = ref[chromosome]
        for start in range(0, len(ref_sequence)-read_length, spacing):
            if start % 1000000 == 0:
                logging.info("On chromosome %s, pos %d/%d" % (chromosome, start, len(ref_sequence)))

            end = start + read_length
            seq = ref_sequence[start:end]
            print("%s %d %d %s" % (chromosome, start, 0, seq))


def main():
    run_argument_parser(sys.argv[1:])


def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='Graph Read Simulator',
        prog='graph_read_simulator',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers()

    # Store alignments
    store = subparsers.add_parser("vg_path_to_obg_interval")
    store.add_argument("path_file_name")
    store.add_argument("out_file_name")
    store.set_defaults(func=vg_path_to_obg_interval_wrapper)

    # make_haplotype_paths
    command = subparsers.add_parser("make_haplotype_paths")
    arguments = "graph_file_name, linear_ref_path_file_name, haplotype0_file_name, haplotype1_file_name, out_base_name, chromosome"
    for arg in arguments.split(", "):
        command.add_argument(arg)
    command.set_defaults(func=make_haplotype_paths_wrapper)

    # prepare simulation
    command = subparsers.add_parser("prepare_simulation")
    command.add_argument("--chromosome", "-c", required=True)
    command.add_argument("--haplotype", "-a", help="Either 0 or 1", type=int, required=True)
    command.add_argument("--vcf", "-v", required=True)
    command.add_argument("--reference", "-r", required=True)
    command.add_argument("--base-output-name", "-o", required=False, )
    command.set_defaults(func=prepare_simulation)

    # simulate reads
    command = subparsers.add_parser("simulate_reads")
    command.add_argument("chr_haplotype", help="String of chromosome and haplotype separated by space")
    command.add_argument("coverage", type=float)
    command.add_argument("--read_length", "-r", type=int, default=150, required=False)
    command.add_argument("--snv_prob", "-s", type=float, default=0.01, required=False)
    command.add_argument("--deletion_prob", "-d", type=float, default=0.001, required=False)
    command.add_argument("--insertion_prob", "-i", type=float, default=0.001, required=False)
    command.add_argument("--data-base-name", "-D", help="Base file name for data, from prepare_simulation")
    command.set_defaults(func=simulate_reads_new_wrapper)

    # simulate reads using multiple threads
    command = subparsers.add_parser("simulate_reads_multithread")
    command.add_argument("--chromosomes", "-c", help="Comma-separated list of chromosomes")
    command.add_argument("--data-base-name", "-D", help="Base file name for data, from prepare_simulation")
    command.add_argument("--coverage", "-C", type=float)
    command.add_argument("--read_length", "-r", type=int, default=150, required=False)
    command.add_argument("--snv_prob", "-s", type=float, default=0.01, required=False)
    command.add_argument("--deletion_prob", "-d", type=float, default=0.001, required=False)
    command.add_argument("--insertion_prob", "-i", type=float, default=0.001, required=False)
    command.set_defaults(func=simulate_reads_new_wrapper)


    # simulate from linear ref
    command = subparsers.add_parser("simulate_from_linear_ref")
    command.add_argument("chromosome")
    command.add_argument("coverage", type=float)
    command.add_argument("-f", "--reference-fasta", required=True)
    command.add_argument("--read_length", "-r", type=int, default=150, required=False)
    command.add_argument("--snv_prob", "-s", type=float, default=0.01, required=False)
    command.add_argument("--deletion_prob", "-d", type=float, default=0.001, required=False)
    command.add_argument("--insertion_prob", "-i", type=float, default=0.001, required=False)
    command.set_defaults(func=simulate_from_linear_ref)


    command = subparsers.add_parser("sample_from_linear_ref")
    command.add_argument("-f", "--reference-fasta", required=True)
    command.add_argument("-s", "--spacing", required=False, type=int, default=5)
    command.add_argument("--read_length", "-r", type=int, default=150, required=False)
    command.set_defaults(func=sample_from_linear_ref)

    # chip-seq
    def simulate_chip_seq(args):
        from .chip_seq import ChipSeqSimulator
        simulator = ChipSeqSimulator(args.chromosome, 0, args.chromosome_size, args.n_peaks, args.read_length, args.fragment_length,
                                     average_fragments_per_peak=args.average_fragments_per_peak,
                                     average_fragments_per_peak_std=args.average_fragments_per_peak_std,
                                     snv_error_prob=args.snv_prob,
                                     deletion_error_prob=args.deletion_prob,
                                     insertion_error_prob=args.insertion_prob,
                                     noise_coverage=args.noise_coverage)
        simulator.simulate()

    command = subparsers.add_parser("simulate_chip_seq")
    command.add_argument("--chromosome", "-c", required=True)
    command.add_argument("--n_peaks", "-n", type=int, default=210, required=False)
    command.add_argument("--fragment_length", "-f", type=int, default=210, required=False)
    command.add_argument("--read_length", "-r", type=int, default=76, required=False)
    command.add_argument("--snv_prob", "-s", type=float, default=0.01, required=False)
    command.add_argument("--deletion_prob", "-d", type=float, default=0.001, required=False)
    command.add_argument("--insertion_prob", "-i", type=float, default=0.001, required=False)
    command.add_argument("--chromosome_size", "-z", type=float, default=1000, required=False)
    command.add_argument("--average_fragments_per_peak", "-a", type=int, default=15, required=False)
    command.add_argument("--average_fragments_per_peak_std", "-A", type=int, default=3, required=False)
    command.add_argument("--noise-coverage", "-N", type=float, default=0.1, required=False)
    command.set_defaults(func=simulate_chip_seq)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    # assign ids
    command = subparsers.add_parser("assign_ids")
    command.add_argument("truth_file_name")
    command.add_argument("fasta_file_name")
    command.set_defaults(func=assign_ids_wrapper)

    def simulate_population_vcf(args):
        np.random.seed(1)
        from .vcf_simulation import GenotypeMatrixSimulator
        ref = str(Fasta(args.reference)["1"])
        genotype_matrix = GenotypeMatrixSimulator(args.n_individuals*2, args.n_variants, correlation_rate=0.95, mutation_rate=0.002)
        simulated_genotype_matrix = genotype_matrix.subsample(args.n_individuals)
        logging.info("Subsampled genotype matrix to %d individuals" % args.n_individuals)

        simulator = VcfSimulator(ref, simulated_genotype_matrix)
        variants = simulator.simulate()
        variants.to_vcf_file(args.out_file_name, args.header_file)
        logging.info("Wrote population to %s" % args.simulate_individual)

        if args.simulate_individual is not None:
            individual = VcfSimulator.numeric_genotypes_to_literal(genotype_matrix.simulate_single_individual_not_in_matrix())
            for i, variant in enumerate(variants):
                variant.set_genotypes([individual[i]])

            Variants(variants).to_vcf_file(args.simulate_individual)
            logging.info("Wrote individual to %s" % args.simulate_individual)


    # Simulate a population vcf
    command = subparsers.add_parser("simulate_population_vcf")
    command.add_argument("-r", "--reference", required=True)
    command.add_argument("-n", "--n-variants", default=50, type=int, required=False)
    command.add_argument("-i", "--n-individuals", default=50, type=int, required=False)
    command.add_argument("-I", "--simulate-individual", required=False, help="If specified, simulate an individual, and write to this file")
    command.add_argument("-o", "--out-file-name")
    command.add_argument("-H", "--header-file", required=False, help="Use vcf header from this file")
    command.set_defaults(func=simulate_population_vcf)

    def simulate_individual_vcf(args):
        import obgraph
        from obgraph.variants import VcfVariants
        from .vcf_simulation import make_simulated_individual_genotypes_from_genotype_matrix, Variant, Variants
        np.random.seed(args.random_seed)
        logging.info("Using random seed %d" % args.random_seed)
        genotype_matrix = obgraph.genotype_matrix.GenotypeMatrix.from_variants(
            obgraph.variants.VcfVariants.from_vcf(args.population_vcf))
        genotype_matrix = genotype_matrix.matrix
        simulated_genotypes = make_simulated_individual_genotypes_from_genotype_matrix(genotype_matrix)

        population_variants = VcfVariants.from_vcf(args.population_vcf)
        individual_variants = []
        for i, variant in enumerate(population_variants):
            individual_variants.append(Variant(1, variant.position, variant.ref_sequence, variant.variant_sequence, [simulated_genotypes[i]]))

        Variants(individual_variants).to_vcf_file(args.out_file_name)


    # simulate individual vcf from population vcf
    command = subparsers.add_parser("simulate_individual_vcf")
    command.add_argument("-v", "--population-vcf", required=True)
    command.add_argument("-o", "--out-file-name")
    command.add_argument("-s", "--random-seed", type=int, default=1, required=True)
    command.set_defaults(func=simulate_individual_vcf)

    cmd = subparsers.add_parser("liftover")
    cmd.add_argument("-i", "--input", required=True)
    cmd.add_argument("-c", "--coordinate-map", required=True)
    cmd.add_argument("-o", "--output", required=True)
    cmd.add_argument("-r", "--reverse", required=False, type=bool)
    cmd.set_defaults(func=liftover)



    args = parser.parse_args(args)
    args.func(args)



