import logging
logging.basicConfig(level=logging.INFO)
import argparse
import sys
from .util import make_haplotype_fasta, vg_path_to_obg_interval, make_haplotype_paths
from .simulation import simulate_reads
from .id_assignment import assign_ids


def assign_ids_wrapper(args):
    assign_ids(args.truth_file_name, args.fasta_file_name)


def simulate_reads_wrapper(args):
    chromosome = args.chr_haplotype.split()[0]
    haplotype = args.chr_haplotype.split()[1]
    random_seed = int(haplotype)

    haplotype_interval_file_name = "haplotype_" + chromosome + "__" + haplotype + ".intervalcollection.indexed"
    haplotype_fasta_file_name = "haplotype_" + chromosome + "__" + haplotype + ".fasta"
    haplotype_reference_interval_file_name = "individual_reference_path_" + chromosome + ".intervalcollection.indexed"

    if args.repeat_mask != None:
        repeat_mask_file_name = "rep_chr" + chromosome + ".npy"
    else:
        repeat_mask_file_name = None

    simulate_reads(chromosome, haplotype_fasta_file_name, haplotype_interval_file_name,
                   haplotype_reference_interval_file_name, args.coverage,
                   args.read_length, args.snv_prob, args.deletion_prob, args.insertion_prob,
                   repeat_mask_file_name, random_seed)


def vg_path_to_obg_interval_wrapper(args):
    vg_path_to_obg_interval(args.path_file_name, args.out_file_name)


def make_haplotype_paths_wrapper(args):
    make_haplotype_paths(args.graph_file_name, args.linear_ref_path_file_name, args.haplotype0_file_name,
                         args.haplotype1_file_name, args.out_base_name, args.chromosome)


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


    # Simulate reads
    command = subparsers.add_parser("simulate_reads")
    command.add_argument("chr_haplotype", help="String of chromosome and haplotype separated by space")
    command.add_argument("coverage", type=float)
    command.add_argument("--read_length", "-r", type=int, default=150, required=False)
    command.add_argument("--snv_prob", "-s", type=float, default=0.01, required=False)
    command.add_argument("--deletion_prob", "-d", type=float, default=0.001, required=False)
    command.add_argument("--insertion_prob", "-i", type=float, default=0.001, required=False)
    command.add_argument("--repeat_mask", "-m", type=str, nargs="?", default=None, required=False)
    command.set_defaults(func=simulate_reads_wrapper)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    # assign ids
    command = subparsers.add_parser("assign_ids")
    command.add_argument("truth_file_name")
    command.add_argument("fasta_file_name")
    command.set_defaults(func=assign_ids_wrapper)

    args = parser.parse_args(args)
    args.func(args)



