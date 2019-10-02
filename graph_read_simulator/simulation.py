import sys
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")
from pyfaidx import Fasta
from offsetbasedgraph import NumpyIndexedInterval
from numpy.random import randint
import numpy as np
from simple_read_mutator import Mutator
np.random.seed(1)


def simulate_reads(chromosome, haplotype_fasta_file_name, haplotype_interval_file_name,
                   haplotype_reference_interval_file_name, coverage=15, read_length=150,
                   snv_prob=0.01, deletion_prob=0.001, insertion_prob=0.001):

    ref = Fasta(haplotype_fasta_file_name)
    logging.info("Chromosome length: %s" % len(ref[chromosome]))
    chrom_length = len(ref[chromosome])
    chrom_min = 0
    chrom_max = chrom_length - read_length - 10

    n_reads = int(coverage * chrom_length / read_length)
    logging.info("Will simulate %d reads to get coverage %.3f on chromosome %s" %
    (n_reads, coverage, chromosome))

    logging.info("Reading intervals")
    haplotype_interval = NumpyIndexedInterval.from_file(haplotype_interval_file_name)
    haplotype_nodes = haplotype_interval.nodes_in_interval()

    ref_interval = NumpyIndexedInterval.from_file(haplotype_reference_interval_file_name)
    linear_ref_nodes = ref_interval.nodes_in_interval()
    mutator = Mutator(read_length + 10)

    logging.info("Starting simulation")
    i = 0

    while i < n_reads:

        start = randint(chrom_min, chrom_max)
        end = start + read_length + 10
        seq = str(ref[chromosome][start:end])
        if "n" in seq:
            continue

        i += 1

        haplotype_end_node = haplotype_interval.get_node_at_offset(end)
        haplotype_node_offset = haplotype_interval.get_node_offset_at_offset(start)
        haplotype_node = haplotype_interval.get_node_at_offset(start)
        if haplotype_node not in linear_ref_nodes:
            # On a variant, check if next
            while True:
                haplotype_node += 1
                if haplotype_node in linear_ref_nodes:
                    break

        ref_offset = ref_interval.get_offset_at_node(haplotype_node) + haplotype_node_offset

        # Hacky way to get shared nodes between read and linear ref
        read_haplotype_nodes = set(list(range(haplotype_node, haplotype_end_node + 1))).intersection(haplotype_nodes)
        n_nodes_not_in_linear_ref = len(read_haplotype_nodes) - len(read_haplotype_nodes.intersection(linear_ref_nodes))

        mutated_seq = mutator.mutate_sequence(seq, snv_prob, deletion_prob, insertion_prob)
        mutated_seq = mutated_seq[0:read_length]

        print("%s %d %d %s" % (chromosome, ref_offset, n_nodes_not_in_linear_ref, mutated_seq))

