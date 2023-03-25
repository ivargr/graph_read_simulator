import sys
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")
from pyfaidx import Fasta
from numpy.random import randint
import numpy as np
from simple_read_mutator import Mutator
from Bio.Seq import Seq


class OneToOneCoordinateMap:
    def __init__(self):
        pass

    def convert(self, coordinate):
        return coordinate

class MultiChromosomeCoordinateMap:
    def __init__(self, coordinate_maps: dict):
        self._coordinate_maps = coordinate_maps

    def convert(self, chromosome, haplotype_coordinate, reverse=False):
        if reverse:
            return self._coordinate_maps[chromosome].convert_reverse(haplotype_coordinate)
        return self._coordinate_maps[chromosome].convert(haplotype_coordinate)

    def haplotype_has_variant_between(self, chromosome, start, end):
        return self._coordinate_maps[chromosome].haplotype_coordinate_exists_between(start, end)

    def reference_has_variant_between(self, chromosome, start, end):
        return self._coordinate_maps[chromosome].reference_coordinate_exists_between(start, end)


class CoordinateMap:
    def __init__(self, reference, haplotype):
        self.reference = reference
        self.haplotype = haplotype

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data["reference"], data["haplotype"])

    def convert(self, haplotype_coordinate):
        index = np.searchsorted(self.haplotype, haplotype_coordinate)
        # Find index of closest reference coordinate
        diff_before = haplotype_coordinate - self.haplotype[index-1]
        try:
            diff_after = haplotype_coordinate - self.haplotype[index]
        except IndexError:
            logging.error("Index %d failed. N indexes in haplotype: %d. Haplotype coordinate: %d" % (index, len(self.haplotype), haplotype_coordinate))
            raise

        if abs(diff_before) < abs(diff_after):
            index -= 1
            diff = diff_before
        else:
            diff = diff_after

        corresponding_ref_coordinate = self.reference[index]
        adjusted_ref_coordinate = corresponding_ref_coordinate + diff

        return adjusted_ref_coordinate

    def convert_reverse(self, reference_coordinate):
        index = np.searchsorted(self.reference, reference_coordinate)
        # Find index of closest reference coordinate
        diff_before = reference_coordinate - self.reference[index - 1]
        diff_after = reference_coordinate - self.reference[index]

        if abs(diff_before) < abs(diff_after):
            index -= 1
            diff = diff_before
        else:
            diff = diff_after

        corresponding_haplotype_coordinate = self.haplotype[index]
        adjusted_haplotype_coordinate = corresponding_haplotype_coordinate + diff

        return adjusted_haplotype_coordinate

    def reference_has_variant_between(self, start, end):
        return self.reference_coordinate_exists_between(start, end)

    def haplotype_has_variant_between(self, start, end):
        return self.haplotype_coordinate_exists_between(start, end)

    def haplotype_coordinate_exists_between(self, start, end):
        index_start = np.searchsorted(self.haplotype, start)
        index_end = np.searchsorted(self.haplotype, end)

        if index_end > index_start:
            return True

        return False

    def reference_coordinate_exists_between(self, start, end):
        index_start = np.searchsorted(self.reference, start)
        index_end = np.searchsorted(self.reference, end)

        if index_end > index_start:
            return True

        return False

    def __str__(self):
        return str(self.reference) + "\n" + str(self.haplotype)


def simulate_reads(chromosome, haplotype, coverage=150, read_length=150, snv_prob=0.01, deletion_prob=0.001,
                       insertion_prob=0.001, random_seed=1, data_base_name=""):
    haplotype_fasta_file_name = "%schromosome%s_haplotype%s_reference.fasta" % (data_base_name, chromosome, haplotype)
    coordinate_map = CoordinateMap.from_file("%scoordinate_map_chromosome%s_haplotype%s.npz" % (data_base_name, chromosome, haplotype))


    np.random.seed(random_seed)
    ref = Fasta(haplotype_fasta_file_name)[chromosome]

    _simulate(ref, chromosome, coordinate_map, coverage, read_length, deletion_prob, insertion_prob, snv_prob)


def _simulate(ref, chromosome, coordinate_map, coverage, read_length=150, deletion_prob=0, insertion_prob=0, snv_prob=0,
              include_reverse_complements=False, skip_reads_with_missing_bases=True):
    chrom_length = len(ref)
    chrom_min = 0
    chrom_max = chrom_length - read_length - 10
    n_reads = int(0.5 * coverage * chrom_length / read_length)
    logging.info("Will simulate %d reads to get coverage %.3f on chromosome %s" % (n_reads, coverage, chromosome))
    logging.warning("Coverage is halved, because assuming one wants half coverage on each haplotype")
    mutator = Mutator(read_length + 10)
    logging.info("Starting simulation")
    i = 0
    while i < n_reads:

        start = randint(chrom_min, chrom_max)
        end = start + read_length + 10
        seq = str(ref[start:end])

        if np.random.randint(0, 2) == 1:
            # reverse complement
            seq = str(Seq(seq).reverse_complement())

        if skip_reads_with_missing_bases and "n" in seq or "N" in seq:
            continue

        i += 1

        ref_offset = coordinate_map.convert(start)
        mutated_seq = mutator.mutate_sequence(seq, snv_prob, deletion_prob, insertion_prob)
        mutated_seq = mutated_seq[0:read_length]

        n_nodes_not_in_linear_ref = 0
        if coordinate_map.haplotype_coordinate_exists_between(start, end):
            n_nodes_not_in_linear_ref += 1

        print("%s %d %d %s" % (chromosome, ref_offset, n_nodes_not_in_linear_ref, mutated_seq))

