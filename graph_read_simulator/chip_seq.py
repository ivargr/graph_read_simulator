import logging
import numpy as np
from .simulation import CoordinateMap
from pyfaidx import Fasta
from Bio.Seq import Seq
from simple_read_mutator import Mutator
np.random.seed(1)

class ChipSeqSimulator:
    def __init__(self, chromosome, chromosome_start, chromosome_end, n_peaks, read_length=76, fragment_length=210,
                 average_fragments_per_peak=15, average_fragments_per_peak_std=5, out_file_base_name="chip_seq",
                 snv_error_prob=0.01, insertion_error_prob=0.001, deletion_error_prob=0.001, noise_coverage=0.1):

        self._mutator = Mutator(read_length + 10)

        self._haplotype_sequences = [
            Fasta("chromosome%s_haplotype%s_reference.fasta" % (chromosome, "0"))[chromosome],
            Fasta("chromosome%s_haplotype%s_reference.fasta" % (chromosome, "1"))[chromosome]
        ]

        self._coordinate_maps = [
            CoordinateMap.from_file("coordinate_map_chromosome%s_haplotype0.npz" % (chromosome)),
            CoordinateMap.from_file("coordinate_map_chromosome%s_haplotype1.npz" % (chromosome))
        ]

        self._chromosome = chromosome
        self._chromosome_start = chromosome_start
        self._chromosome_end = chromosome_end
        self._n_peaks = n_peaks
        self._fragment_length = fragment_length
        self._read_length = read_length
        self._average_fragments_per_peak = average_fragments_per_peak
        self._average_fragments_per_peak_std = average_fragments_per_peak_std
        self._binding_site_size = 10  # assumed number of base pairs for binding
        self._out_file_base_name = out_file_base_name + "_" + self._chromosome
        self._read_positions_file = None
        self._reads_file = None
        self._peaks_file = None
        self._snv_error_prob = snv_error_prob
        self._deletion_error_prob = deletion_error_prob
        self._insertion_error_prob = insertion_error_prob
        self._read_number = 0
        self._noise_coverage = noise_coverage

    def simulate_read_from_fragment(self, fragment_position_reference, haplotype, direction):
        fragment_position_haplotype = self._coordinate_maps[haplotype].convert_reverse(fragment_position_reference)

        fragment_sequence = self._haplotype_sequences[haplotype][fragment_position_haplotype:fragment_position_haplotype+self._fragment_length]
        assert len(fragment_sequence) == self._fragment_length
        effective_read_length = self._read_length + 10  # Add a bit so that we can mutate the sequence and allow it to become shorter
        if direction == 1:
            read_start = fragment_position_haplotype
            sequence = str(fragment_sequence[0:effective_read_length])
            read_position_reference = fragment_position_reference
        else:
            read_start = fragment_position_haplotype + self._fragment_length - self._read_length
            read_position_reference = fragment_position_reference + self._fragment_length - self._read_length
            sequence = str(Seq(str(fragment_sequence[-effective_read_length:])).reverse_complement())

        sequence = sequence.lower()
        if "N" in sequence or "n" in sequence:
            return

        mutated_seq = self._mutator.mutate_sequence(sequence, self._snv_error_prob, self._deletion_error_prob, self._insertion_error_prob)
        mutated_seq = mutated_seq[0:self._read_length]

        n_nodes_not_in_linear_ref = 0
        if self._coordinate_maps[haplotype].haplotype_coordinate_exists_between(read_start, read_start+self._read_length):
            n_nodes_not_in_linear_ref += 1

        self._reads_file.writelines([">%d\n%s\n" % (self._read_number, mutated_seq)])
        self._read_positions_file.writelines(["%d\t%s\t%d\t.\t.\t.\t.\t%d\n" % (self._read_number, self._chromosome, read_position_reference, n_nodes_not_in_linear_ref)])
        self._read_number += 1

    def simulate(self):
        self._read_positions_file = open(self._out_file_base_name + ".pos", "w")
        self._reads_file = open(self._out_file_base_name + ".fa", "w")
        self._peaks_file = open(self._out_file_base_name + "_peaks" + ".bed", "w")

        peak_positions = np.sort(np.random.randint(self._chromosome_start + self._fragment_length, self._chromosome_end-self._fragment_length, self._n_peaks))

        for i, peak_position in enumerate(peak_positions):
            if i % 1000 == 0:
                logging.info("%d peaks simulated" % i)

            self._peaks_file.writelines(["%s\t%d\t%d\n" % (self._chromosome, peak_position, peak_position + self._binding_site_size)])

            possible_fragment_start_start = peak_position - self._fragment_length + self._binding_site_size
            possible_fragment_start_end = peak_position

            for haplotype in [0, 1]:
                n_fragments = int(np.random.normal(self._average_fragments_per_peak, self._average_fragments_per_peak_std))
                for fragment in range(n_fragments):
                    fragment_start = np.random.randint(possible_fragment_start_start, possible_fragment_start_end)
                    direction = np.random.choice([1, -1])
                    self.simulate_read_from_fragment(fragment_start, haplotype, direction)

        # Add som noise
        n_noise_fragments = int(self._chromosome_end * self._noise_coverage / self._fragment_length)
        logging.info("Will simulate %d noise reads" % n_noise_fragments)
        positions = np.sort(np.random.randint(self._chromosome_start + self._fragment_length, self._chromosome_end-self._fragment_length, n_noise_fragments))
        for position in positions:
            direction = np.random.choice([1, -1])
            haplotype = np.random.choice([0, 1])
            self.simulate_read_from_fragment(position, haplotype, direction)

        self._reads_file.close()
        self._peaks_file.close()
        self._read_positions_file.close()
