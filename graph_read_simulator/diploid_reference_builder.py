from pyfaidx import Fasta
import gzip
import logging
import numpy as np
from obgraph.variants import get_variant_type, VcfVariant
from .simulation import CoordinateMap


class DiploidReferenceBuilder:
    def __init__(self, reference_file_name, vcf_file_name, chromosome, haplotype=0, base_output_name=""):
        logging.info("Reading fasta %s" % reference_file_name)
        self.reference = Fasta(reference_file_name)[chromosome]
        self.chromosome = chromosome
        self.haplotype = haplotype
        self.vcf_file_name = vcf_file_name
        self._reference_coordinates = [1]
        self._haplotype_coordinates = [1]
        self._variant_positions = []
        self.base_output_name = base_output_name

    def build(self):
        haplotype_sequence = []  # Stores a list of sequences that are to be joined at the end
        current_haplotype_coordinate = 1
        previous_reference_coordinate = 1
        warned = False

        is_gzipped = False
        f = open(self.vcf_file_name)
        if self.vcf_file_name.endswith(".gz"):
            is_gzipped = True
            f = gzip.open(self.vcf_file_name)

        prev_line = None
        for line in f:
            if is_gzipped:
                line = line.decode("utf-8")

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM") or line.startswith("#chrom") or line.startswith("#Chrom"):
                sample = line.split()[9]
                logging.info("Will create two linear references for individual %s" % sample)
                continue

            #logging.info("LINE: %s" % line)
            prev_line = line
            l = line.split()
            chrom = l[0]
            if chrom != self.chromosome:
                continue

            genotype = l[9].split(":")[0]

            # We phase all genotypes that are not phased
            genotype = genotype.replace("/", "|")

            if genotype == "0|0" or genotype == ".|." or genotype == ".":
                continue

            try:
                allele = genotype.split("|")[self.haplotype]
            except IndexError:
                logging.error("Could not parse genotype %s" % genotype)
                logging.error("Original vcf line: \n%s" % line)
                raise

            if allele == "0":
                # Has the reference, do nothing
                continue

            ref_sequence = l[3]
            ref_position = int(l[1])

            vcf_variant = VcfVariant.from_vcf_line(line, dont_encode_chromosome=True)
            variant_sequence = vcf_variant.variant_sequence  #  we want variant sequence with leading base for insertions/deletions since we add reference sequence up until the variant start pos
            if "," in variant_sequence:
                variant_sequence = variant_sequence.split(",")[0]
                if not warned:
                    logging.warning("VCF has variant %s with multiple alleles. Not supported. Choosing first allele now." % vcf_variant)
                    logging.warning("There might be more such variants, will only warn for first.")
                    warned = True

            # First add all reference sequence from last time we added
            #print("Getting prev ref seq from %d to %d" % (previous_reference_coordinate, ref_position))
            if previous_reference_coordinate == ref_position:
                new_reference_sequence = ""
            elif previous_reference_coordinate > ref_position:
                logging.error("Ref position in vcf is smaller than where we are at (%d). Is this an overlapping variant? Skipping" % ref_position)
                logging.error(line)
                logging.error("Prev line: %s" % prev_line)
                continue
            elif "*" in variant_sequence:
                logging.info("SKipping variant with *")
                logging.info(line)
                continue
            else:

                try:
                    new_reference_sequence = str(self.reference[previous_reference_coordinate-1:ref_position-1])
                except ValueError:
                    logging.error("Could not get sequence from %d to %d" % (previous_reference_coordinate, ref_position))
                    raise
                except OSError:
                    logging.error("Processed line %s" % line)
                    logging.error("Could not get sequence from %d to %d" % (previous_reference_coordinate-1, ref_position-1))
                    raise

            haplotype_sequence.append(new_reference_sequence)
            #print("Adding new reference sequence: %s" % new_reference_sequence)
            #print("Adding %d to haplotype coordinate (from new ref seq)" % len(new_reference_sequence))
            current_haplotype_coordinate += len(new_reference_sequence)

            # Store ref and haplotype coordinate mappings here
            self._reference_coordinates.append(ref_position)
            #print("Adding haplotype coordinate %d" % current_haplotype_coordinate)
            self._haplotype_coordinates.append(current_haplotype_coordinate)

            # Now add new haplotype sequence for this variant and increase reference and haplotype coordinates
            previous_reference_coordinate += len(ref_sequence) + len(new_reference_sequence)
            #variant_sequence = l[4].split(",")[int(allele) - 1]  # Gets the sequence for this haplotype
            #print("Adding haplotype sequence %s" % variant_sequence)
            #print("Previous ref coordinat is now: %d" % previous_reference_coordinate)

            haplotype_sequence.append(variant_sequence)
            #print("Adding %d to haplotype coordinate (from haplotype sequence)" % len(haplotype_sequence))
            current_haplotype_coordinate += len(variant_sequence)

            #l = len(''.join(haplotype_sequence))
            #print(current_haplotype_coordinate, l)
            #assert current_haplotype_coordinate-1 == l


        logging.info("Joining haplotype sequence")
        full_haplotype_sequence = ''.join(haplotype_sequence)
        #print(full_haplotype_sequence)
        #print(self._reference_coordinates)
        #print(self._haplotype_coordinates)

        # add coordinates for end of haplotypes as well
        self._reference_coordinates.append(len(self.reference))
        self._haplotype_coordinates.append(len(full_haplotype_sequence))


        logging.info("Saving")
        coordinate_map_fil_name = "%scoordinate_map_chromosome%s_haplotype%s.npz" % (self.base_output_name, self.chromosome, self.haplotype)
        np.savez(coordinate_map_fil_name, reference=np.array(self._reference_coordinates), haplotype=np.array(self._haplotype_coordinates))
        logging.info("Saved coordinate map to %s" % coordinate_map_fil_name)
        with open("%schromosome%s_haplotype%d_reference.fasta" % (self.base_output_name, self.chromosome, self.haplotype), "w") as f:
            f.writelines([">%s\n" % self.chromosome, full_haplotype_sequence + "\n"])

        logging.info("Saved")
        return CoordinateMap(self._reference_coordinates, self._haplotype_coordinates), full_haplotype_sequence

