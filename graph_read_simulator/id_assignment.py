import sys
import logging


def assign_ids(truth_file_name, fasta_file_name):

    truth_file = open(truth_file_name, "w")
    fasta_file = open(fasta_file_name, "w")


    for i, line in enumerate(sys.stdin):
        if i % 100000 == 0:
            logging.info("%d processed" % i)

        data = line.split()
        identifyer = "%09d" % i
        chromosome = data[0]
        ref_offset = data[1]
        n_nodes_not_in_linear_ref = data[2]
        sequence = data[3]

        fasta_file.writelines([">%s\n" % identifyer, "%s\n" % sequence])

        truth_file.writelines(["%s %s %s 150 0 0 150 %s %s 0\n" % (identifyer, chromosome, ref_offset, n_nodes_not_in_linear_ref, n_nodes_not_in_linear_ref)])

    truth_file.close()
    fasta_file.close()

    logging.info("Done")
