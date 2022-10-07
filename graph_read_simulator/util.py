import logging
import sys



def make_haplotype_fasta(chromosome, haplotype, data_dir):
    from offsetbasedgraph import IntervalCollection, Graph, SequenceGraph, Interval

    s = SequenceGraph.from_file(data_dir + "giab_chr" + chromosome + ".nobg.sequences")
    print("Getting interval")
    interval = list(
        IntervalCollection.from_file(data_dir + "haplotype_" + chromosome + "__" + haplotype + ".intervalcollection",
                                     text_file=True).intervals)[0]

    print("Getting sequence")
    sequence = s.get_interval_sequence(interval)
    print("Writing to file")
    f = open(data_dir + "giab_chr" + chromosome + "_haplotype" + haplotype + ".fasta", "w")
    f.write(">seq\n%s\n" % sequence)
    f.close()
    f.close()


def vg_path_to_obg_interval(path_file_name, out_file_name):
    from offsetbasedgraph import IntervalCollection, Graph, SequenceGraph, Interval
    from pyvg.vgobjects import Alignment
    from pyvg.conversion import get_json_lines
    from pyvg.conversion import vg_json_file_to_intervals

    json_objects = get_json_lines(path_file_name)

    alignments = (Alignment.from_json(json_object) for json_object in json_objects)
    intervals = []
    for alignment in alignments:
        path = alignment.path
        interval = path.to_obg()
        intervals.append(interval)
        chrom = path.name
        start_node = path.mappings[0].node_id()  # [m.node_id() for m in path.mappings]

        logging.info("Processing chromosome %s with start node %d" % (chrom, start_node))

        with open("chr%s_start_node.txt" % chrom, "w") as f:
            f.write(str(start_node))

        file_name = out_file_name.split(".")[0] + "_" + chrom + "." + out_file_name.split(".")[-1]
        IntervalCollection([interval]).to_file(file_name, text_file=True)
        logging.info("Number of files in interval for chrom %s: %d" % (chrom, len(interval.region_paths)))
        logging.info("Wrote path as obg interval to %s" % file_name)



def make_haplotype_paths(graph_file_name, linear_ref_path_file_name, haplotype0_file_name, haplotype1_file_name, out_base_name, chromosome):
    # Make a linear reference fasta and interval and haplotypes fasta and intervals
    from offsetbasedgraph import IntervalCollection, Graph, SequenceGraph, Interval
    from pyvg.vgobjects import Alignment
    from pyvg.conversion import get_json_lines
    from pyvg.conversion import vg_json_file_to_intervals

    chrom = chromosome
    graph = Graph.from_file(graph_file_name)
    sequence_graph = SequenceGraph.from_file(graph_file_name + ".sequences")

    linear_ref = IntervalCollection.from_file(linear_ref_path_file_name, text_file=True)
    linear_ref = list(linear_ref.intervals)[0]
    linear_ref_nodes = set(linear_ref.region_paths)

    # Write linear ref fasta to file
    linear_ref_seq = sequence_graph.get_interval_sequence(linear_ref)
    out_file = open("linear_ref_" + chrom + ".fasta", "w")
    out_file.writelines([">%s\n" % chrom])
    out_file.writelines([linear_ref_seq + "\n"])
    out_file.close()
    logging.info("Wrote linear ref sequence. N nodes in linear ref: %d" % len(linear_ref_nodes))

    haplotype_nodes = [set(), set()]  # For haplotype 0 and 1
    for haplotype in [0, 1]:
        haplotype_file_name = haplotype0_file_name
        if haplotype == 1:
            haplotype_file_name = haplotype1_file_name

        intervals = vg_json_file_to_intervals(haplotype_file_name, graph)

        for interval in intervals:
            for node in interval.region_paths:
                haplotype_nodes[haplotype].add(node)

    logging.info("N nodes in haplotype 0: %d" % len(haplotype_nodes[0]))
    logging.info("N nodes in haplotype 0 that are also in linear ref: %d" % len(haplotype_nodes[0].intersection(linear_ref_nodes)))
    logging.info("N nodes in haplotype 1: %d" % len(haplotype_nodes[1]))

    # Traverse graph to get full correct haplotype intervals
    first_nodes = graph.get_first_blocks()
    assert len(first_nodes) == 1
    logging.info("N nodes in graph: %d" % len(graph.blocks))

    for haplotype in [0, 1]:
        logging.info("Traversing haplotype %d" % haplotype)

        nodes = []
        node = first_nodes[0]
        nodes_in_haplotype = haplotype_nodes[haplotype]
        nodes_in_haplotype = set(range(0, max(linear_ref_nodes))).difference(linear_ref_nodes)
        logging.info("There are %d haplotype nodes" % len(nodes_in_haplotype))

        assert len(nodes_in_haplotype) > 0, "There are no haplotype nodes. Check that haplotype json files are not empty"

        n_haplotype_nodes = 0
        i = 0
        while True:

            nodes.append(node)
            if i % 50000 == 0:
                logging.info("#%d nodes traversed. On node %d" % (i, node))
            i += 1

            next_nodes = set(graph.adj_list[node])

            if len(next_nodes) == 0:
                logging.info("Reached end node %d with 0 edges" % node)
                break

            next_on_haplotype = next_nodes.intersection(nodes_in_haplotype)
            if len(next_on_haplotype) == 1:
                n_haplotype_nodes += 1
                next_node = list(next_on_haplotype)[0]
                assert next_node != node
                node = next_node
            elif len(next_on_haplotype) == 0:
                logging.debug("No new haplotype node from %d. Will follow reference" % node)
                # Choose reference with lowest id to avoid deletion
                node = min(list(next_nodes.intersection(linear_ref_nodes)))
            else:
                # logging.warning("There is a deletion from node %d. Choosing lowest node id as next to avoid deletion." % node)
                # This means more than one next node is on haplotype. Choose the one with lowest id to avoid taking deletion
                node = min(list(next_on_haplotype))

        logging.info("Found %d nodes. %d on haplotype" % (len(nodes), n_haplotype_nodes))
        haplotype_interval = Interval(0, graph.blocks[nodes[-1]].length(), nodes, graph)
        print("Path length: %d" % haplotype_interval.length())

        file_base_name = out_base_name + "_" + str(haplotype)
        IntervalCollection([haplotype_interval]).to_file(file_base_name + ".intervalcollection", text_file=True)

        sequence = sequence_graph.get_interval_sequence(haplotype_interval)
        out_file = open(file_base_name + ".fasta", "w")
        out_file.writelines([">%s\n" % chrom])
        out_file.writelines([sequence + "\n"])
        out_file.close()
        logging.info("Wrote fasta sequence to %s" % file_base_name + ".fasta")



