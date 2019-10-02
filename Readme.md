
# Graph read simulator
This is a simple collection of scripts for simulating reads for a diploid genome by using vg and genome graphs

# Prepare graphs and vcfs

# Simulate the reads

```bash
coverage=15
cat haplotypes.txt | parallel --line-buffer -j 2 "graph_read_simulator simulate_reads {} $coverage" | graph_read_simulator assign_ids positions.tsv simulated_reads.fa

```

