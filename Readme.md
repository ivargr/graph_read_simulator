
# Graph read simulator
This is a simple collection of scripts for simulating reads for a diploid genome by using vg and genome graphs



## Running
```bash
graph_read_simulator simulate_reads -s 0.01 CHROMOSOME_HAPLOTYPE COVERAGE
```


## Complex running (using GNU parallel, e.g. on a cluster)

### Step 1: Make config files with chromosomes and haplotypes
*chromosomes.txt* (list of chromosome names, one on each line):
```
21
22
```

*haplotypes.txt* (list of chromosome names paired with haplotype numbers (either 0 or 1)):
```
21 0
21 1
22 0
22 1
```

### Step 2: Init some config variables and run the prepare data script
This will create all the graphs, intervals etc needed for simulating reads

```bash
population_vcf="population.vcf.gz"
individual_vcf="individual.vcf.gz"
linear_ref_fasta="ref.fa"
chromosomes="21,22"
n_chromosomes=2
individual_ID="HG002"

scripts/graph_read_simulator_prepare_data $population_vcf $individual_vcf $linear_ref_fasta $chromosomes $n_chromosomes $individual_ID

```

### Step 3: Simulate the reads

```bash
coverage=15
cat haplotypes.txt | parallel --line-buffer -j 2 "graph_read_simulator simulate_reads -s 0.01 {} $coverage" | graph_read_simulator assign_ids positions.tsv simulated_reads.fa
```

Note: The simulate_reads command can be run with different parameters, run `graph_read_simulator simulate_reads` to see available options.


