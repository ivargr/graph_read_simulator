#!/usr/bin/env bash

population_vcf="population.vcf.gz"
individual_vcf="individual.vcf.gz"
linear_ref_fasta="ref.fa"
chromosomes="21,22"
n_chromosomes=2
individual_ID="HG002"


# Make one graph per chromosome (remember --alt-paths to get variants as paths. We don't chop up nodes to smaller than 1000 bp, since we are not going to map to this graph)
cat "chromosomes.txt" | parallel -j $n_chromosomes "time vg construct -C -R {} -r $linear_ref_fasta -v $individual_vcf -t 1 -m 1000 --alt-paths > individual_chr{}.vg"

# Node id conversion
vg ids -j $(for chromosome in $(echo $chromosomes | tr "," "\n"); do echo individual_chr$chromosome.vg; done)

# Index the graph (make a gbwt index that will include all phased variants as paths. We will use this gbwt index later to pull out paths containing phased variants)
vg index -x individual.xg -G individual.gbwt -v individual.vcf.gz $(for chromosome in $(echo $chromosomes | tr "," "\n"); do echo individual_chr$chromosome.vg; done)


# Get all haplotype paths in gam for each chromosome and each haplotype
# Also, get haplotype paths in vg format
# (need the vg paths for later removing all nodes not part of paths to make graphs for single haplotypes, which will be used for read simulation)
for chromosome in $(echo $chromosomes | tr "," "\n")
    do
    echo "Chromosome $chromosome"
	vg paths --gbwt individual.gbwt --extract-vg -x individual.xg -Q _thread_${individual_ID}_${chromosome}_0 > haplotype_${chromosome}_0.paths &
	vg paths --gbwt individual.gbwt --extract-vg -x individual.xg -Q _thread_${individual_ID}_${chromosome}_1 > haplotype_${chromosome}_1.paths &

	vg paths --gbwt individual.gbwt --extract-gam -x individual.xg -Q _thread_${individual_ID}_${chromosome}_0 > haplotype_${chromosome}_0.gam &
	vg paths --gbwt individual.gbwt --extract-gam -x individual.xg -Q _thread_${individual_ID}_${chromosome}_1 > haplotype_${chromosome}_1.gam &
done
wait

# Extract a linear reference path through each chromosome
vg paths -X -x individual.xg > individual_reference_paths.gam
vg view -aj individual_reference_paths.gam > individual_reference_paths.json
graph_read_simulator vg_path_to_obg_interval individual_reference_paths.json individual_reference_path.intervalcollection


# Convert haplotype paths and vg graphs to json
for chromosome in $(echo $chromosomes | tr "," "\n")
    do
    vg view -aj haplotype_${chromosome}_0.gam > haplotype_${chromosome}_0.json &
    vg view -aj haplotype_${chromosome}_1.gam > haplotype_${chromosome}_1.json &
    vg view -Vj individual_chr$chromosome.vg > individual_chr$chromosome.json &
done
wait

# Create ob graph for each chromosome
for chromosome in $(echo $chromosomes | tr "," "\n")
    do
    graph_peak_caller create_ob_graph individual_chr$chromosome.json &
done
wait

# Note: The haplotype paths are not complete. They may have gaps in them
# Thus, we want to traverse the giab graph, folloing the haplotype paths, but fill in with reference paths where possible
for chromosome in $(echo $chromosomes | tr "," "\n")
    do
	graph_read_simulator make_haplotype_paths individual_chr$chromosome.nobg individual_reference_path_$chromosome.intervalcollection haplotype_${chromosome}_0.json haplotype_${chromosome}_1.json haplotype_${chromosome}_ $chromosome &
	#python3 ../make_linear_reference_and_interval.py giab_chr$chromosome.nobg giab_reference_path_$chromosome.intervalcollection haplotype_${chromosome}_0.json haplotype_${chromosome}_1.json haplotype_${chromosome}_ $chromosome &
done
wait

cat haplotype_*__0.fasta >> haplotype0.fasta
cat haplotype_*__1.fasta >> haplotype1.fasta
# Verify by grep ">" haplotype0.fasta


# Index reference and haplotype intervals (we need them later to convert coordinates from haplotype interval offset to linear ref offset)
for chromosome in $(echo $chromosomes | tr "," "\n")
    do
	graph_peak_caller index_interval -g individual_chr$chromosome.nobg individual_reference_path_$chromosome.intervalcollection &
	graph_peak_caller index_interval -g individual_chr$chromosome.nobg haplotype_${chromosome}__0.intervalcollection &
	graph_peak_caller index_interval -g individual_chr$chromosome.nobg haplotype_${chromosome}__1.intervalcollection &
done


# Make modified graphs only containing the haplotypes, will be used to simulate reads from
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"
for chromosome in $(echo $chromosomes | tr "," "\n")
	do
	vg mod -D giab_chr$chromosome.vg > giab_chr${chromosome}_haplotype0.tmp && \
	cat haplotype_${chromosome}_0.paths >> giab_chr${chromosome}_haplotype0.tmp && \
	vg mod -N giab_chr${chromosome}_haplotype0.tmp | vg mod -D - > giab_chr${chromosome}_haplotype0.vg && \
	vg index -x giab_chr${chromosome}_haplotype0.xg giab_chr${chromosome}_haplotype0.vg &

	vg mod -D giab_chr$chromosome.vg > giab_chr${chromosome}_haplotype1.tmp && \
	cat haplotype_${chromosome}_1.paths >> giab_chr${chromosome}_haplotype1.tmp && \
	vg mod -N giab_chr${chromosome}_haplotype1.tmp | vg mod -D - > giab_chr${chromosome}_haplotype1.vg && \
	vg index -x giab_chr${chromosome}_haplotype1.xg giab_chr${chromosome}_haplotype1.vg &
done


# Create a whole genome graph with only reference, used by vg to annotate positions on
cat giab_chr?.vg giab_chr??.vg > giab_all_chromosomes.vg
# Extract reference paths in vg format for each chromosome (bit hackish since vg only allows fetching by prefix, not by actual name)
> reference_paths.vg
chromosomes="1,2,3,4,5,6,7,8,9,X"
for chromosome in $(echo $chromosomes | tr "," "\n")
	do
	vg paths -Q $chromosome -V -x giab.xg >> reference_paths.vg
done

# Verify by vg paths -L -v reference_paths.vg
# Remove all paths from giab_all_chromosomes.vg and add these paths:
vg mod -D giab_all_chromosomes.vg > tmp.vg
cat reference_paths.vg >> tmp.vg
vg mod -N tmp.vg > giab_reference.vg
vg index -x giab_reference.xg giab_reference.vg


# Simulate by running
mkdir simulation_dir
nohup ./vg_benchmark.sh results_3m hg19_chr1-Y.fa None /data/bioinf/whole_genome/wg1.6 data/giab_chr20_haplotype0 data/giab_chr20_haplotype1 data/giab data/giab_reference 75 "--forward-only -n 1500000 -e 0.01 -i 0.002 -l 150" 2358792 150 "" > log.txt &

