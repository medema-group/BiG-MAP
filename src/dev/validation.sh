#! /bin/bash


for fastani_setting in  0.7; do #1.0 0.9 0.8
    # Creating data for geneclusters of anti- and gutSMASH
    #python3 ~/Github/Modules/2.MODULE/metaclust.genecluster.py -D /mnt/scratch/pasca005/CGR_reference_genomes/gutsmash_output/GCA_00348** -O /mnt/scratch/berg266/validation/bowtie_settings/data/$fastani_setting.fastani/ -t $fastani_setting

    #python3 ~/Github/Modules/2.MODULE/metaclust.genecluster.py -D /mnt/scratch/pasca005/CGR_reference_genomes/antismash_output/GCA_003433* -O /mnt/scratch/berg266/validation/bowtie_settings/data/$fastani_setting.fastani/ -t $fastani_setting

    

    for s in fast-local ; do  #very-fast-local sensitive-local very-sensitive-local
	# Mapping against the simulated reads
	#python3 ~/Github/Modules/3.MODULE/metaclust.map.py -R /mnt/scratch/berg266/validation/bowtie_settings/data/$fastani_setting.fastani/metaclust.GCFs_DNA_reference.fna -I1 /mnt/scratch/berg266/validation/bowtie_settings/data/grinder-reads_1.fasta -I2 /mnt/scratch/berg266/validation/bowtie_settings/data/grinder-reads_2.fasta -O /mnt/scratch/berg266/validation/bowtie_settings/$s/ -f True -s $s

	# Validating for recall and specificity
	#python3 ~/Github/Modules/VALIDATION/metaclust.validation.py -C validate -R /mnt/scratch/berg266/validation/bowtie_settings/data/grinder-reads.fa -J /mnt/scratch/berg266/validation/bowtie_settings/data/$fastani_setting.fastani/absolute_cluster_locs.json -S /mnt/scratch/berg266/validation/bowtie_settings/$s/bowtie2-map-results/grinder-reads.sam -O /mnt/scratch/berg266/validation/bowtie_settings/$s/ -F /mnt/scratch/berg266/validation/bowtie_settings/data/$fastani_setting.fastani/fastani.GCFheaders.json

	#mv /mnt/scratch/berg266/validation/bowtie_settings/$s/metrics.csv /mnt/scratch/berg266/validation/bowtie_settings/$s/$fastani_setting.$s.metrics.csv
	#mv /mnt/scratch/berg266/validation/bowtie_settings/$s/$fastani_setting.$s.metrics.csv /mnt/scratch/berg266/validation/bowtie_settings/results/
	#rm -r /mnt/scratch/berg266/validation/bowtie_settings/$s/*
    done
done







