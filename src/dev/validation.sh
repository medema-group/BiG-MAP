#! /bin/bash

# Create a directory named $validation_folder with:
# - validation_folder/
# | fast/
# | fast-local/
# | ground_truth/
# | results/
# | sensitive/
# | sensitive-local/
# | validation.sh
# | very-fast/
# | very-fast-local/
# | very-sensitive/
# | very-sensitive-local/
# | - data/
#   | 0.1.fastani/
#   | 0.2.fastani/
#   | 0.3.fastani/
#   | 0.4.fastani/
#   | 0.5.fastani/
#   | 0.6.fastani/
#   | 0.7.fastani/
#   | 0.8.fastani/
#   | 0.99.fastani/
#   | 0.9.fastani/
#   | genomes.fna
#   | grinder-ranks.txt
#   | grinder-reads_1.fasta
#   | grinder-reads_2.fasta
#   | grinder-reads.fa
#   | options_grinder.txt



validation_folder='gutSMASH_PROT_EXPONENTIAL_0.05'

for fastani_setting in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99; do #0.99
    # Creating data for geneclusters of anti- and gutSMASH
    #python3 ~/Github/Modules/2.MODULE/metaclust.genecluster.py -D /mnt/scratch/pasca005/CGR_reference_genomes/gutsmash_output/GCA_00348** -O /mnt/scratch/berg266/validation/$validation_folder/data/$fastani_setting.fastani/ -t $fastani_setting

    #python3 ~/Github/Modules/2.MODULE/metaclust.genecluster.py -D /mnt/scratch/pasca005/CGR_reference_genomes/antismash_output/GCA_00348** -O /mnt/scratch/berg266/validation/$validation_folder/data/$fastani_setting.fastani/ -t $fastani_setting



    

    
    for s in very-fast ; do  # fast sensitive very-sensitive very-fast-local fast-local sensitive-local very-sensitive-local
	# Mapping against the simulated reads
	python3 ~/Github/Modules/3.MODULE/metaclust.map.py -R /mnt/scratch/berg266/validation/$validation_folder/data/$fastani_setting.fastani/metaclust.GCFs_DNA_reference.fna -I1 /mnt/scratch/berg266/validation/$validation_folder/data/grinder-reads_1.fasta -I2 /mnt/scratch/berg266/validation/$validation_folder/data/grinder-reads_2.fasta -O /mnt/scratch/berg266/validation/$validation_folder/$s/ -f True -s $s

	# Validating for recall and specificity
	python3 ~/Github/Modules/VALIDATION/metaclust.validation.py -C validate -R /mnt/scratch/berg266/validation/$validation_folder/data/grinder-reads.fa -J /mnt/scratch/berg266/validation/$validation_folder/data/$fastani_setting.fastani/absolute_cluster_locs.json -S /mnt/scratch/berg266/validation/$validation_folder/$s/bowtie2-map-results/grinder-reads.sam -O /mnt/scratch/berg266/validation/$validation_folder/$s/ -F /mnt/scratch/berg266/validation/$validation_folder/data/$fastani_setting.fastani/fastani.GCFheaders.json -g /mnt/scratch/berg266/validation/$validation_folder/ground_truth/ground_truth -s $fastani_setting

	mv /mnt/scratch/berg266/validation/$validation_folder/$s/genecluster.metrics.csv /mnt/scratch/berg266/validation/$validation_folder/results/$fastani_setting.$s.genecluster.metrics.csv
	mv /mnt/scratch/berg266/validation/$validation_folder/$s/housekeepinggene.metrics.csv /mnt/scratch/berg266/validation/$validation_folder/results/$fastani_setting.$s.housekeepinggene.metrics.csv
	
	rm -r /mnt/scratch/berg266/validation/$validation_folder/$s/*/
    done
done
