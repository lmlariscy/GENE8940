#!/bin/bash
#SBATCH --job-name=grinder_sim_reads
#SBATCH --partition=glenn_p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=40gb
#SBATCH --export=NONE
#SBATCH --time=90:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=mandyh@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load Grinder/0.5.4-GCCcore-8.3.0-Perl-5.30.0

cd /file/path/to/working/directory

##Will need to do this with each genome of interest

genomes='/file/path/to/ref/genomes'
## this is the genomes that you will be using to simulate reads from
input='/file/path/to/input/primers'
##inputs will be fasta files with the two primers that will essentially create the PCR
output='/file/path/to/sim_amp_outputs'
final_output='/scratch/mandyh/WISER_InSilico/new_variant_detect/final_sim_amp_reads'

##Here we looping through the primers and generating fastqs
## i loops through primer sets
## x loops through number of reads
## y loops through the different read lengths
## z loops through number of replicates - here we make 5
## these loops are puttings these attributes into the output name so that you can call upon them in the next steps
for i in {2..113}
do
for x in {5000,10000}
do
for y in {50,100,150}
do
for z in {1..5}
do
grinder -reference_file $genomes/Omicron_BA_1_1_gisaid_hcov-19_2023_04_04_18.fasta -forward_reverse $input/V4.1_POS$i\.fas -total_reads $x -read_dist $y -insert_dist 350 -mate_orientation FR -mutation_dist poly4 3e-3 3.3e-10 -qual_levels 38 10 -fastq_output 1 -output_dir $output/ -base_name OmBA_1_1_POS$i\_V4.1_PE$y_$x\_$z
done
done
done
done

##here we combine the fastq files from different primer positions into one fastq file by number of reads, read length and replicate
for x in {5000,10000}
do
for y in {50,100,150}
do
for z in {1..5}
do
cat $output/OmBA_1_1_POS$i\_V4.1_PE$y_$x\_$z\-reads.fastq > $final_output/OmBA_1_1_V4.1_PE$y_$x\_$z\-reads.fastq
done
done
done