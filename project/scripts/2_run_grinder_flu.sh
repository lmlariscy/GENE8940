#!/bin/bash
#SBATCH --job-name=simulate_reads_flu	               
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=4		                        
#SBATCH --mem=40gb			                            
#SBATCH --time=2:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL                      

#set working directory
cd /work/gene8940/lml38336/project/references/flu

#load modules
module load Grinder/0.5.4-foss-2022a

##Will need to do this with each genome of interest

###RSV
genomes='/work/gene8940/lml38336/project/references/flu/gisaid'
## this is the genomes that you will be using to simulate reads from
input='/work/gene8940/lml38336/project/references/flu/primers'
##inputs will be fasta files with the two primers that will essentially create the PCR
output='/work/gene8940/lml38336/project/references/flu/simulations'
final_output='/work/gene8940/lml38336/project/references/flu/final_sim_amp_reads'


##Here we looping through the primers and generating fastqs
## i loops through primer sets
## x loops through number of reads
## y loops through the different read lengths
## z loops through number of replicates - here we make 5
## these loops are puttings these attributes into the output name so that you can call upon them in the next steps
for i in {1..2}
do
for x in {5000,10000}
do
for y in {50,100,150}
do
for z in {1..5}
do
grinder -reference_file $genomes/gisaid_epiflu_sequences.fasta -forward_reverse $input/flu_$i\.fas -total_reads $x -read_dist $y -insert_dist 350 -mate_orientation FR -mutation_dist poly4 3e-3 3.3e-10 -qual_levels 38 10 -fastq_output 1 -output_dir $output/ -base_name flu_$i\_PE$y\_$x\_$z
done
done
done
done

##here we combine the fastq files from different primer positions into one fastq file by number of reads, read length and replicate
for i in {1..2}
do
for x in {5000,10000}
do
for y in {50,100,150}
do
for z in {1..5}
do
cat $output/flu_$i\_PE$y\_$x\_$z\-reads.fastq > $final_output/flu_PE$y\_$x\_$z\-reads.fastq
done
done
done
done