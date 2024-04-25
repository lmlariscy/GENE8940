#!/bin/bash
#SBATCH --job-name=grinder_sim_reads_rsv	            # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=12			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=90:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lml38336/log.%j			# Location of standard output and error log files (replace lml38336 with your myid)
#SBATCH --mail-user=lml38336@uga.edu                    # Where to send mail (replace lml38336 with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                      # Mail events (BEGIN, END, FAIL, ALL)

#load modules
module load Grinder/0.5.4-GCCcore-8.3.0-Perl-5.30.0

##Will need to do this with each genome of interest

###RSV
genomes='/work/gene8940/lml38336/project/references/rsv/gisaid'
## this is the genomes that you will be using to simulate reads from
input='/work/gene8940/lml38336/project/rsv/primers'
##inputs will be fasta files with the two primers that will essentially create the PCR
output='/work/gene8940/lml38336/project/rsv/simulations'
final_output='/work/gene8940/lml38336/project/rsv/simulations/final_sim_amp_reads'


##Here we looping through the primers and generating fastqs
## i loops through primer sets
## x loops through number of reads
## y loops through the different read lengths
## z loops through number of replicates - here we make 5
## these loops are puttings these attributes into the output name so that you can call upon them in the next steps
for i in {1..6}
do
for x in {5000,10000}
do
for y in {50,100,150}
do
for z in {1..5}
do
grinder -reference_file $genomes/EPI_ISL_19057389.fasta -forward_reverse $input/rsv_$i\.fasta -total_reads $x -read_dist $y -insert_dist 350 -mate_orientation FR -mutation_dist poly4 3e-3 3.3e-10 -qual_levels 38 10 -fastq_output 1 -output_dir $output/ -base_name rsv_a$i\_PE$y_$x\_$z
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
cat $output/rsv_a$i\_PE$y_$x\_$z\-reads.fastq > $final_output/rsv_a_PE$y_$x\_$z\-reads.fastq
done
done
done