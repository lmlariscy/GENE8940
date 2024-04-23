#!/bin/bash
#SBATCH --job-name=Combining_mix
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=1gb
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=mandyh@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#Set directory
cd /scratch/mandyh/WISER_InSilico/new_variant_detect/final_sim_amp_reads

#define directory
output='/file/path/to/mix_sim_reads'


##Here you combine the files that you made from single genomes
## Looping through read length ($y) and reps($z)

##This is also where you can add in wastewater reads.
##You will need to math out proportions by number of reads

for y in {50,100,150}
do
for z in {1..5}
do
cat Alpha_V4.1_PE$y\_2500_$z\-reads.fastq Beta_V4.1_PE$y\_2500_$z\-reads.fastq Delta_V4.1_PE$y\_2500_$z\-reads.fastq wh1_V4.1_PE$y\_2500_$z\-reads.fastq OmBA_1_V4.1_PE$y\_475000_$z\-reads.fastq OmBA_1_1_V4.1_PE$y\_5000_$z\-reads.fastq > $output/M1_V4.1_PE$y\_$z\-reads.fastq
cat Alpha_V4.1_PE$y\_2500_$z\-reads.fastq Beta_V4.1_PE$y\_2500_$z\-reads.fastq Delta_V4.1_PE$y\_2500_$z\-reads.fastq wh1_V4.1_PE$y\_2500_$z\-reads.fastq OmBA_1_V4.1_PE$y\_470000_$z\-reads.fastq OmBA_1_1_V4.1_PE$y\_10000_$z\-reads.fastq > $output/M2_V4.1_PE$y\_$z\-reads.fastq
cat Alpha_V4.1_PE$y\_2500_$z\-reads.fastq Beta_V4.1_PE$y\_2500_$z\-reads.fastq Delta_V4.1_PE$y\_2500_$z\-reads.fastq wh1_V4.1_PE$y\_2500_$z\-reads.fastq OmBA_1_V4.1_PE$y\_455000_$z\-reads.fastq OmBA_1_1_V4.1_PE$y\_25000_$z\-reads.fastq > $output/M3_V4.1_PE$y\_$z\-reads.fastq
cat Alpha_V4.1_PE$y\_2500_$z\-reads.fastq Beta_V4.1_PE$y\_2500_$z\-reads.fastq Delta_V4.1_PE$y\_2500_$z\-reads.fastq wh1_V4.1_PE$y\_2500_$z\-reads.fastq OmBA_1_V4.1_PE$y\_355000_$z\-reads.fastq OmBA_1_1_V4.1_PE$y\_125000_$z\-reads.fastq >$output/M4_V4.1_PE$y\_$z\-reads.fastq
cat Alpha_V4.1_PE$y\_2500_$z\-reads.fastq Beta_V4.1_PE$y\_2500_$z\-reads.fastq Delta_V4.1_PE$y\_2500_$z\-reads.fastq wh1_V4.1_PE$y\_2500_$z\-reads.fastq OmBA_1_V4.1_PE$y\_230000_$z\-reads.fastq OmBA_1_1_V4.1_PE$y\_250000_$z\-reads.fastq > $output/M5_V4.1_PE$y\_$z\-reads.fastq
cat Alpha_V4.1_PE$y\_2500_$z\-reads.fastq Beta_V4.1_PE$y\_2500_$z\-reads.fastq Delta_V4.1_PE$y\_2500_$z\-reads.fastq wh1_V4.1_PE$y\_2500_$z\-reads.fastq OmBA_1_V4.1_PE$y\_105000_$z\-reads.fastq OmBA_1_1_V4.1_PE$y\_375000_$z\-reads.fastq >$output/M6_V4.1_PE$y\_$z\-reads.fastq
done
done