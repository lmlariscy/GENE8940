#!/bin/bash
#SBATCH --job-name=Kallisto_quant_psuedobam
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=4	                        
#SBATCH --mem=40gb			                            
#SBATCH --time=4:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL  

##kallisto will make output folders for each sample containing these files: abundance.h5, abundance.tsv, run_info.json and pseudoalignments.bam

module load kallisto/0.48.0-gompi-2022a

#Set directory
cd /work/gene8940/lml38336/project/kallisto

#define directory
index='/work/gene8940/lml38336/project/references/all'
input='/work/gene8940/lml38336/project/spiked_reads'
output='/work/gene8940/lml38336/project/kallisto'

if [ ! -d $output ]
then
    mkdir -p $output
fi

#for loop through mixtures, read lengths, and replicates
for x in {1,2}
do
for y in {50,100,150}
do
for z in {1..5}
do
kallisto quant -t 6 -b 100 --pseudobam -i $index/sequences.kallisto_idx -o $output/M$x\_PE$y\_$z $input/M$x\_PE$y\_$z\_R1.fq $input/M$x\_PE$y\_$z\_R2.fq
done
done
done
