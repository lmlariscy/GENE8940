#!/bin/bash
#SBATCH --job-name=split_sim_reads           
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=10		                        
#SBATCH --mem=40gb			                            
#SBATCH --time=90:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL  

#define directory
input='/work/gene8940/lml38336/project/references/mixed_sim_reads'
output='/work/gene8940/lml38336/project/references/mixed_split_reads'

if [ ! -d $output ]
then
    mkdir -p $output
fi

#load module
module load seqtk/1.3-GCC-11.2.0

#split paired end reads into R1 and R2
for y in {50,100,150}
do
for x in {5000,10000}
do
for z in {1..5}
do
seqtk seq -1 $input/mixed_variants_PE$y\_$x\_$z\-reads.fastq > $output/mixed_variants_PE$y\_$x\_$z\_R1.fq
seqtk seq -2 $input/mixed_variants_PE$y\_$x\_$z\-reads.fastq > $output/mixed_variants_PE$y\_$x\_$z\_R2.fq
done
done
done