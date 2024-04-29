#!/bin/bash
#SBATCH --job-name=combine_sim_reads           
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=10		                        
#SBATCH --mem=40gb			                            
#SBATCH --time=90:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL  

#define directories
covid='/work/gene8940/lml38336/project/references/covid/final_sim_amp_reads'
flu='/work/gene8940/lml38336/project/references/flu/final_sim_amp_reads'
rsv='/work/gene8940/lml38336/project/references/rsv/final_sim_amp_reads'

output='/work/gene8940/lml38336/project/references/mixed_sim_reads'

#if output directory doesn't exist, create it
if [ ! -d $output ]
then
    mkdir -p $output
fi

##Here you combine the files that you made from single genomes
## Looping through read length ($y) and reps($z)

##This is also where you can add in wastewater reads.
##You will need to math out proportions by number of reads

for y in {50,100,150}
do
for x in {5000,10000}
do
for z in {1..5}
do
cat $covid/covid_PE$y\_$x\_$z\-reads.fastq $flu/flu_PE$y\_$x\_$z\-reads.fastq $rsv/rsv_PE$y\_$x\_$z\-reads.fastq > $output/mixed_variants_PE$y\_$x\_$z\-reads.fastq
done
done
done 