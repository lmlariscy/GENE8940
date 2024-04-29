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
sim_reads='/work/gene8940/lml38336/project/references/mixed_split_reads'
ww_reads='/work/gene8940/lml38336/project/wastewater/raw_data/subset_reads'
output='/work/gene8940/lml38336/project/spiked_reads'

#if output directory doesn't exist, create it
if [ ! -d $output ]
then
    mkdir -p $output
fi

#combine simulated reads and wastewater reads
for y in {50,100,150}
do
for z in {1..5}
do 
cat $sim_reads/mixed_variants_PE$y\_5000_$z\_R1.fq $ww_reads/ww_485000_R1.fq > $output/M1_PE$y\_$z\_R1.fq
cat $sim_reads/mixed_variants_PE$y\_5000_$z\_R2.fq $ww_reads/ww_485000_R2.fq > $output/M1_PE$y\_$z\_R2.fq
cat $sim_reads/mixed_variants_PE$y\_10000_$z\_R1.fq $ww_reads/ww_470000_R1.fq > $output/M2_PE$y\_$z\_R1.fq
cat $sim_reads/mixed_variants_PE$y\_10000_$z\_R2.fq $ww_reads/ww_470000_R2.fq > $output/M2_PE$y\_$z\_R2.fq
done
done

