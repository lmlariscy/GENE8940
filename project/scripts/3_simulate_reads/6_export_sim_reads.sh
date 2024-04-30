#!/bin/bash
#SBATCH --job-name=export_sim_reads           
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=6                       
#SBATCH --mem=24gb			                            
#SBATCH --time=4:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL  

#define directories
covid='/work/gene8940/lml38336/project/references/covid/final_sim_amp_reads'
flu='/work/gene8940/lml38336/project/references/flu/final_sim_amp_reads'
rsv='/work/gene8940/lml38336/project/references/rsv/final_sim_amp_reads'

output='/home/lml38336/GENE8940/project/simulations'

for y in {50,100,150}
do
scp $covid/covid_PE$y\_5000_1-reads.fastq $output/covid_PE$y\_5000-reads.fastq
scp $flu/flu_PE$y\_5000_1-reads.fastq $output/flu_PE$y\_5000-reads.fastq
scp $rsv/rsv_PE$y\_5000_1-reads.fastq $output/rsv_PE$y\_5000-reads.fastq
done

