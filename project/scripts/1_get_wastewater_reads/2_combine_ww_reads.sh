#!/bin/bash
#SBATCH --job-name=combine_ww_reads             
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=10		                        
#SBATCH --mem=40gb			                            
#SBATCH --time=90:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL  

#set working directory
cd /work/gene8940/lml38336/project/wastewater/raw_data/fastq

input='/work/gene8940/lml38336/project/wastewater/raw_data/fastq'
output='/work/gene8940/lml38336/project/wastewater/raw_data/fastq/unzip'
final_output='/work/gene8940/lml38336/project/wastewater/raw_data/fastq/combined'

#unzip all files
for i in {0..3}
do 
for x in {1..2}
do
gunzip -c $input/SRR1516481$i\_$x\.fastq.gz > $output/SRR1516481$i\_$x\.fastq
done
done

#combine ww files by PE read 1 and 2
for i in {0..3}
do 
for x in {1..2}
do
cat $output/SRR1516481$i\_$x\.fastq > $final_output/ww_R$x\.fq
done
done
