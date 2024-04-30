#!/bin/bash
#SBATCH --job-name=Kallisto_abundance_export
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=4	                        
#SBATCH --mem=24gb			                            
#SBATCH --time=2:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL  

#set directory
cd /work/gene8940/lml38336/project/kallisto

kallisto='/work/gene8940/lml38336/project/kallisto'
output='/home/lml38336/GENE8940/project/kallisto_output/abundances'

#copy abundance.tsv to github repo
for x in {1,2}
do
for y in {50,100,150}
do 
for z in {1..5}
do
scp $kallisto/M$x\_PE$y\_$z\/abundance.tsv $output/M$x\_PE$y\_$z\_abundance.tsv
done
done
done
