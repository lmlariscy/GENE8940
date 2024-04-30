#!/bin/bash
#SBATCH --job-name=Kallisto_indexing
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=10		                        
#SBATCH --mem=40gb			                            
#SBATCH --time=2:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL 

#load module
module load kallisto/0.48.0-gompi-2022a

#Set directory
cd /work/gene8940/lml38336/project/references/all
ref_seq='covid_flu_rsv_refseq.fasta'


kallisto index -i sequences.kallisto_idx $ref_seq