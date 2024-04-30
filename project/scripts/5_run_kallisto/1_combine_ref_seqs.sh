#!/bin/bash
#SBATCH --job-name=combine_ref_seqs         
#SBATCH --partition=batch		                        
#SBATCH --ntasks=1			                            
#SBATCH --cpus-per-task=4	                        
#SBATCH --mem=24gb			                            
#SBATCH --time=2:00:00  		                       
#SBATCH --output=/work/gene8940/lml38336/log.%j			
#SBATCH --mail-user=lml38336@uga.edu                    
#SBATCH --mail-type=BEGIN,END,FAIL  

#define paths to reference genomes
covid='/work/gene8940/lml38336/project/references/covid/gisaid/EPI_ISL_19070571.fasta'
flu='/work/gene8940/lml38336/project/references/flu/gisaid/gisaid_epiflu_sequences.fasta'
rsv='/work/gene8940/lml38336/project/references/rsv/gisaid/EPI_ISL_2584506.fasta'

#define output 
output='/work/gene8940/lml38336/project/references/all'

#if output directory doesn't exist, create it
if [ ! -d $output ]
then
    mkdir -p $output
fi

#combine into one file

cat $covid $flu $rsv > $output/covid_flu_rsv_refseq.fasta