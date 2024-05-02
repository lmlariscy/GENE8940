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

cd /home/lml38336/GENE8940/project/ref_seqs

#define paths to reference genomes
covid='covid/EPI_ISL_19070571.fasta'
covid1='covid/EPI_ISL_15358343.fasta'
covid2='covid/EPI_ISL_5713952.fasta'
flu='flu/gisaid_epiflu_sequences.fasta'
flu1='flu/gisaid_epiflu_sequence_offtarget.fasta'
rsv='rsv/EPI_ISL_2584506.fasta'
rsv1='rsv/EPI_ISL_2584486.fasta'
rsv2='rsv/EPI_ISL_18980314.fasta'

#define output 
output='/work/gene8940/lml38336/project/references/all'

#if output directory doesn't exist, create it
if [ ! -d $output ]
then
    mkdir -p $output
fi

#combine into one file

cat $covid $covid1 $covid2 $flu $flu1 $rsv $rsv1 $rsv2 > $output/covid_flu_rsv_refseq.fasta