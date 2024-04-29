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
input='/work/gene8940/lml38336/project/wastewater/raw_data/fastq/combined'
output='/work/gene8940/lml38336/project/wastewater/raw_data/subset_reads'

#load module
module load seqtk/1.3-GCC-11.2.0

##This is where you subset reads from the SRA ww reads
##You will need to math out the proportion of wastewater reads that you will need in your simulated mixtures
## the -s flag is the seed number - leave it as is
## the number after the fq is the number of reads that you are subsetting
## the wastewater fq files will probably have a different name than below
## you might have to combine ww fq files before you can subset. 

seqtk sample -s 10 $input/ww_R1.fq 485000 > $output/ww_485000_R1.fq
seqtk sample -s 10 $input/ww_R2.fq 485000 > $output/ww_485000_R2.fq
seqtk sample -s 10 $input/ww_R1.fq 470000 > $output/ww_470000_R1.fq
seqtk sample -s 10 $input/ww_R2.fq 470000 > $output/ww_470000_R2.fq
