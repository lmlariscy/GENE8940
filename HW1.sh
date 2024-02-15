#!/bin/bash
#SBATCH --job-name=homework_1  		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lml38336/log.%j			# Location of standard output and error log files (replace lml38336 with your myid)
#SBATCH --mail-user=lml38336@uga.edu                    # Where to send mail (replace lml38336 with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/lml38336/homework_1" 

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi       

#download E. coli GFF file
URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.58.gff3.gz"      
curl -s $URL > /work/gene8940/lml38336/homework_1/ecoli_MG1655.gff

#count protein-coding sequences
grep -c "CDS" /work/gene8940/lml38336/homework_1/ecoli_MG1655.gff > /work/gene8940/lml38336/homework_1/results.txt