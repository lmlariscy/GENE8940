#!/bin/bash
#SBATCH --job-name=BLAST-test		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=6	                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=30gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lml38336/log.%j			# Location of standard output and error log files (replace lml38336 with your myid)
#SBATCH --mail-user=lml38336@uga.edu                    # Where to send mail (replace lml38336 with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/lml38336/project"    
REF_DIR="/work/gene8940/lml38336/project/references"             

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# Set GISAID username and password
USERNAME="lmlariscy"
PASSWORD="Cowsgo_teehee1"

# URL to GISAID's data search API
GISAID_API="https://www.gisaid.org/epiflu-applications/sequenceDownloadsApi"

# Specify the accession number you want to download
ACCESSION_NUMBER="EPI_ISL_1234567"  # Replace with the actual accession number you want

# Build the curl command with authentication and search parameters
curl -u "$USERNAME:$PASSWORD" \
     -X POST \
     -d "accessionId=$ACCESSION_NUMBER" \
     -d "fileFormat=fasta" \
     "$GISAID_API"

for i in EPI_ISL_18914328 
do
curl -u "$USERNAME:$PASSWORD" -X POST -d "accessionId=$i" -d "fileFormat=fasta" "$GISAID_API" > $REF_DIR/${i}.fasta
done