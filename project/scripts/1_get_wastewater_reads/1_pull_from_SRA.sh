#!/bin/bash
#SBATCH --job-name=pull_ww_sra		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lml38336/log.%j			# Location of standard output and error log files (replace lml38336 with your myid)
#SBATCH --mail-user=lml38336@uga.edu                    # Where to send mail (replace lml38336 with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                      # Mail events (BEGIN, END, FAIL, ALL)

#This script pulls data from SRA using the accession number. 
#Ideal to avoid transferring large files on to local computer
# first step fetches the data and the second step splits it

module load SRA-Toolkit/3.0.3-gompi-2022a

#define directory
output='/work/gene8940/lml38336/project/wastewater/raw_data'
output2='/work/gene8940/lml38336/project/wastewater/raw_data/fastq'

#if output directory doesn't exist, create it
if [ ! -d $output ]
then
    mkdir -p $output
fi

#if output directory doesn't exist, create it
if [ ! -d $output2 ]
then
    mkdir -p $output2
fi

##Pull data in a loop
set -ueo pipefail
SAMPLES="SRR15164810
SRR15164811
SRR15164812
SRR15164813
"
i=1

for i in $SAMPLES
do
prefetch -O $output $i
fastq-dump --split-files --gzip $output/$i\/$i\.sra -O $output2
done