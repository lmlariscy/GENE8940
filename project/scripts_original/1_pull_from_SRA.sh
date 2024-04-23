#!/bin/bash
#SBATCH --job-name=Pull_data_SRA
#SBATCH --partition=glenn_p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=50gb
#SBATCH --export=NONE
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=mandyh@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#This script pulls data from SRA using the accession number. 
#Ideal to avoid transferring large files on to local computer
# first step fetches the data and the second step splits it

module load SRA-Toolkit/3.0.3-gompi-2022a

#Set working directory
cd /file/path/to/working/directory

#define directory
output='file/path/to/Raw_data'
output2='file/path/tp/Raw_data/fastq'

##Pull data in a loop
set -ueo pipefailx
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