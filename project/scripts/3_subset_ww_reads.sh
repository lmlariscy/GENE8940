#!/bin/bash
#SBATCH --job-name=subset_reads
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=2gb
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=mandyh@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load seqtk/1.3-GCC-11.3.0

cd /work/gene8940/lml38336/project/wastewater/raw_data

#define directory
input='/work/gene8940/lml38336/project/wastewater/raw_data/fastq/combined'
output='/work/gene8940/lml38336/project/wastewater/raw_data/subset_reads'

##This is where you subset reads from the SRA ww reads
##You will need to math out the proportion of wastewater reads that you will need in your simulated mixtures
## the -s flag is the seed number - leave it as is
## the number after the fq is the number of reads that you are subsetting
## the wastewater fq files will probably have a different name than below
## you might have to combine ww fq files before you can subset. 

echo "x: $x, i: $i, z: $z"
seqtk sample -s 10 $input/ww_R1.fq 495000 > $output/ww_495000_R1.fq
seqtk sample -s 10 $input/ww_R2.fq 495000 > $output/ww_495000_R2.fq
seqtk sample -s 10 $input/ww_R1.fq 490000 > $output/ww_490000_R1.fq
seqtk sample -s 10 $input/ww_R2.fq 490000 > $output/ww_490000_R2.fq
