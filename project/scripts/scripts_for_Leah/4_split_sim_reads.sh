#!/bin/bash
#SBATCH --job-name=Splitting_reads
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=5gb
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=mandyh@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load seqtk/1.2-foss-2019b

cd /scratch/mandyh/WISER_InSilico/new_variant_detect

#define directory
input='/file/path/to/mix_sim_reads'
output='/file/path/to/Split_mix_reads'


for i in {1..6}
do
for y in {50,100,150}
do
for z in {1..5}
do
echo "x: $x, i: $i, z: $z"
seqtk seq -1 $input/M$i\_V4.1_PE$y\_$z\-reads.fastq > $output/M$i\_V4.1_PE$y\_$z\_R1.fq
seqtk seq -2 $input/M$i\_V4.1_PE$y\_$z\-reads.fastq > $output/M$i\_V4.1_PE$y\_$z\_R2.fq
done
done
done