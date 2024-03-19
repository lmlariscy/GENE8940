#!/bin/bash
#SBATCH --job-name=homework_2  		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lml38336/log.%j			# Location of standard output and error log files (replace lml38336 with your myid)
#SBATCH --mail-user=lml38336@uga.edu                    # Where to send mail (replace lml38336 with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/lml38336/homework_2" 

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi   

#download genome sequence and unzip
FASTA="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa.gz"
curl -s $FASTA | gunzip -c > $OUTDIR/ecoli_MG1655.fasta

#download GFF3 annotation and unzip
GFF="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.58.chromosome.Chromosome.gff3.gz"
curl -s $GFF | gunzip -c > $OUTDIR/ecoli_MG1655.gff

#load modules
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
module load SAMtools/1.14-GCC-11.2.0
module load ucsc/443

#convert GFF3 to BED
convert2bed --input=gff < $OUTDIR/ecoli_MG1655.gff > $OUTDIR/ecoli_MG1655.bed

#filter BED file for only CDS regions
grep "ID=CDS" $OUTDIR/ecoli_MG1655.bed > $OUTDIR/ecoli_MG1655_CDS.bed

#create a genome index file for BEDtools
samtools faidx $OUTDIR/ecoli_MG1655.fasta
cut -f1,2 $OUTDIR/ecoli_MG1655.fasta.fai > $OUTDIR/ecoli_MG1655_lengths.txt

#use CDS BED file to create complementary BED intervals for non-CDS regions
bedtools complement -i $OUTDIR/ecoli_MG1655_CDS.bed \ 
-g $OUTDIR/ecoli_MG1655_lengths.txt \
> $OUTDIR/ecoli_MG1655_intergenic.bed

#generate fasta file for all CDS regions
bedtools getfasta -fi $OUTDIR/ecoli_MG1655.fasta \
-bed $OUTDIR/ecoli_MG1655_CDS.bed \
-fo $OUTDIR/ecoli_MG1655_CDS.fasta

#generate fasta file for all non-CDS regions
bedtools getfasta -fi $OUTDIR/ecoli_MG1655.fasta \
-bed $OUTDIR/ecoli_MG1655_intergenic.bed \
-fo $OUTDIR/ecoli_MG1655_intergenic.fasta

#compute GC content for CDS regions
faCount -summary $OUTDIR/ecoli_MG1655_CDS.fasta | awk '{print ($4 + $5)}' | tail -n 1 \
> $OUTDIR/ecoli_MG1655_CDS_GC.txt

#compute GC content for non-CDS regions
faCount -summary $OUTDIR/ecoli_MG1655_intergenic.fasta | awk '{print ($4 + $5)}' | tail -n 1 \
> $OUTDIR/ecoli_MG1655_intergenic_GC.txt
