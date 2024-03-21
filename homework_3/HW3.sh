#!/bin/bash
#SBATCH --job-name=homework_3  		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=6		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=30gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lml38336/log.%j			# Location of standard output and error log files (replace lml38336 with your myid)
#SBATCH --mail-user=lml38336@uga.edu                    # Where to send mail (replace lml38336 with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/lml38336/homework_3"

#load modules
module load canu/2.2-GCCcore-11.2.0
module load SPAdes/3.15.5-GCC-11.3.0
module load QUAST/5.2.0-foss-2022a
module load BWA/0.7.17-GCCcore-11.3.0
module load MUMmer/4.0.0rc1-GCCcore-11.3.0

#load PacBio data
LONG="/work/gene8940/instructor_data/ecoli_p6_25x.filtered.fastq.gz"
cp $LONG $OUTDIR/ecoli_pacbio.fastq.gz 

#load Illumina data
SHORT1="/work/gene8940/instructor_data/s_6_1.fastq.gz"
cp $SHORT1 $OUTDIR/ecoli_illumina1.fastq.gz

SHORT2="/work/gene8940/instructor_data/s_6_2.fastq.gz"
cp $SHORT2 $OUTDIR/ecoli_illumina2.fastq.gz

#load MG1655 reference genome
REF="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa.gz"
curl -s $REF | gunzip -c > $OUTDIR/GCA_000005845.fna

#assemble e. coli genome (PacBio) with canu
canu -p ecoli -d $OUTDIR/canu_ecoli genomeSize=4.8m useGrid=false -pacbio-raw \ $OUTDIR/ecoli.fastq.gz

#assemble e. coli genome (Illumina) with SPAdes
spades.py -t 10 -k 21,33,55,77 --isolate --memory 24 --pe1-1 \ $OUTDIR/ecoli_illumina1.fastq.gz --pe1-2 $OUTDIR/ecoli_illumina2.fastq.gz -o $OUTDIR/spades_ecoli

#use QUAST to get quality statistics on canu assembly
quast.py -o $OUTDIR/quast_canu_ecoli -t 10 -r $OUTDIR/GCA_000005845.fna \ $OUTDIR/canu_ecoli/ecoli.contigs.fasta 

#use QUAST to get quality statistics on spades assembly
quast.py -o $OUTDIR/quast_spades_ecoli -t 10 -r $OUTDIR/GCA_000005845.fna \ $OUTDIR/spades_ecoli/scaffolds.fasta

#create mummer plot for canu
nucmer -t 10 $OUTDIR/GCA_000005845.fna $OUTDIR/canu_ecoli/ecoli.contigs.fasta -p canu_ecoli 
delta-filter -1 canu_ecoli.delta > canu_ecoli_filter.delta
mummerplot --size large -layout --color -f --png canu_ecoli_filter.delta -p \ canu_ecoli

#create mummer plot for spades
nucmer -t 10 $OUTDIR/GCA_000005845.fna $OUTDIR/spades_ecoli/contigs.fasta -p spades_ecoli 
delta-filter -1 spades_ecoli.delta > spades_ecoli_filter.delta
mummerplot --size large -layout --color -f --png spades_ecoli_filter.delta -p \ spades_ecoli