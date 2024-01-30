#!/bin/bash
#SBATCH --job-name=BLAST-test		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lml38336/log.%j			# Location of standard output and error log files (replace lml38336 with your myid)
#SBATCH --mail-user=lml38336@uga.edu                    # Where to send mail (replace lml38336 with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/lml38336/blast"                  # replace lml38336 in the following line with your myid
QUERY="/home/lml38336/GENE8940/sample.fasta"            # replace lml38336 in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load blast module
module load BLAST+/2.13.0-gompi-2022a

# run blast against local copy of NCBI nucleotide database
#blastn -num_threads 2 -query $QUERY -db /db/ncbiblast/nt/06042020/nt -out $OUTDIR/sample.fa.blastn.${SLURM_JOB_ID}.tsv -outfmt 6 -max_target_seqs 2

# same command as above, but split over multiple lines for improved readability
blastn -num_threads 4 \
       -query $QUERY \
       -db /db/ncbiblast/nt/06042020/nt \
       -out $OUTDIR/sample.fa.blastn.${SLURM_JOB_ID}.tsv \
       -outfmt 0 \
       -max_target_seqs 2
