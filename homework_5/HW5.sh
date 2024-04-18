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
OUTDIR="/work/gene8940/lml38336/homework_5"
KALLISTO="/work/gene8940/lml38336/homework_5/kallisto"               

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

if [ ! -d $KALLISTO ]
then
    mkdir -p $KALLISTO
fi

#load modules
module load kallisto/0.48.0-gompi-2022a

#download E. coli CDS fasta file
CDS="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
curl -s $CDS | gunzip -c > $OUTDIR/GCA_000005845_CDS.fa

#create kallisto index file from CDS file
kallisto index -i $OUTDIR/ecoli_MG1655_cds.fa.idx $OUTDIR/ecoli_MG1655_cds.fa

#pseudoalign RNA-seq reads and get transcript abundance estimations
for i in SRR5344681 SRR5344682 SRR5344683 SRR5344684
do
  kallisto quant -t 10 -b 100 -i $OUTDIR/ecoli_MG1655_cds.fa.idx -o $KALLISTO/$i /work/gene8940/instructor_data/${i}_1.fastq.gz /work/gene8940/instructor_data/${i}_2.fastq.gz
done

#activate conda environment and execute R script
source activate sleuth
R --no-save < /home/lml38336/GENE8940/homework_5/homework5.r
