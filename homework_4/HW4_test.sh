
#set output directory variable
OUTDIR="/work/gene8940/lml38336/homework_4"                 

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#load modules
module load SRA-Toolkit/3.0.1-centos_linux64
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load BCFtools/1.15.1-GCC-11.3.0

#download ecoli illumina reads from NCBI
prefetch -O $OUTDIR/ SRR8082143
prefetch -O /work/gene8940/lml38336/homework_4/ SRR8082143

#validate prefetched NCBI data
vbd-validate $OUTDIR/ SRR8082143
vdb-validate /work/gene8940/lml38336/ SRR8082143

#extract all NCBI paired-end reads
fastq-dump --split-files --gzip $OUTDIR/SRR8082143 -O $OUTDIR/SRR8082143
fastq-dump --split-files --gzip /work/gene8940/lml38336/homework_4/SRR8082143 -O /work/gene8940/lml38336/homework_4/SRR8082143

#download reference genome from ensemble
FASTA="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa.gz"
mkdir $OUTDIR/refseq
curl -s $FASTA | gunzip -c > $OUTDIR/refseq/ecoli_MG1655.fasta

#construct BWA index file for reference genome
bwa index refseq/ecoli_MG1655.fasta
bwa index $OUTDIR/refseq/ecoli_MG1655.fasta

#map illumina reads to reference genome (create .sam file)
bwa mem -t 6 $OUTDIR/refseq/ecoli_MG1655.fasta $OUTDIR/SRR8082143/SRR8082143_1.fastq.gz $OUTDIR/SRR8082143/SRR8082143_2.fastq.gz \
 > $OUTDIR/SRR8082143/SRR8082143.sam
bwa mem -t 6 refseq/ecoli_MG1655.fasta SRR8082143/SRR8082143_1.fastq.gz SRR8082143/SRR8082143_2.fastq.gz > SRR8082143/SRR8082143.sam

#convert .sam to .bam
samtools view $OUTDIR/SRR8082143/SRR8082143.sam -O BAM -o $OUTDIR/SRR8082143/SRR8082143.bam
samtools view SRR8082143/SRR8082143.sam -O BAM -o SRR8082143/SRR8082143.bam

#sort .bam file 
samtools sort --threads 6 $OUTDIR/SRR8082143/SRR8082143.bam -o $OUTDIR/SRR8082143/SRR8082143.sorted.bam
samtools sort --threads 6 SRR8082143/SRR8082143.bam -o SRR8082143/SRR8082143.sorted.bam

#create index for sorted .bam file
samtools index -@ 6 $OUTDIR/SRR8082143/SRR8082143.sorted.bam
samtools index -@ 6 SRR8082143/SRR8082143.sorted.bam

#compute genotype likelihoods and create mpileup
bcftools mpileup -Oz --threads 6 --min-MQ 60 -f $OUTDIR/refseq/ecoli_MG1655.fasta $OUTDIR/SRR8082143/SRR8082143.sorted.bam > $OUTDIR/SRR8082143/SRR8082143.sorted.mpileup.vcf.gz
bcftools mpileup -Oz --threads 6 --min-MQ 60 -f refseq/ecoli_MG1655.fasta SRR8082143/SRR8082143.sorted.bam > SRR8082143/SRR8082143.sorted.mpileup.vcf.gz

#call variants from mpileup
bcftools call -Oz -m -v --threads 6 --ploidy 1 $OUTDIR/SRR8082143/SRR8082143.sorted.mpileup.vcf.gz > $OUTDIR/SRR8082143/SRR8082143.sorted.mpileup.call.vcf.gz
bcftools call -Oz -m -v --threads 6 --ploidy 1 SRR8082143/SRR8082143.sorted.mpileup.vcf.gz > SRR8082143/SRR8082143.sorted.mpileup.call.vcf.gz

#filter variants calls (remove variants with quality score less than 40 and less than 10 mapped reads)
bcftools filter -Oz -e 'QUAL<40 || DP<10' $OUTDIR/SRR8082143/SRR8082143.sorted.mpileup.call.vcf.gz > $OUTDIR/SRR8082143/SRR8082143.sorted.mpileup.call.filter.vcf.gz
bcftools filter -Oz -e 'QUAL<40 || DP<10' SRR8082143/SRR8082143.sorted.mpileup.call.vcf.gz > SRR8082143/SRR8082143.sorted.mpileup.call.filter.vcf.gz

#index gzipped filtered variant calls (for viewing purposes)
bcftools index $OUTDIR/SRR8082143/SRR8082143.sorted.mpileup.call.filter.vcf.gz
bcftools index SRR8082143/SRR8082143.sorted.mpileup.call.filter.vcf.gz