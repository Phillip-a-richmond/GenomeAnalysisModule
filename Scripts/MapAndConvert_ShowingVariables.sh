# Load Conda
source  /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/miniconda3/etc/profile.d/conda.sh
conda  activate /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/AnnotateVariantsEnvironment

# Arrange your data
cp /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/Files/Sample1*fastq /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Sherlock/

# bwa command
# bwa mem  <genomeFastaIndex>  <R1.fastq>  <R2.fastq>  >   <out.sam>
bwa mem /mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa Sample1_R1.fastq Sample1_R2.fastq  >  Sample1.sam

# samtools view, sort, index
# samtools view  -b  <in.sam>   -o  <out.bam>
samtools view -b  Sample1.sam   -o  Sample1.bam

# samtools sort  <in.bam>  -o  <out.sorted.bam>
samtools sort Sample1.bam  -o  Sample1.sorted.bam

# samtools index  <in.sorted.bam>
samtools   index   Sample1.sorted.bam


## Redone with variables instead
# Variables are first defined (no spaces allowed!)
WORKING_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Sherlock/
SAMPLE=Sample1
GENOME_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
FASTQR1=${SAMPLE}_R1.fastq
FASTQR2=${SAMPLE}_R2.fastq

bwa mem $GENOME_INDEX $FASTQR1 $FASTQR2 > ${SAMPLE}.sam

samtools view -b ${SAMPLE}.sam  -o ${SAMPLE}.bam

samtools sort ${SAMPLE}.bam -o ${SAMPLE}.sorted.bam

samtools index ${SAMPLE}.sorted.bam


## Variables and line-continuations
# For BWA mem, I'm adding a read group tag. this will help us later

bwa mem $GENOME_INDEX \
	$FASTQR1 \
	$FASTQR2 \
	-R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" \
	> ${SAMPLE}.sam

samtools view \
	-b ${SAMPLE}.sam  \
	-o ${SAMPLE}.bam

samtools sort \
	${SAMPLE}.bam \
	-o ${SAMPLE}.sorted.bam

samtools index \
	${SAMPLE}.sorted.bam


