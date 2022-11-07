# Load Conda
source  /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda  activate  GenomeAnalysis


### Version 1, same as we did in class #3
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


### Version 2, Redone with variables instead
# Variables are first defined (no spaces allowed!)
WORKING_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Sherlock/
SAMPLE=Sample1
GENOME_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
FASTQ_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/Files/
FASTQR1=${FASTQ_DIR}/${SAMPLE}_R1.fastq
FASTQR2=${FASTQ_DIR}/${SAMPLE}_R2.fastq

bwa mem $GENOME_INDEX $FASTQR1 $FASTQR2 > ${SAMPLE}.sam

samtools view -b ${SAMPLE}.sam  -o ${SAMPLE}.bam

samtools sort ${SAMPLE}.bam -o ${SAMPLE}.sorted.bam

samtools index ${SAMPLE}.sorted.bam



### Version 3, Variables and line-continuations, adding a couple parameters to bwa mem
# Note, the variable definitions from above are still fine

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

samtools flagstat \
	${SAMPLE}.sorted.bam 
