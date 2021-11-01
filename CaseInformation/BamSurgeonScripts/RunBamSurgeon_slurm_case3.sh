#!/bin/bash

#SBATCH --partition=dev_q

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 60 Gb of RAM for the whole job
#SBATCH --mem=80G

## Using 16 CPUs
#SBATCH --cpus-per-task=10

## Running for a max time of 48 hours
#SBATCH --time=48:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

##########
# Set up #
##########

NSLOTS=10
REFGENOME=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
PICARDJAR=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/opt/BamsurgeonEnvironment/share/picard-2.23.8-0/picard.jar
VCF2BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/GeneBreaker/BenchmarkingTransition/FullSimulation/reformatSimToBamSurgeon.py
GENEBREAKER=/mnt/common/WASSERMAN_SOFTWARE/GeneBreaker/

cd /mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/


# Case 3

###########
# Proband #
###########
# Activate BamSurgeon
BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
source $BAMSURGEON/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BAMSURGEON/opt/BamsurgeonEnvironment

# Define variables
OUTBAM=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_proband.bam
OUTBAM_SORTED=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_proband.sorted.bam
OUTFQ1=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_proband_R1.fastq
OUTFQ2=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_proband_R2.fastq
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR2304565_GRCh38.sorted.bam
INVCFGZ=${GENEBREAKER}/TrainingScenarios/MALT1_GRCh38_AutosomalRecessiveHomozygous_Female/proband_PathoVar.vcf.gz
INVCF=${GENEBREAKER}/TrainingScenarios/MALT1_GRCh38_AutosomalRecessiveHomozygous_Female/proband_PathoVar.vcf
VARFILE=Case3_proband_bamsurgeon_varfile.tsv

# Create varfile
gunzip -c $INVCFGZ > $INVCF

python $VCF2BAMSURGEON \
	-I $INVCF \
	-O $VARFILE

# Mutate Bam
addindel.py \
	-v $VARFILE \
	-f $INBAM \
	-r $REFGENOME \
	-o $OUTBAM \
	--picardjar $PICARDJAR \
	--force \
	--mindepth 5 \
	-p $NSLOTS

# Samtools sort & index
samtools sort -n $OUTBAM -@ $NSLOTS -o $OUTBAM_SORTED
#samtools index $OUTBAM_SORTED

samtools sort $OUTBAM -@ $NSLOTS -o ${OUTBAM[::-4]}_position.sorted.bam

# Convert to fastq
samtools fastq \
	-@ $NSLOTS \
	-1 $OUTFQ1 -2 $OUTFQ2 \
	-0 /dev/null \
	-s /dev/null \
	-n $OUTBAM_SORTED

gzip $OUTFQ1
gzip $OUTFQ2

#######
# Dad #
#######
# Activate BamSurgeon
BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
source $BAMSURGEON/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BAMSURGEON/opt/BamsurgeonEnvironment

# Define variables
OUTBAM=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_father.bam
OUTBAM_SORTED=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_father.sorted.bam
OUTFQ1=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_father_R1.fastq
OUTFQ2=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_father_R2.fastq
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR1955499_GRCh38.sorted.bam
INVCFGZ=${GENEBREAKER}/TrainingScenarios/MALT1_GRCh38_AutosomalRecessiveHomozygous_Female/father_PathoVar.vcf.gz
INVCF=${GENEBREAKER}/TrainingScenarios/MALT1_GRCh38_AutosomalRecessiveHomozygous_Female/father_PathoVar.vcf
VARFILE=Case3_father_bamsurgeon_varfile.tsv

# Create varfile
gunzip -c $INVCFGZ > $INVCF

python $VCF2BAMSURGEON \
	-I $INVCF \
	-O $VARFILE

# Mutate Bam
addindel.py \
	-v $VARFILE \
	-f $INBAM \
	-r $REFGENOME \
	-o $OUTBAM \
	--picardjar $PICARDJAR \
	--force \
	--mindepth 5 \
	-p $NSLOTS

# Samtools sort & index
samtools sort -n $OUTBAM -@ $NSLOTS -o $OUTBAM_SORTED
#samtools index $OUTBAM_SORTED

# Convert to fastq
samtools fastq \
	-@ $NSLOTS \
	-1 $OUTFQ1 -2 $OUTFQ2 \
	-0 /dev/null \
	-s /dev/null \
	-n $OUTBAM_SORTED


gzip $OUTFQ1
gzip $OUTFQ2



#######
# Mom #
#######
# Activate BamSurgeon
BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
source $BAMSURGEON/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BAMSURGEON/opt/BamsurgeonEnvironment

# Define variables
OUTBAM=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_mother.bam
OUTBAM_SORTED=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_mother.sorted.bam
OUTFQ1=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_mother_R1.fastq
OUTFQ2=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case3_mother_R2.fastq
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR1955435_GRCh38.sorted.bam
INVCFGZ=${GENEBREAKER}/TrainingScenarios/MALT1_GRCh38_AutosomalRecessiveHomozygous_Female/mother_PathoVar.vcf.gz
INVCF=${GENEBREAKER}/TrainingScenarios/MALT1_GRCh38_AutosomalRecessiveHomozygous_Female/mother_PathoVar.vcf
VARFILE=Case3_mother_bamsurgeon_varfile.tsv

# Create varfile
gunzip -c $INVCFGZ > $INVCF

python $VCF2BAMSURGEON \
	-I $INVCF \
	-O $VARFILE

# Mutate Bam
addindel.py \
	-v $VARFILE \
	-f $INBAM \
	-r $REFGENOME \
	-o $OUTBAM \
	--picardjar $PICARDJAR \
	--force \
	--mindepth 5 \
	-p $NSLOTS

# Samtools sort & index
samtools sort -n $OUTBAM -@ $NSLOTS -o $OUTBAM_SORTED
#samtools index $OUTBAM_SORTED

# Convert to fastq
samtools fastq \
	-@ $NSLOTS \
	-1 $OUTFQ1 -2 $OUTFQ2 \
	-0 /dev/null \
	-s /dev/null \
	-n $OUTBAM_SORTED


gzip $OUTFQ1
gzip $OUTFQ2
