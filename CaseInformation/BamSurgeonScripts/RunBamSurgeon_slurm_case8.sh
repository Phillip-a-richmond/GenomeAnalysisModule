#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 60 Gb of RAM for the whole job
#SBATCH --mem=350G

## Using 16 CPUs
#SBATCH --cpus-per-task=40

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

NSLOTS=78
REFGENOME=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
PICARDJAR=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/opt/BamsurgeonEnvironment/share/picard-2.23.8-0/picard.jar
VCF2BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/GeneBreaker/BenchmarkingTransition/FullSimulation/reformatSimToBamSurgeon.py
GENEBREAKER=/mnt/common/WASSERMAN_SOFTWARE/GeneBreaker/

cd /mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/


# Creating case 8
# Case 8

# Dad
# not needed
cp /mnt/common/OPEN_DATA/POLARIS_RAW/ERR1955420_1.fastq.gz /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/CaseFiles/Case8/Case8_father_R1.fastq.gz
cp /mnt/common/OPEN_DATA/POLARIS_RAW/ERR1955420_2.fastq.gz /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/CaseFiles/Case8/Case8_father_R2.fastq.gz

# Mom
# not needed
cp /mnt/common/OPEN_DATA/POLARIS_RAW/ERR1955443_1.fastq.gz /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/CaseFiles/Case8/Case8_mother_R1.fastq.gz
cp /mnt/common/OPEN_DATA/POLARIS_RAW/ERR1955443_2.fastq.gz /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/CaseFiles/Case8/Case8_mother_R2.fastq.gz



# Proband
# Activate BamSurgeon
BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
source $BAMSURGEON/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BAMSURGEON/opt/BamsurgeonEnvironment

# Define variables
OUTBAM=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case8_proband.bam
OUTBAM_SORTED=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case8_proband.sorted.bam
OUTFQ1=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case8_proband_R1.fastq.gz
OUTFQ2=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Case8_proband_R2.fastq.gz
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR2304556_GRCh38.sorted.bam
INVCFGZ=${GENEBREAKER}/TrainingScenarios/PTPN11_GRCh38_AutosomalDominantDeNovo_Male/proband_PathoVar.vcf.gz
INVCF=${GENEBREAKER}/TrainingScenarios/PTPN11_GRCh38_AutosomalDominantDeNovo_Male/proband_PathoVar.vcf
VARFILE=Case8_proband_bamsurgeon_varfile.tsv

# Create varfile
#gunzip -c $INVCFGZ > $INVCF
#
#python $VCF2BAMSURGEON \
#	-I $INVCF \
#	-O $VARFILE
#
## Mutate Bam
#addsnv.py \
#	-v $VARFILE \
#	-f $INBAM \
#	-r $REFGENOME \
#	-o $OUTBAM \
#	--picardjar $PICARDJAR \
#	--force \
#	--mindepth 5 \
#	-p 40
#
# Samtools sort & index
samtools sort -n $OUTBAM -@ $NSLOTS -o $OUTBAM_SORTED

# Convert to fastq
samtools fastq \
	-@ 40 \
	-1 $OUTFQ1 -2 $OUTFQ2 \
	-0 /dev/null \
	-s /dev/null \
	-c 9 \
	-n $OUTBAM_SORTED
