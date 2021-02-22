# Activate BamSurgeon
BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
source $BAMSURGEON/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BAMSURGEON/opt/BamsurgeonEnvironment
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR1955507_GRCh38.sorted.bam
OUTBAM=ERR1955507_GRCh38_PlusABCD1.sorted.bam
REFGENOME=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
VARFILE=ABCD1_varfile.txt
PICARDJAR=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/opt/BamsurgeonEnvironment/share/picard-2.23.8-0/picard.jar

addsnv.py \
	-v $VARFILE \
	-f $INBAM \
	-r $REFGENOME \
	-o $OUTBAM \
	--picardjar $PICARDJAR \
	--force \
	--mindepth 5 \
	-p 16

