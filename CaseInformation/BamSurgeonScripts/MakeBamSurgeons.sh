VCF2BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/GeneBreaker/BenchmarkingTransition/FullSimulation/reformatSimToBamSurgeon.py

GENEBREAKER=/mnt/common/WASSERMAN_SOFTWARE/GeneBreaker/


# Case 1
# dad
invcfgz=${GENEBREAKER}/TrainingScenarios/WAS_GRCh38_XLinkedRecessiveHemizygous_Male/father_PathoVar.vcf.gz
invcf=${GENEBREAKER}/TrainingScenarios/WAS_GRCh38_XLinkedRecessiveHemizygous_Male/father_PathoVar.vcf
gunzip -c $invcfgz > $invcf

python $VCF2BAMSURGEON \
	-I $invcf \
	-O Case1_father_bamsurgeon_varfile.tsv

# mom
invcfgz=${GENEBREAKER}/TrainingScenarios/WAS_GRCh38_XLinkedRecessiveHemizygous_Male/mother_PathoVar.vcf.gz
invcf=${GENEBREAKER}/TrainingScenarios/WAS_GRCh38_XLinkedRecessiveHemizygous_Male/mother_PathoVar.vcf
gunzip -c $invcfgz > $invcf

python $VCF2BAMSURGEON \
	-I $invcf \
	-O Case1_mother_bamsurgeon_varfile.tsv

# proband
invcfgz=${GENEBREAKER}/TrainingScenarios/WAS_GRCh38_XLinkedRecessiveHemizygous_Male/proband_PathoVar.vcf.gz
invcf=${GENEBREAKER}/TrainingScenarios/WAS_GRCh38_XLinkedRecessiveHemizygous_Male/proband_PathoVar.vcf
gunzip -c $invcfgz > $invcf

python $VCF2BAMSURGEON \
	-I $invcf \
	-O Case1_proband_bamsurgeon_varfile.tsv


