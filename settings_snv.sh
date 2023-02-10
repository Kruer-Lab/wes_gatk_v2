# General Parameters
NCORE=32
BUILD="GRCh38"

#pipeline directory
PIPELINE="/media/sheetalshetty/Data2/Pipeline_V2/"

# Samples, data and pedigrees
ALL_SAMPLES_FASTQ_DIR="/media/sheetalshetty/Data2/Pipeline_V2/Data/Fastqs"
DATA="$PIPELINE/Data/"
INPUTS="${DATA}/Inputs/"
ALL_SAMPLES_OUTPUT_DIR="/media/sheetalshetty/Data2/Pipeline_V2/Data/Outputs"
ANNOVAR_DIR="/media/sheetalshetty/Data2/Pipeline_V2/Tools/Annovar"
PEDIGREES_DIR="/media/sheetalshetty/Data2/Pipeline_V2/Data/Pedigrees"


# Reference genome (downloaded from gs://gcp-public-data--broad-references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy)
REF="${INPUTS}/Ref/Homo_sapiens_assembly38_noALT_noHLA_noDecoy/v1/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"
DBSNP="${INPUTS}/Bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
INTERVAL_LIST="${INPUTS}/Interval_List/Exome-IDT_V1V2_span50bp.bed"
#GFF="/media/sheetalshetty/Data2/Pipeline_V2/Data/Inputs/Bundle/gencode.v42.basic.annotation.gff3.gz"
GFF="/media/sheetalshetty/Data2/Pipeline_V2/Data/Inputs/Bundle/Homo_sapiens.GRCh38.108.chr_added.gff3"
#GFF="/media/sheetalshetty/Data2/Pipeline_V2/Data/Inputs/Bundle/Homo_sapiens.GRCh38.108.chr.gff3"
#Tools
GATK="/media/sheetalshetty/Data2/Pipeline_V2/Tools/GATK/gatk-4.2.6.1/gatk"
IGVTOOLS="/media/sheetalshetty/Data2/Pipeline_V2/Tools/IGV/IGV_2.14.1/igvtools"
#SNPEFF="/media/sheetalshetty/Data2/Pipeline_V2/Tools/snpEff/snpEff.jar"
#VEP="/media/sheetalshetty/Data2/Pipeline_V2/Tools/VEP/ensembl-vep/vep"

# FILTER CUTOFFS 
# Denovo variants
export MAP_QUAL_DEN=59	# mapping quality 
export COV_PROB_DEN=10	# proband coverage for denovo variants
export MAD_PROB_DEN=5	#proband minor allele depth for denovo variant
export MA_FRAC1_PROB_DEN=0.2	# proband minor allele fraction threhsod 1 for denovo variant
export MA_FRAC2_PROB_DEN=0.28	# proband minor allele fraction threshold 2 for denovo variant
export MAD_LINE_DEN=10	# alernate allele margin for denovo variant
export MAF_DEN=0.0005	# minor allele frequency threshod for denovo variants
export COV_PAR_DEN=10	# parent coverage for denovo variants
export MAD_PROB_DEN=5 	#proband minor allele depth for denovo variants
export MAD_FRAC_PAR_DEN=0.035  #parent minor allele depth fraction for denovo variants

# Recessive variants
export MAP_QUAL_REC=59 # mapping quality for recessive variants
export GT_QUAL_REC=20.00	# genotype quality for recessive variants
export MAF_REC=0.005 	# minor allele frequency for recessive variants
export COV_PROB_REC=8	# proband coverage for recessive variants
export CADD_THR=20.00	# minimum CADD score 

# Dominant variants
export MAP_QUAL_DOM=59 # mapping quality for dominant variants
export GT_QUAL_DOM=20	# genotype quality
export COV_PROB_DOM=8	# proband coverage for dominant variants
export MAF_DOM=0.0005	# minor allele frequency for dominant variants
export MA_FRAC1_PROB_DOM=0.2	# proband minor allele fraction threhsod 1 for dominant variant
export MA_FRAC2_PROB_DOM=0.28	# proband minor allele fraction threshold 2 for dominant variant
export MAD_LINE_DOM=10	# alernate allele margin for dominant variant


#Xlink variants
export MAF_XLINK=0.0005	# minor allele frequency for dominant variants










