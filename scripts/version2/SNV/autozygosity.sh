#!/bin/bash
source settings_snv.sh


CHILD_DIR=$1
#TRIO_ID=$(echo $CHILD_DIR | awk -F "-" '{print $1}')
TRIO_ID=$(echo $CHILD_DIR | awk -F "-" '{print $1"-"$2}')
MOTHER_DIR=${TRIO_ID}-001-U
FATHER_DIR=${TRIO_ID}-002-U
#MOTHER_DIR=${TRIO_ID}-001_G-U
#FATHER_DIR=${TRIO_ID}-002_G-U
#CHILD_DIR=${TRIO_ID}-003-A

# End of Section1

# Section2: uncomment this section and comment the section1 if the file formats DO NOT FOLLOW our trio file name convention
#echo -e "Mother directory name:"
#read MOTHER_DIR
#echo -e "Father directory name:"
#read FATHER_DIR
#echo -e "Child directory name:"
#read CHILD_DIR

# End of section2

echo ${ALL_SAMPLES_OUTPUT_DIR}

OUTPUTS_CHILD="${ALL_SAMPLES_OUTPUT_DIR}/${CHILD_DIR}"
SAMPLE_ID_CHILD=$(ls ${OUTPUTS_CHILD}/*.raw.g.vcf | awk -F".raw.g.vcf" '{print $1}'| awk -F"/" '{print $(NF)}')
CHILD="${OUTPUTS_CHILD}/${SAMPLE_ID_CHILD}"


OUTPUTS_MOTHER="${ALL_SAMPLES_OUTPUT_DIR}/${MOTHER_DIR}"
SAMPLE_ID_MOTHER=$(ls ${OUTPUTS_MOTHER}/*.raw.g.vcf | awk -F".raw.g.vcf" '{print $1}'| awk -F"/" '{print $(NF)}')
MOTHER=${OUTPUTS_MOTHER}/${SAMPLE_ID_MOTHER}


OUTPUTS_FATHER="${ALL_SAMPLES_OUTPUT_DIR}/${FATHER_DIR}"
SAMPLE_ID_FATHER=$(ls ${OUTPUTS_FATHER}/*.raw.g.vcf | awk -F".raw.g.vcf" '{print $1}' | awk -F"/" '{print $(NF)}')
FATHER=${OUTPUTS_FATHER}/${SAMPLE_ID_FATHER}

OUTPUTS_TRIO="${ALL_SAMPLES_OUTPUT_DIR}/${CHILD_DIR}_Trio"
TRIO="${OUTPUTS_TRIO}/${SAMPLE_ID_CHILD}"

BAM_CHILD="${OUTPUTS_CHILD}/${SAMPLE_ID_CHILD}.bqrecal.reads.bam"
BAM_MOTHER="${OUTPUTS_MOTHER}/${SAMPLE_ID_MOTHER}.bqrecal.reads.bam"
BAM_FATHER="${OUTPUTS_FATHER}/${SAMPLE_ID_FATHER}.bqrecal.reads.bam"

#CHILD_SAMPLE_NAME=$(echo $CHILD_DIR | awk -F "_" '{print $1"_"$2}')
#MOTHER_SAMPLE_NAME=$(echo $MOTHER_DIR | awk -F "_" '{print $1"_"$2}')
#FATHER_SAMPLE_NAME=$(echo $FATHER_DIR | awk -F "_" '{print $1"_"$2}')

CHILD_SAMPLE_NAME=$CHILD_DIR
MOTHER_SAMPLE_NAME=$MOTHER_DIR 
FATHER_SAMPLE_NAME=$FATHER_DIR

PEDIGREE="${PEDIGREES_DIR}/pedigree.${CHILD_SAMPLE_NAME}.p2.ped"
SAMPLE_LIST=($CHILD_SAMPLE_NAME $MOTHER_SAMPLE_NAME $FATHER_SAMPLE_NAME)
SAMPLE_LIST_SORTED=( $(
for SAMP in "${SAMPLE_LIST[@]}"
do
     echo "$SAMP"
done | sort) )

echo "Sorted trio: ${SAMPLE_LIST_SORTED[@]}"

IDX=0
for SAMP in "${SAMPLE_LIST_SORTED[@]}"	
do
	if [ "$CHILD_SAMPLE_NAME" == "$SAMP" ] 
	then		
		let "CHILD_IDX=$IDX+9"
	elif [ "$MOTHER_SAMPLE_NAME" == "$SAMP" ]
	then
		let "MOTHER_IDX=$IDX+9"
	elif [ "$FATHER_SAMPLE_NAME" == "$SAMP" ]
	then
		let "FATHER_IDX=$IDX+9"
	else
		echo "Something wrong with sample names!!"
		exit
	fi
	let "IDX++"
done

export CHILD_IDX
export MOTHER_IDX
export FATHER_IDX







BEAGLE=/media/sheetalshetty/Data2/Pipeline_V2/Tools/Beagle/beagle.22Jul22.46e.jar
#GVCF="${ALL_SAMPLES_OUTPUT_DIR}/$1/$1.vcf"
GVCF="/media/sheetalshetty/Data2/Pipeline_V2/Data/Outputs/KRU-F1116-003-A/KRU-F1116-003-A-DNA-01_168_025_S33_L004.raw.g.vcf"
#$GATK --java-options "-Xmx4g" GenotypeGVCFs -R $REF --variant ${CHILD}.raw.g.vcf -O ${CHILD}.vcf
#OUTDIR="${DATA}/Autozygosity_Output/"
PLINK=/media/sheetalshetty/Data2/Pipeline_V2/Tools/Plik/plink


#ls ${CHILD}.vcf
#ls $BEAGLE
ls $REF
echo $CHILD_SAMPLE_NAME

#java -Xmx16g -jar $BEAGLE gt=${TRIO}.trio.raw.vcf out=${TRIO} 



# Convert the phased VCF file to a PLINK binary file format
#$PLINK --vcf ${TRIO}.vcf --make-bed --out ${TRIO}

# Identify ROHs using PLINK
####$PLINK --bfile ${TRIO}.fam --homozyg --homozyg-window-kb 500 --homozyg-window-het 1 --homozyg-window-snp 50 --out ${TRIO}.homozygosity


#$PLINK --bfile ${TRIO}.bed --homozyg --out ${TRIO}.roh

# Extract ROHs from the PLINK output file
#grep -v '#' homozygosity.hom | awk '{print $1,$4,$5}' > rohs.bed





# Input file containing genotype data
#genotype_file="${TRIO}.trio.raw.vcf"

# Output file for homozygous regions
#output_file="${TRIO}.roh.bed"

# Minimum length of homozygous region
#min_length=100

# Call homozygous regions using Beagle
#####java -jar $BEAGLE gt=$genotype_file out=$output_file homozygous=true minlength=$min_length

#java -jar $BEAGLE gt=$genotype_file out=$output_file 

# Extract regions of homozygosity from output file
#awk '$6 == "HOM" {print $1"\t"$2"\t"$3}' $output_file > $output_file.homozygous

#ls "Phased VCFFFFFF: ${TRIO}.phased.vcf"

#java -Xmx16g -jar $BEAGLE gt=${TRIO}.trio.raw.vcf out=${TRIO}.phased
# Create a BED file with ROH regions
#ls ${TRIO}.phased.vcf
#tabix -p vcf ${TRIO}.phased.vcf.gz
#bcftools query -f '%CHROM\t%POS\t%POS\n' ${TRIO}.phased.vcf.gz | awk '$3-$2 > 100 { print $0 }' > ${TRIO}.roh.bed

#bcftools +fill-tags -O b -o ${TRIO}.roh.bed ${TRIO}.phased.vcf.gz

#bcftools roh -m 100 -q 0.05 -R ${TRIO}.phased.vcf.gz -O v -o ${TRIO}.output.vcf

#bcftools roh -G -i ${TRIO}.phased.vcf.gz -o ${TRIO}.out.txt




# Convert the VCF file to a Beagle-compatible format
#java -Xmx2g -jar $BEAGLE gt=${TRIO}.trio.raw.vcf out=${TRIO}.bgl

# Phase the data using Beagle
#java -Xmx2g -jar $BEAGLE gt=${TRIO}.bgl.vcf.gz out=${TRIO}.phased nthreads=4

# Extract regions of homozygosity from the phased data
#java -Xmx2g -jar $BEAGLE gt=${TRIO}.phased.vcf.gz roh=true out=${TRIO}.roh

# Filter the regions of homozygosity to remove any regions likely to be errors or artifacts
#java -Xmx2g -jar beagle.28Sep18.793.jar gt=genetic_data_roh.vcf.gz filter=true out=genetic_data_filtered_roh

# Analyze the filtered regions of homozygosity
#java -Xmx2g -jar beagle.28Sep18.793.jar gt=genetic_data_filtered_roh.vcf.gz stats=true out=genetic_data_roh_stats



#$PLINK --vcf ${TRIO}.phased.vcf.gz --make-bed --out ${TRIO}.phased
#plink --bfile phased --homozyg --out homozygosity

# Analyze the BEAGLE output and extract regions of homozygosity
# (Example code using awk, this will vary based on the format of the output file)
#awk '$5 == 1 { print $0 }' ${TRIO}.phased.vcf > ${TRIO}.homozygosity_regions.txt


#java -Xmx4g -jar $BEAGLE gt=${TRIO}.phased.vcf out=${TRIO}.beagle_output
#grep -P '\t1\t1\t' ${TRIO}.beagle_output.vcf > ${TRIO}.homozygosity.txt


#$PLINK --vcf ${TRIO}.trio.raw.vcf --homozyg --homozyg-window-kb 1000 --homozyg-density 200 --homozyg-kb 500 --out ${TRIO}.output
$PLINK --bfile ${TRIO}.trio.raw.vcf --homozyg --homozyg-kb 500 --homozyg-snp 250 --out ${TRIO}.ROH