#!/bin/bash
source settings_snv.sh
# Section1: uncomment this section and comment the section2 if the file formats FOLLOW our trio file name convention
#echo -e "Trio Child ID (e.g M_F071-003-A):"
#read CHILD_DIR
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

ls ${CHILD}.bqrecal.bam
ls $REF
# First, index the BAM file
#samtools index ${CHILD}.bqrecal.bam

# Next, extract the read depth for the X and Y chromosomes
#samtools depth -r chrX ${CHILD}.bqrecal.bam > ${TRIO}.x_depth.txt
#samtools depth -r chrY ${CHILD}.bqrecal.bam > ${TRIO}.y_depth.txt


#bcftools stats -F $REF -r chrX ${CHILD}.bqrecal.bam > x_depth.txt
#bcftools stats -F $REF -r chrY ${CHILD}.bqrecal.bam > y_depth.txt

# Now, compute the ratio of X to Y read depth
#bcftools stats -F ${TRIO}.x_depth.txt -I ${TRIO}.y_depth.txt > ${TRIO}.stats.txt


#y_coverage=$(samtools depth -a -r chrY -l $INTERVAL_LIST ${CHILD}.bqrecal.bam | awk '{sum+=$3} END {print sum}')

# Count the number of bases in the Y chromosome within the intervals
#y_length=$(cat $INTERVAL_LIST | awk '{sum+=$3-$2+1} END {print sum}')



# extract the read depth for the Y chromosome
#samtools depth -r chrY ${CHILD}.bqrecal.bam -b $INTERVAL_LIST> ${TRIO}.read_depth.txt

# use awk to calculate the Y chromosome read depth
#awk '{if ($1=="chrY") {y_depth+=$3}} END {if (y_depth > 5) {print "Male"} else {print "Female"}}' ${TRIO}.read_depth.txt



bedtools coverage -a $INTERVAL_LIST -b ${CHILD}.bqrecal.bam -mean  -chrom Y | awk '{sum+=$NF} END {print sum/NR}'
