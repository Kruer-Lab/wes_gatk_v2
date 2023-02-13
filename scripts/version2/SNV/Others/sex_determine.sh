
#!/bin/bash
source settings_snv.sh


# First, index the BAM file
samtools index sample.bam

# Next, extract the read depth for the X and Y chromosomes
samtools depth -r X sample.bam > x_depth.txt
samtools depth -r Y sample.bam > y_depth.txt

# Now, compute the ratio of X to Y read depth
bcftools stats -F x_depth.txt -I y_depth.txt > stats.txt