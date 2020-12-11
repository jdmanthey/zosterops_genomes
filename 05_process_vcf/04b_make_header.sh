# make a vcf file that contains just the header of the vcf files (all are in common)

grep '#' CHRLGE22.recode.vcf > header.vcf
