#!/bin/bash
# Brendan Reardon
# Van Allen lab, Dana-Farber Cancer Institute
# build_exac.sh, Version 1.4.0
# Last updated 5 April 2018

v=1.4
outname='exac.lite-pass'$v'-r1.txt'

if [ -f $1 ]
then
	exac_vcf=$1
else
	exac_vcf=~/storage/exac/ExAC.r1.sites.vep.vcf
fi

if [ -f $2 ]
then
	reference_genome=$2
else
	reference_genome=~/storage/hg19/Homo_sapeisn_assembly19.fasta
fi

if [ -f $3 ]
then
	gatk=$3
else
	gatk=~/software/gatk/GenomeAnalysisTK-3.8.jar
fi

java -jar $gatk -T VariantsToTable -R $reference_genome -V $exac_vcf -o $outname -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AF -F AC -F AC_AFR -F AC_AMR -F AC_EAS -F AC_FIN -F AC_NFE -F AC_OTH -F AC_SAS -F AN -F AN_AFR -F AN_AMR -F AN_EAS -F AN_FIN -F AN_NFE -F AN_OTH -F AN_SAS
