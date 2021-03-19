#!/bin/bash
tag=$1;
folder=$2;

#cd $folder;
tabix -f stitch.Herato1603.3450000.3550000.vcf.gz 
../../scripts/add_PL.sh stitch.Herato1603.3450000.3550000.vcf.gz
echo $folder;
bcftools reheader -s ../../sourceData/sample.names stitch.Herato1603.3450000.3550000.PL.vcf.gz > stitch.Herato1603.3450000.3550000.PL.reheadered.vcf.gz
tabix stitch.Herato1603.3450000.3550000.PL.reheadered.vcf.gz;
perl ../../scripts/check_truth_ARG.pl $tag
if [ ! -e tmp ]; then mkdir tmp; fi; 
