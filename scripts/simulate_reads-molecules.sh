#!/bin/bash
tag=$1;

nSNP=`bcftools query -f "%POS\n" sourceData/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz | wc -l`;

echo "[simulate_reads-molecules.sh] Generating summaries of haplotype frequencies and haplotype diversity statistics...
";

perl scripts/extract_hapFreq.pl $tag
perl scripts/block_switch.pl $tag &
perl scripts/uniq_haps.pl $tag > simHaps/PC062_merged_Herat1603_3.45Mb.simBlock.$tag.hapDiv
#

echo "[simulate_reads-molecules.sh] Generating truth set VCFs
";
cd simHaps
cut -f 2- PC062_merged_Herat1603_3.45Mb.simBlock.$tag.haps | sed 's/./&\t/g' | datamash transpose | awk '{str=$1"/"$2;for(i=3;i<NF;i+=2) {str=str"\t"$i"/"$(i+1)}; print str}'  | awk 'NF > 1' > PC062_merged_Herat1603_3.45Mb.simBlock.$tag.GT
cat ../sourceData/vcf.header <(paste ../sourceData/header_cols.vcf PC062_merged_Herat1603_3.45Mb.simBlock.$tag.GT) | bcftools view - -Oz -o PC062_merged_Herat1603_3.45Mb.simBlock.$tag.GT.vcf.gz
tabix PC062_merged_Herat1603_3.45Mb.simBlock.$tag.GT.vcf.gz
cd ..

#
echo "[STITCH_simulate.sh] Simulating reads...
";

for i in `seq 0 967`; do dwgsim -e 0.001 -E 0.005 -d 250 -s 150 -H -o 1 -C 10 -1 150 -2 150 -x sourceData/focal_region.bed altFa/PC062_merged_Herat1603_3.45Mb.simBlock.hap_$i.$tag.fa simFastq/hap_$i.$tag 2> log/$tag.dwgsim.stderr.log  ; done
echo "[STITCH_simulate.sh] Placing simulated reads against canonical reference assembly...
";
if [ ! -e simBam/$tag ]; then mkdir simBam/$tag; fi
for i in `seq 0 967`; do bwa mem -O 20 -E 20 -B 1 sourceData/helera1_demo_Herato1603.fa simFastq/hap_$i.$tag.bwa.read1.fastq.gz simFastq/hap_$i.$tag.bwa.read2.fastq.gz -O BAM | samtools sort - -O BAM -o simBam/$tag/hap_$i.$tag.bam; done 2> log/$tag.dwgsim.stderr.log
####

cd simBam
for i in `seq 0 967`; do echo $i; perl ../scripts/simMols.pl $i $tag; done
for i in `seq 0 2 967`; do odd=`echo "$i+1" | bc`; echo $i"	"$odd; samtools merge -c -p -f $tag/hap_$i.linkedReads.diploid.$tag.bam $tag/hap_$i.linkedReads.haploid.$tag.bam $tag/hap_$odd.linkedReads.haploid.$tag.bam; samtools index $tag/hap_$i.linkedReads.diploid.$tag.bam; done
##
cd ..
ls simBam/$tag/hap_*.linkedReads.diploid.$tag.bam | sort -k 1,1V > simBam/bamfiles.linkedReads.$tag.968.list

