#!/bin/bash

tag=$1;

if [ ! -e simBam/${tag}_subsample ]; then mkdir simBam/${tag}_subsample; fi

cd simBam;
for i in `seq 0 2 967`; do samtools view $tag/hap_$i.linkedReads.diploid.$tag.bam | sed 's/BX:Z:/\nBX:Z:/' | grep BX:Z | sort | uniq > ${tag}_subsample/hap_$i.linkedReads.diploid.$tag.bx; done
for frac in 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75; do 
	echo "[mol_subsample.sh] Subsampling at ${frac}x...
";
	perl ../scripts/mol_subsample.pl $tag $frac
	ls ${tag}_subsample/hap_$i.linkedReads.diploid.$tag.$frac.bam | sort -k 1,1V > bamfiles.linkedReads.$tag.$frac.968.list
done
rm hap*.linkedReads.diploid.$tag.bx
cd ../
