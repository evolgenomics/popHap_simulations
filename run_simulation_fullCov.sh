#!/bin/bash
tag=$1
if [ -n "$2" ]; then nInd=$2; else nInd=484; fi
if [ -n "$3" ]; then rho=$3; else rho=10; fi
if [ -n "$4" ]; then clust_threshold=$4; else clust_threshold=0.01; fi

nHap=`$nInd+$nInd | bc`;

echo "***NOTE: The entire pipeline with $nInd samples takes a long time to run. Consider running this in a screen or similar virtual terminals.***";

echo "[run_mspms.sh] Simulate using mspsm $nHap haplotypes, 13505 SNPs, at $rho rho using the tag \"$tag\"...
";

scripts/run_mspms.sh example $nHap 13505 $rho 1

echo "[run_all.sh] Generating truth haploytpes...
";

perl scripts/simulate_haplotypes_ARG.pl $tag 1 $nHap $clust_threshold

echo "[run_all.sh] Starting simulation...
";

scripts/simulate_reads-molecules.sh $tag

scripts/run_STITCH.sh $tag linkedReads bamfiles.linkedReads.$tag.968.list --use_bx_tag=TRUE  
scripts/run_STITCH.sh $tag shortReads bamfiles.linkedReads.$tag.968.list

for datatype in linked short; do 
	scripts/post_STITCH.sh STITCH_out/STITCH_968_${datatype}Reads_$tag; 
	scripts/run_hapcut2.sh STITCH_out/STITCH_968_${datatype}Reads_$tag; 
done; 
