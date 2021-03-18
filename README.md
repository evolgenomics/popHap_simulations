# Simulations
The purpose of this repo is to allow users to simulate our "population haplotyping" pipeline under a range of diversity conditions.
The focal example we use here is from the *Heliconius erato* system. 

Broadly speaking, this simulation scheme creates a set of diploid individuals with distinct haplotypes and allows recombination between them. It then simulates short reads (with rare sequencing errors) against each haplotype and map these against a canonical reference genome assembly. A separate Perl script would simulate molecules (each marked by the SAM tag "BX" for beadTag) and choose from among the set of mapped reads (as a BAM file) to form linked reads. Then pairs of BX-tagged haploid BAM files would be merged into diploid BAMs. 

At this point the set of BAM files can be optionally subsampled by molecules (`BX` tag), to reflect varying sequencing coverage.

The entire set of `BX`-tagged, diploid BAM files (optionally subsampled) would then be processed by STITCH to perform population-level statistical phasing, genotyping and imputation. The specific parameters here can be modified, but are set here to replicate those used for the empirical study.

Please note, that here we *do not* include an extra *de novo* SNP calling step, because we are mainly interested in the genotyping and phasing performance. Should any user be interested in that, common tools such as samtools or GATK, among others, can be included here to generate a SNP set for STITCH. Instead, we simply feed STITCH with the truth set of SNPs.

At each of the following STITCH and HAPCUT2 steps, there is the option to include, or ignore `BX` tag. This is to specifically test the contribution of linkage information to haplotype reconstruction/genotyping/imputation and phasing.

The resulting VCF file generated by STITCH would then be processed to add `PL` tag, move the `GT` tag into the `PG` ("prior genotype") tag and instead populate the `GT` tag with all forced calls. In other words, "high confidence" calls originally emitted by STITCH are stored in the `PG` tag, whereas `GT` is now the most likely genotype calls, regardless of confidence. 

Additionally, each individual would be subjected to molecular phasing at its heterozygous sites. Helper scripts would then merge all individual phased VCFs into a single phased VCF file.

These phased genotypes are then compared by Perl scripts against the truth genotypes for genotyping and phasing accuracy. Genotype accuracy is reported as global accuracy (defined as proportion of identical gentoype calls against true genotypes) and within each non-recombinant segment (according to the tree sequence as generated by `mspms`). Phasing accuracy is scored at high-confidence sites (using the `PG` tag), as a percentage of correctly phased genotypes out of all high-confidence sites. Switch errors are quantified as switches between adjacent windows (e.g., 10 SNPs). This latter code works, but lacks efficiency and needs some refinement. 

Overall summaries were aggregated by bash scripts.

A note on the focal region: for the purpose of this simulation, we will be creating a 1 Mbp region, chosen from a putatively neutral, background location.

# Dependencies
To run this simulation, you'll need:
- [msprime](http://https://github.com/tskit-dev/msprime)
- [dwgsim](https://github.com/nh13/DWGSIM)
- [bwa](https://github.com/lh3/bwa)
- [samtools](https://github.com/samtools/samtools)
- [tabix](https://github.com/samtools/tabix.git)
- [bcftools](https://github.com/samtools/bcftools)
- [STITCH](https://github.com/rwdavies/STITCH)
- [HAPCUT2](https://github.com/vibansal/HapCUT2)
- [datamash](https://www.gnu.org/software/datamash/)

# Quick start
To help you get familiar with the simulation scheme, I have written a quick wrapper script to download the relevant sample data and basic assembly, followed by an example `run_all.sh` script.
```
git clone https://github.com/evolgenomics/popHap_simulations.git
cd popHap_simulations
./setup.sh
./run_all.sh
```

# Pipeline map
Here is a quick overview of the pipelines used in the paper:
![Pipeline](Pipelines_simulations_1.png?raw=true "Simulation pipeline")

