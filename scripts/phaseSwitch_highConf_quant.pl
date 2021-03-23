use warnings;
use strict;
use Text::Levenshtein qw(distance); 

my $tag = $ARGV[0];
my $linked = $ARGV[1];
system("bcftools query -f \"%CHROM\\t%POS\\n\" stitch.Herato1603.3450000.3550000.PL.HAPCUT2.$linked.vcf.gz > phased_seg_sites.pos");
my @pos = map {chomp; $_ } `cut -f 2 phased_seg_sites.pos`;

my %snp;
my @phased_calls = map{chomp; $_} `bcftools query -f "[%GT\t]\n"  stitch.Herato1603.3450000.3550000.PL.HAPCUT2.$linked.vcf.gz | sed 's/\\t\$//;s/.\\/./9\\|9/g;s/|/\\t/g;s/\\./9\\\t9/g' | datamash transpose | sed 's/\\t//g'`;# 50000.550000.PL.reheadered.vcf.gz`;
my @phased_calls_hiConf = map{chomp; $_} `bcftools query -f "[%PG\t]\n"  stitch.Herato1603.3450000.3550000.PL.HAPCUT2.$linked.vcf.gz | sed 's/\\.\\/\\./9\\\t9/g;s/\\t\$//;s/\\//\\t/g;s/\\./9\\\t9/g' | datamash transpose | sed 's/\\t//g'`;# 50000.550000.PL.reheadered.vcf.gz`;

my @truth=map{chomp; $_} `bcftools query -f "[%GT\t]\n" ../../simHaps/PC062_merged_Herat1603_3.45Mb.simBlock.$tag.GT.vcf.gz -r Herato1603:3450000-3550000 -T phased_seg_sites.pos | sed 's/\\t\$//;s/\\//\\t/g' | datamash transpose | sed 's/\\t//g' `;

my $window = 10;#$ARGV[0];

open(PHASE_CORRECT, ">stitch.Herato1603.3450000.3550000.PL.HAPCUT2.$linked.phaseCorrect.out");#merged.$tag.PL.HAPCUT2.$linked.phasedSwitch.out");
my @switch;
for (my $hap = 0; $hap < $#phased_calls; $hap+= 2) {
	my @dirs;
	my $last_dir = -1;
	$phased_calls[$hap]=~s/\./9/g;
	$phased_calls[$hap+1]=~s/\./9/g;


	my $last_call_0="";
	my $last_call_1="";
	my $last_truth_0="";
	my $last_truth_1="";
	
	my @correct_phased=(0,0,0,0);
	my $het_pos =0;
	foreach my $idx (0..length($phased_calls[$hap])) {
		if (substr($phased_calls_hiConf[$hap], $idx, 1) ne substr($phased_calls_hiConf[$hap+1], $idx, 1)) {# && substr($phased_calls[$hap], $idx, 1) ne substr($phased_calls[$hap+1], $idx, 1)) {
			$correct_phased[0]++ if (substr($phased_calls[$hap], $idx, 1) eq substr($truth[$hap], $idx, 1) && substr($phased_calls[$hap+1], $idx, 1) eq substr($truth[$hap+1], $idx, 1) && (substr($phased_calls[$hap], $idx, 1) ne 9));
			$correct_phased[1]++ if (substr($phased_calls[$hap+1], $idx, 1) eq substr($truth[$hap], $idx, 1) && substr($phased_calls[$hap], $idx, 1) eq substr($truth[$hap+1], $idx, 1) && (substr($phased_calls[$hap+1], $idx, 1) ne 9) );
			$het_pos++;
		}
	}
	if ($correct_phased[1] > $correct_phased[0]) {
		$correct_phased[0]=$correct_phased[1];# = ($correct_phased[2],$correct_phased[3]);
	}
	if ($het_pos > 0) {
		print PHASE_CORRECT "$hap\t".(($correct_phased[0])/($het_pos))."\t$correct_phased[0]\t$het_pos\n";#(length($phased_calls[$hap])*2)."\n";
	} else {
		print PHASE_CORRECT "$hap\tNA\t0\t0\n";#(length($phased_calls[$hap])*2)."\n";
	}
}
close (PHASE_CORRECT);
