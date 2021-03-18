use warnings;
use strict;

my $tag = $ARGV[0];
my $frac = $ARGV[1];

foreach my $sample (0..483) {
	print "[mol_subsample.pl] Handling sample - $sample\n";
	my $i = $sample * 2;
	my $folder = $tag."_subsample";
	my $molNum = `cat $folder/hap_$i.linkedReads.diploid.$tag.bx | wc -l`;
	chomp($molNum);
	$molNum = int($molNum * $frac + 0.5); 
	my %bx = map {chomp;$_ => 1} `shuf -n $molNum $folder/hap_$i.linkedReads.diploid.$tag.bx`;
	open (OUT, " | samtools view - -O BAM -o $folder/hap_$i.linkedReads.diploid.$tag.$frac.bam");
	open (IN, "samtools view $tag/hap_$i.linkedReads.diploid.$tag.bam -h |");		
	while (<IN>) {
		print OUT $_ if (/^@/);
		if (/(BX:Z:\S+)/) {
			if (exists($bx{$1})) {
				print OUT $_; 
			}
		}
	}
	close (OUT);
	system("samtools index $folder/hap_$i.linkedReads.diploid.$tag.$frac.bam");
}
