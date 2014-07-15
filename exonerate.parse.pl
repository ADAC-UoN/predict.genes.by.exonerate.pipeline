#!/usr/bin/perl 
use warnings;
use strict;

use Getopt::Long;


# Richard Emes University of Nottingham 2010
my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
(C) R D Emes University of Nottingham 2010

parse exonerate

USAGE:
-f	exonerate file
-o	output file
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
" ;


my ($file, $output);

GetOptions(
        'f|fasta:s'     => \$file,
        'o|output:s'     => \$output,	
         );


if( ! defined $file) {
print "$usage\nWARNING: Cannot proceed without fasta file to process\n\n"; exit;
}
if( ! defined $output) {
print "$usage\nWARNING: Cannot proceed without exonerate output file\n\n"; exit;
}

####################################################



#~ PARSE OUTPUT
#####

open FILE, "$file";
open OUT, ">$output";
print OUT "Protein Seq\tscore\tquery_start\tquery_end\tquery_length\ttarget_start\ttarget_end\ttarget_length\tStrand\tIntrons\n";

my $query = "N_A";
my $score = "N_A";
my $t_end = "";
my $t_start = "";
my $q_end = "";
my $q_start = "";
my $strand = "";
my $query_length = 0;
my $target_length = 0;
my $intron_count = 0;
my $flag_hit = 0;
while (<FILE>)
{
chomp $_;
if ($_ =~ /^C4/){}
elsif ($_ =~ /^----/){}
elsif ($_ =~ /^-- completed exonerate analysis/) {}
elsif ($_ =~ /^Command line\:/){}
elsif ($_ =~ /^Hostname/){}


elsif ($_ =~ /^\s+Query\:\s+(.*)/){$query = $1;  $query =~ s/\s+$//;}
elsif ($_ =~ /^\s+Raw score\:\s+(\d+)/) {$score = $1}
elsif ($_ =~ /^\s+Query range\:\s+(\d+)\s+\-\>\s+(\d+)\s*/)
		{
		$q_start = $1;
		$q_end = $2;
		$query_length = ($q_end - $q_start);
		}
elsif ($_ =~ /^\s+Target range\:\s+(\d+)\s+\-\>\s+(\d+)\s*/)
		{
		$t_start = $1;
		$t_end = $2;
		$target_length = ($t_end - $t_start);
		if
		($t_end < $t_start)
			{
			$strand = "-";
			my $left = $t_end;
			my $right = $t_start;
			$t_end = $right;
			$t_start = $left;
			$target_length = ($t_end - $t_start);
			}
		else 
		{
		$strand = "+";
		}
		
		}
elsif ($_ =~ /> T/ || /\:\s+Targ/)
	{
	$intron_count++;
	}

elsif ($_ =~ /^\>/)
	{
	print OUT "$query\t$score\t";
	print OUT "$q_start\t$q_end\t$query_length\t$t_start\t$t_end\t$target_length\t$strand\t$intron_count\n";
	$query = "N_A"; $score = "N_A";
	$intron_count = 0;
	}


else {}
}
