#!/usr/bin/perl 
use warnings;
use strict;

use Getopt::Long;


# Richard Emes University of Nottingham 2010
my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
(C) R D Emes University of Nottingham 2010

predict genes from genomic DATA using exonerate

HAVE YOU CHECKED THE EXONERATE COMMANDS

USAGE:
-f	fasta file (protein sequence of genes of interest)
-g 	genomic DNA file
-o	exonerate output file
-p	parsed exonerate fasta file
-b	best hit only (y/n)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
" ;


my ($file, $output, $genomic, $parsed, $best);

GetOptions(
        'f|fasta:s'     => \$file,
        'o|output:s'     => \$output,	
	'g|genomic:s'	=> \$genomic,
	'p|parsed:s'	=> \$parsed,
	'b|best:s'	=> \$best,	
         
	 );


if( ! defined $file) {
print "$usage\nWARNING: Cannot proceed without fasta file to process\n\n"; exit;
}
if( ! defined $output) {
print "$usage\nWARNING: Cannot proceed without exonerate output file\n\n"; exit;
}
if( ! defined $genomic) {
print "$usage\nWARNING: Cannot proceed without genomic DNA file\n\n"; exit;
}
if( ! defined $parsed) {
print "$usage\nWARNING: Cannot proceed without parsed fasta output file name\n\n"; exit;
}

if( ! defined $best)
        {
        print "\n$usage\nWARNING: Cannot proceed without best hit only y/n \n\n"; exit;
        }

####################################################

### RUN EXONERATE
if ($best eq "Y" || $best eq "y")
{
system "~/tools/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --ryo \">%ti (%tab - %tae)\n%tcs\n\" --showcigar no --showquerygff no --showvulgar no --showalignment yes --useaatla no $file $genomic > $output";
}

elsif ($best eq "N" || $best eq "n")
{
system "~/tools/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --ryo \">%ti (%tab - %tae)\n%tcs\n\" --showcigar no --showquerygff no --showvulgar no --showalignment yes --useaatla no $file $genomic > $output";
}

else { print "\n$usage\nWARNING: Cannot proceed without best hit only y/n \n\n"; exit;}

#~ PARSE OUTPUT
#####

open FILE, "./$output";
open PARSED, ">$parsed";
open SUMMARY, ">$parsed\.summary";
print SUMMARY "Parsed From\t$output\nGenomic DNA\t$genomic\nProtein Seq\t$file\nquery\tscore\tquery_start\tquery_end\tquery_length\ttarget_start\ttarget_end\tStrand\n";

my $query = "N_A";
my $score = "N_A";
my $t_end = "";
my $t_start = "";
my $q_end = "";
my $q_start = "";
my $strand = "";
my $query_length = "";


while (<FILE>)
{
chomp $_;
if ($_ =~ /^C4/){}
elsif ($_ =~ /^----/){}
elsif ($_ =~ /^-- completed exonerate analysis/) {}
elsif ($_ =~ /^Command line\:/){}
elsif ($_ =~ /^Hostname/){}


elsif ($_ =~ /^\s+Query\:\s+(.*?)$/){$query = $1}
elsif ($_ =~ /^\s+Raw score\:\s+(\d+)/) {$score = $1}
elsif ($_ =~ /^\s+Query range\:\s+(\d+)\s+\-\>\s+(\d+)\s*/)
		{
		$q_start = $1;
		$q_end = $2;
		$query_length = ($q_end - $q_start);
		#~ print "$_\n$q_start\n$q_end\n";exit;
		}
elsif ($_ =~ /^\s+Target range\:\s+(\d+)\s+\-\>\s+(\d+)\s*/)
		{
		$t_start = $1;
		$t_end = $2;
		if
		($t_end < $t_start)
			{
			$strand = "-";
			my $left = $t_end;
			my $right = $t_start;
			$t_end = $right;
			$t_start = $left;
			}
		else 
		{
		$strand = "+";
		}
		
		}

elsif ($_ =~ /^\>/)
	{
	print PARSED"$_\t$query\t$score\n"; 
	print SUMMARY "$query\t$score\t";
	print SUMMARY "$q_start\t$q_end\t$query_length\t$t_start\t$t_end\t$strand\n";
	$query = "N_A"; $score = "N_A";
	}
elsif ($_ =~ /^[A-Za-z0-9\-]/)
	{
	print PARSED "$_\n";
	
	}

else {}
}
