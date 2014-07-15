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
-s	protein sequence file to use as template for predictions [fasta format]
-d	directory of of DNA or chromosome files [should be named x.fa]
-p	prefix for gene predictions
-o	minimum size of overlap (amino acids of protein seq)
-i	maximum number of introns allowed in predicted genes
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
" ;


my ($seq, $dir, $prefix, $overlap, $max_introns);

GetOptions(
        's|seq:s'     => \$seq,
	'd|dir:s'     => \$dir,
        'p|prefix:s'     => \$prefix,
        'o|overlap:s'     => \$overlap,
	'i|prefix:s'     => \$max_introns,
	 );


if( ! defined $seq) {
print "$usage\nWARNING: Cannot proceed without peptide sequence file to process\n\n"; exit;
}
if( ! defined $dir) {
print "$usage\nWARNING: Cannot proceed without directory of chromosomes to search\n\n"; exit;
}
if( ! defined $prefix) {
print "$usage\nWARNING: Cannot proceed without prefix name\n\n"; exit;
}

if( ! defined $overlap) {
print "$usage\nWARNING: Cannot proceed without minimum overlap\n\n"; exit;
}
if( ! defined $max_introns) {
print "$usage\nWARNING: Cannot proceed without max number of introns\n\n"; exit;
}



####################################################



#collect all data folders to be analyzed
if ($dir =~ /\/$/){}
else {$dir = $dir."/";}


opendir(FILES, $dir)  || die "can't open directory ($!)";

my @tmp = sort readdir(FILES);
my @folders = ();
closedir(FILES);
foreach(@tmp)
        {
        if(/\.fa$/)
                {push(@folders, $_);
                }
        }
@tmp=();


my $time = scalar localtime;
print "\nStarting Exonerate\t$time\n";


foreach (@folders)
	{
	chomp $_;
	my $filename = $_;
	my $file2process = "$dir"."$filename";
	system "perl exonerate.predict.genes.pl -f $seq -g $file2process -o $seq\.$filename\.ex -p $seq\.$filename\.ex\.fa -b N";
	print ".";
	}
print "\n";
$time = scalar localtime;
print "Finished Exonerate\t$time\n";


$time = scalar localtime;
print "\nParsing Exonerate\t$time\n";

foreach (@folders)
	{
	chomp $_;
	my $filename = $_;
	my $infile = "$seq\.$filename"."\.ex";
	system "perl exonerate.parse.pl -f $infile -o $infile\.parsed";
	print ".";
	}
print "\n";
$time = scalar localtime;
print "Finished Parsing Exonerate\t$time\n";




$time = scalar localtime;
print "\nFinding non-overlapping best hits\t$time\n";

foreach (@folders)
	{
	chomp $_;
	my $filename = $_;
	my $infile = "$seq\.$filename\.ex\.parsed";
	system "perl find.non.overlapping.exonerate.predictions.pl -f $infile -o $infile\.nonoverlapping -p $prefix -c $filename -i $max_introns -m $overlap";
	print ".";
	}
print "\n";
$time = scalar localtime;
print "Finished Finding non-overlapping predicted genes\t$time\n";


$time = scalar localtime;
print "\nParse hit list\t$time\n";
my $uniq_number = 1;
open TABLE, ">$seq\.predictions.table";
print TABLE "Unique_Name\tChr_prediction_name\tChr\tProtein Seq\tscore\tquery_start\tquery_end\tquery_length\ttarget_start\ttarget_end\ttarget_length\tStrand\tIntrons\n";
foreach (@folders)
	{
	chomp $_;
	my $filename = $_;
	my $infile = "$seq\.$filename\.ex\.parsed\.nonoverlapping";
	open FILE, "<$infile";
		while (<FILE>)
		{
		chomp $_;
		my $line = $_;
		if ($line =~ /^Name/){}
		else {
			print TABLE "$prefix\.$uniq_number\t$line\n";
			$uniq_number++;
			}
		}
	close FILE;
	}

$time = scalar localtime;
print "Finished Producing hit list table\t$time\n";
close TABLE;

$time = scalar localtime;
print "\nGenerate Fasta File of predictions\t$time\n";
open POSITIONS_TABLE, ">tmp.table";
open INFILE, "$seq\.predictions.table";
<INFILE>; #remove header

while (<INFILE>)
	{
	chomp $_;
	my $line = $_;
	my @data = split '\t', $line;
	print POSITIONS_TABLE "$data[0]\t$data[2]\t$data[8]\t$data[9]\t$data[11]\n";
	}
	close POSITIONS_TABLE;

system "split -l 50 tmp.table split.";

opendir(SPLITS, "./")  || die "can't open directory ($!)";
my @splits = sort readdir(SPLITS);
my @split_files = ();

closedir(SPLITS);
foreach(@splits)
        {
        if(/split\.*/)
                {push(@split_files, $_);
                }
        }
@splits = ();

foreach (@split_files)
{
chomp $_;

system "perl list2.subsequence.pl -f $dir -c $_ -o $_\.splits.tojoin";
}

system "cat *.splits.tojoin > $seq\.predictions.fa";
system "rm split*";


$time = scalar localtime;
print "Finished Generate Fasta File of predictions\t$time\n";
close INFILE;


## Cleanup
system "rm tmp.table";
system "mkdir Exonerate_prediction_files";
system "mv *.ex *.ex.fa *.summary *.parsed *.parsed.nonoverlapping Exonerate_prediction_files/";
exit;
