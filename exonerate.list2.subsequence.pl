#!/usr/bin/perl 
use warnings;
use strict;


use Getopt::Long;


# Richard Emes University of Nottingham 2010
my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
(C) R D Emes University of Nottingham 2010

pull multiple subsequences from multiple fasta files

USAGE:
-f	path to directory containing fasta files of type name.fa
-c	co-ordinates file name, chrosome, start, end, strand (-/+)
-o	output file
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
" ;


my ($file, $coord, $out);

GetOptions(
        'f|fasta:s'     => \$file,
        'c|coorrd:s'   => \$coord,
	'o|output:s'   => \$out,	
                 );


if( ! defined $file) {
print "$usage\nWARNING: Cannot proceed without directory of fasta files to process\n\n"; exit;
}
if( ! defined $coord) {
print "$usage\nWARNING: Cannot proceed coordinates file\n\n"; exit;
}
if( ! defined $out) {
print "$usage\nWARNING: Cannot proceed without output file\n\n"; exit;
}
####################################################

if ($file =~ /(.*?)\/$/){$file = $1;}
open OUT, ">$out";

# read in and sort coordinates file into chrosomes
open COORD, $coord;
my @list_of_fasta = (); # will hold a unique list of fasta files to process
my %fasta_lookup;

while (<COORD>)
{
chomp $_;
my $line = $_;
my @data = split '\t', $line;
my $fasta_file = $data[1];
push @{$fasta_lookup{$fasta_file}}, $line;
push @list_of_fasta, $fasta_file;
}

my @uniq_list = uniq_array(@list_of_fasta);
@list_of_fasta = @uniq_list;

foreach (@list_of_fasta)
{
chomp $_;
my $fasta = $_;
my @working = @{$fasta_lookup{$fasta}};

## read in fasta file 
my $fasta_sequence;
{
	open FASTA, "$file\/$fasta";
	{
	local $/ = '>'; 
	<FASTA>;                                             # throw away the first line 'cos will only contain ">"

	while (<FASTA>) 
		{	
		chomp $_;
		my ($seq_id, @sequence) = split "\n";            # split the fasta input into Id and sequence
		$fasta_sequence = join '',@sequence;          # reassembles the sequence
		@sequence = ();
		}
	}
	close FASTA;
}

foreach (@working) 
{
chomp $_;
my @data2 = split '\t', $_;
my $name = $data2[0];
my $chr = $data2[1];
my $start = $data2[2];
my $start_print = $start;
#~ $start--;
my $end = $data2[3];
my $strand = $data2[4];

my $length = ($end-$start);

my $seq = substr($fasta_sequence, $start, $length); #$seq, start, length of substring)

if ($strand eq "-")
	{
	my $rev_comp = reverse $seq;
	$rev_comp =~ tr/NACGTacgtn/NTGCAtgcan/;
	$seq = $rev_comp;
	}


print OUT "\>$name\_$chr\_$start_print\_$end\n$seq\n";

}


}

close OUT;
exit;

##### SUBROUTINES
sub uniq_array{
##### make a unique list from the @genes array
my @in = @_;
my %seen = ();
my @uniq = ();
foreach (@in)
{chomp $_;
unless ($seen{$_}){
$seen{$_} = 1;
push (@uniq, $_);
}
}

return @uniq;
}
