#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;


# Richard Emes University of Nottingham 2010
my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
R D Emes University of Nottingham 2010
K Brown University of Nottingham 2014

predict genes from target chromosomes that could encode proteins of interest.
Spawns multiple jobs (will take time be ready).

USAGE:
-s  protein sequence file to use as template for predictions [single or multiple fasta format]
-t  target chromosome or contig files [multiple fasta format]
-p  prefix for gene predictions
-o  minimum size of overlap (amino acids of protein seq)
-i  maximum number of introns allowed in predicted genes
-j  max number of exonerate jobs to spawn

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
" ;


my ($seq, $target_fa, $prefix, $overlap, $max_introns,$jobs);

GetOptions(
        's|seq:s'     => \$seq,
	       't|tar:s'     => \$target_fa,
        'p|prefix:s'     => \$prefix,
        'o|overlap:s'     => \$overlap,
	       'i|introns:s'     => \$max_introns,
	 	      'j|jobs:s'     => \$jobs,
    );


if( ! defined $seq) {
print "$usage\nWARNING: Cannot proceed without peptide sequence file to process\n\n"; exit;
}
if( ! defined $target_fa) {
print "$usage\nWARNING: Cannot proceed without target fasta to search\n\n"; exit;
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
if( ! defined $jobs) {
print "$usage\nWARNING: Cannot proceed without max number of jobs\n\n"; exit;
}


####################################################


#
# #filter contigs if too short to hold sequence of interest
# # generate lookup so can parse files easier if long names.
# open SEQLOOKUP, ">$prefix.Target.sequence.lookup.txt";
# open PASSEDTARGETFA, ">$prefix.Target.passed.filter.fa";
# open TMPTABTARGETFA, ">$prefix.Target.passed.filter.tab.tmp";
#
# my $time = scalar localtime;
# print "\n$time Filtering target fasta\n";
# open FASTA, $target_fa;
# my $seq_tally = 1;
# my $passed_count = 0;
# {
# local $/ = '>';
# <FASTA>;                                  # throw away the first line 'cos will only contain ">"
#
# while (<FASTA>) {
#     chomp $_;
#     my ($seq_id, @sequence) = split "\n";            # split the fasta input into Id and sequence
#     my $sequence = join '',@sequence;          # reassembles the sequence
#         my $length = length $sequence;
#
#         if ($length >= 3*$overlap)
#         {
#         print PASSEDTARGETFA "\>$seq_tally\n$sequence\n";
#         print SEQLOOKUP "$seq_tally\t$seq_id\n";
#         print TMPTABTARGETFA "$seq_tally\t$sequence\n";
#         $passed_count++;
#       }
#
#
# $seq_tally++;
# }
# }
#
#
# ###################################################################################
# ### SPLIT CONTIG FASTA AND RUN EXONERATE
# ###################################################################################
#
#
# $time = scalar localtime;
# print "$time Splitting files\n";
#
# # split fasta into multiple files to run exonerate on
# my $maxlines = int(($passed_count/$jobs)+0.5);
# system "split -l $maxlines $prefix.Target.passed.filter.tab.tmp splits --additional-suffix=.tmp";
#
# opendir(TMPFILES, "./");
# my @tmp = sort readdir(TMPFILES);
# closedir(TMPFILES);
#
#
# my @fastain = ();
# my $fasta_count = 1;
# foreach(@tmp)
#         {
#         chomp $_;
#         if(/^splits.*\.tmp$/)
#                 {
#                 chomp $_;
#                 my $filenamein = $_;
#                 open FILEIN, $filenamein;
#                 open FILEOUT, ">$filenamein\.fa.tmp";
#                 while (<FILEIN>){
#                 chomp $_;
#                 if ($_ =~ /(.*)\t(.*)/)
#                 {
#                 print FILEOUT "\>$1\n$2\n";
#                 }
#                 }
#
#                 push(@fastain, "$filenamein\.fa.tmp");
#                 }
#         }
# @tmp=();
#
# $time = scalar localtime;
# print "$time Starting Exonerate using nohup to spawn multiple jobs\n";
#
# my @tmpexonerate = ();
# foreach (@fastain)
# 	{
# 	chomp $_;
# 	my $filename = $_;
#   my $output = "$filename\.ex\.tmp";
#   push (@tmpexonerate, $output);
#   system "nohup exonerate --model protein2genome --ryo \">%ti (%tab - %tae)\n%tcs\n\" --showcigar no --showquerygff no --showvulgar no --showalignment yes --useaatla no $seq $filename > $output &";
#   }
#
# my $finished_flag = 0;
# my $exonerate_complete = 0;
#
# while ($finished_flag == 0)
#   {
#   system "tail -n1 *.ex.tmp > tails.tmp";
#   open TAILSTMP, "<tails.tmp";
#   while (<TAILSTMP>)
#     {
#     chomp $_;
#     if ($_ =~ /^-- completed exonerate analysis/) {$exonerate_complete++;}
#     }
#   close TAILSTMP;
#
#   if ($exonerate_complete < $jobs) {sleep 60; $exonerate_complete=0;}
#   elsif ($exonerate_complete == $jobs)
#         {
#         $finished_flag = 1;
#         print "\n";
#         $time = scalar localtime;
#         print "$time Finished Exonerate\n";
#         }
#   }
#
# # join exonerate files
# foreach (@tmpexonerate)
#   {
#   chomp $_;
#   system "cat $_ >> $prefix.ex";
#   }
# # clean up a bit
# system "rm ./*.tmp";
#
# ###################################################################################
# ### PARSE EXOENERATE FILES
# ###################################################################################
# # parse exonerate files
# $time = scalar localtime;
# print "$time Parsing Exonerate\n";
#
# open FILE, "$prefix.ex";
# open OUT, ">$prefix.ex.parsed";
# print OUT "Query\tTarget\tscore\tquery_start\tquery_end\tquery_length\ttarget_start\ttarget_end\ttarget_length\tStrand\tIntrons\n";
#
# my $query = "N_A";
# my $target = "N_A";
# my $score = "N_A";
# my $t_end = "";
# my $t_start = "";
# my $q_end = "";
# my $q_start = "";
# my $strand = "";
# my $query_length = 0;
# my $target_length = 0;
# my $intron_count = 0;
# my $flag_hit = 0;
# while (<FILE>)
# {
# chomp $_;
# if ($_ =~ /^C4/){}
# elsif ($_ =~ /^----/){}
# elsif ($_ =~ /^-- completed exonerate analysis/) {}
# elsif ($_ =~ /^Command line\:/){}
# elsif ($_ =~ /^Hostname/){}
#
#
# elsif ($_ =~ /^\s+Query\:\s+(.*)/){$query = $1;  $query =~ s/\s+$//;}
# elsif ($_ =~ /^\s+Target\:\s+(.*)/){$target = $1; $target =~ s/\s+$//;$target =~ s/\s+\[revcomp\]//;}
# elsif ($_ =~ /^\s+Raw score\:\s+(\d+)/) {$score = $1}
# elsif ($_ =~ /^\s+Query range\:\s+(\d+)\s+\-\>\s+(\d+)\s*/)
# 		{
# 		$q_start = $1;
# 		$q_end = $2;
# 		$query_length = ($q_end - $q_start);
# 		}
# elsif ($_ =~ /^\s+Target range\:\s+(\d+)\s+\-\>\s+(\d+)\s*/)
# 		{
# 		$t_start = $1;
# 		$t_end = $2;
# 		$target_length = ($t_end - $t_start);
# 		if
# 		($t_end < $t_start)
# 			{
# 			$strand = "-";
# 			my $left = $t_end;
# 			my $right = $t_start;
# 			$t_end = $right;
# 			$t_start = $left;
# 			$target_length = ($t_end - $t_start);
# 			}
# 		else
# 		{
# 		$strand = "+";
# 		}
#
# 		}
# elsif ($_ =~ /> T/ || /\:\s+Targ/)
# 	{
# 	$intron_count++;
# 	}
#
# elsif ($_ =~ /^\>/)
# 	{
# 	print OUT "$query\t$target\t$score\t";
# 	print OUT "$q_start\t$q_end\t$query_length\t$t_start\t$t_end\t$target_length\t$strand\t$intron_count\n";
# 	$query = "N_A"; $score = "N_A";
# 	$intron_count = 0;
# 	}
#
#
# else {}
# }
#
#
#
# $time = scalar localtime;
# print "$time Finished Parsing Exonerate\n";
# # #
###################################################################################
### MERGE OVERLAPPING AND FIND NON OVERLAPPING BEST HITS BASED ON EXOENERATE SCORES
###################################################################################

my $time = scalar localtime;
print "$time Finding non-overlapping best hits\n";

open PARSED, "<$prefix.ex.parsed";
my $header = <PARSED>; # drop header
chomp $header;
# quick filter for introns and overlap.
# generate lookup

my @contigswithhits = ();
my %contiglookup;

while (<PARSED>)
  {
    chomp $_;
    my $line = $_;
    my @data = split '\t', $line;
    my $contig = $data[1];
    if ($data[5] =~ /^\d/ && $data[10] =~ /^\d/)
      {
      if ($data[5] >= $overlap && $data[10] <= $max_introns)
        {
        push (@contigswithhits, $contig);
        push(@{$contiglookup{$contig}},$line);
        }
      }
    }

my @uniqcontigswithhits = uniq_array(@contigswithhits);

foreach (@uniqcontigswithhits)
  {
  chomp $_;
  my $contig = $_;
  open CONTIGTMP, ">$contig\.ex\.tmp";

  foreach (@{$contiglookup{$contig}})
        {
        chomp $_;
        print CONTIGTMP "$_\n";
        }
  close CONTIGTMP;


  system "sort -n -k7 $contig\.ex\.tmp > ./hold.tmp"; ## sort numerically on the start position of prediction



   open FILE, "<hold.tmp";
   open OUT, ">$contig\.ex\.nonoverlapping.tmp";

   my $flag = 0;
   my $gene_start = 0;
   my $gene_end = 0;
   my @holding_gene_details = ();

   while (<FILE>)
   {
   	chomp $_;
   	my $line = $_;
   	my @data = split '\t', $line;
   	my $start = $data[6];
  	my $end = $data[7];
  	my $nintrons = $data[10];
  	my $length = $data[5];
  	if (($length >= $overlap) &&($nintrons <= $max_introns))
  		{
  		if ($flag ==0)
  			{
  			$gene_start = $start;
  			$gene_end = $end;
  			$flag++;
  			}

  		if ($start >= $gene_start && $start < $gene_end)
  			{
  			push @holding_gene_details, $line; ### these are the overlapping
  			if ($end > $gene_end) {$gene_end = $end;} # keeps tabs on extreme right side #################
  			}
  		elsif ($start >= $gene_end)
  			{
  			my $current_best = 0;
  			my $best_hit_line = "";

  			foreach (@holding_gene_details)
  					{
  					chomp $_;
  					my $best_line = $_;
  					my @best_data = split  '\t', $best_line;
  					my $best_score = $best_data[2];
  					my $best_intron_count = $best_data[10];
  					my $best_q_length = $best_data[5];
  					if ($best_score >= $current_best && $best_intron_count <= $max_introns && $best_q_length >= $overlap)
  						{
  						$current_best = $best_score;
  						$best_hit_line = $best_line;
  						}
  					}
  			#~ find best in @holding_gene_details

  			if ($best_hit_line =~ /^[A-Za-z0-9\_\-]/)
  				{
  				print OUT "$best_hit_line\n";
  				}
  			@holding_gene_details = ();
  			push @holding_gene_details, $line;
  			$gene_start = $start;
  			$gene_end = $end;
  			}
  	}
  }
  my $current_best = 0;
  my $best_hit_line = "";

  foreach (@holding_gene_details)
  		{
  		chomp $_;
  		my $best_line = $_;
  		my @best_data = split  '\t', $best_line;
  		my $best_score = $best_data[2];
  		my $best_intron_count = $best_data[10];
  		my $best_q_length = $best_data[5];
  		if ($best_score >= $current_best && $best_intron_count <= $max_introns && $best_q_length >= $overlap)
  			{
  			$current_best = $best_score;
  			$best_hit_line = $best_line;
  			}
  		}
  #~ find best in @holding_gene_details

  if ($best_hit_line =~ /^[A-Za-z0-9\_\-]/)
  	{
  	print OUT "$best_hit_line\n";
  	}


  close OUT;
  close FILE;

  system "rm -rf ./hold.tmp nohup.out";
}

###################################################################################
### MERGE BEST HITS BASED ON EXOENERATE SCORES INTO SINGLE FILE AND GIVE UNIQUE ID
###################################################################################

system "cat *\.ex\.nonoverlapping\.tmp > combined.non.overlapping.tmp";

open COMBBEST, "combined.non.overlapping.tmp";
open FINALBEST, ">$prefix\.ex\.parsed\.nonoverlapping";
print FINALBEST "ID\t$header\n";

my @nonoverlapping_lookup = (); # create an array of coordinates for later parsing from fasta file.

my $number = 1;
my %prefixlookup; # keep tabs on named sequence for later additon to fasta file
while (<COMBBEST>)
  {
  chomp $_;
  print FINALBEST "$prefix\.$number\t$_\n";
  my @data = split '\t', $_;
  my $coordinates = "NA";
  if ($data[9] eq "+") {$coordinates = $data[1]." \($data[6] \- $data[7]\)";}
  if ($data[9] eq "-") {$coordinates = $data[1]." \($data[7] \- $data[6]\)";} #order changes for rev strand
  push @nonoverlapping_lookup, $coordinates;
  $prefixlookup{$coordinates} = "$prefix\.$number";
  $number++;
  }


###################################################################################
### PARSE FASTA FROM EXONERATE FILE
###################################################################################

$time = scalar localtime;
print "$time Creating Fasta file of predicted genes\n";

system "grep -v '^ ' $prefix.ex |grep -v '^\$' |grep -v '^-' | grep -v '^C4' | grep -v '^Host' |grep -v '^Command' | grep -v '^%' > $prefix.ex.fasta.tmp";

# turn fasta to tab for quick searching
open ALLFASTA, "$prefix.ex.fasta.tmp";
open ALLFASTATAB, ">$prefix.ex.tab.tmp";
{
local $/ = '>';
<ALLFASTA>;                                             # throw away the first line 'cos will only contain ">"

while (<ALLFASTA>)
        {
        chomp $_;
        my ($seq_id, @sequence) = split "\n";            # split the fasta input into Id and sequence
        my $fasta_sequence = join '',@sequence;          # reassembles the sequence
        print ALLFASTATAB "$seq_id\t$fasta_sequence\n";
      }
}

system "sort $prefix.ex.tab.tmp \| uniq > $prefix.ex.tab.uniq.tmp";

my @nonoverlapping_lookup_uniq = uniq_array(@nonoverlapping_lookup);

foreach (@nonoverlapping_lookup_uniq)
{
chomp $_;
system "grep \"$_\" $prefix.ex.tab.uniq.tmp >> $prefix.ex.tab.uniq.keep.tmp";
}


open KEEPTABIN, "$prefix.ex.tab.uniq.keep.tmp";
open KEEPFASTAOUT, ">$prefix.ex.nonoverlapping.All.predictions.fasta";
while (<KEEPTABIN>){
        chomp $_;
        my $line = $_;
        if ($line =~ /(.*)\t(.*)/)
                {
                print KEEPFASTAOUT "\>$prefixlookup{$1} $1\n$2\n";
                }
        }

# cleanup
system "rm *.tmp";

$time = scalar localtime;
print "$time Finished pipeline\n";




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
