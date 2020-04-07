#!/usr/bin/perl
use strict;
use warnings;



# declaring some variables
my ($dna_seq,$revdna_seq,$seqlen,$TGA,$TAG,$TAA,$TGAr,$TAGr,$TAAr,$Ntuples,$Nrevtuples,$revORFcounter,@revORFlengths,$ORFcounter,@ORFlengths);

# initializing variables
$Ntuples= $Nrevtuples = $revORFcounter = $ORFcounter = 0;
@revORFlengths=();
@ORFlengths=();

my $file = shift @ARGV; # get input filename from
                        # script command-line arguments

# Open input file
open(BED, "< $file");

# Loading input file records

while (<BED>) {

    next if $_ =~ m/^\>/o; # skipping comment lines starting with ">"

    chomp;    # Remove trailing newline char

    $dna_seq .= $_;

    };# while

#Getting reverse strain
$revdna_seq = reverse($dna_seq);
$revdna_seq =~  tr/ACGT/TGCA/;
$seqlen = length($dna_seq);

loop_str(0);
loop_str(1);
loop_str(2);

loop_rev(0);
loop_rev(1);
loop_rev(2);

# Printing results for straigh strain
print STDOUT "Forward strain total TGA = $TGA\n";
print STDOUT "Forward strain total TAG = $TAG\n";
print STDOUT "Forward strain total TAA = $TAA\n\n";

# Printing results for reverse strain
print STDOUT "Reverse strain total TGArev = $TGAr\n";
print STDOUT "Reverse strain total TAGrev = $TAGr\n";
print STDOUT "Reverse strain total TAArev = $TAAr\n\n";


my $TGAtotal = $TGA + $TGAr;
my $TAGtotal = $TAG + $TAGr;
my $TAAtotal = $TAA + $TAAr;
my $totalstopcodons = $TGAtotal + $TAGtotal + $TAAtotal;
my $totaltuples = 2*$seqlen/3*3; # 2xseqlen because of the forward and reverse strains, /3 because of the codons, x3 because of reading frame
				

#Getting the frequencies
my $TGAfreqall = $TGAtotal/$totaltuples*100;
my $TAGfreqall = $TAGtotal/$totaltuples*100;
my $TAAfreqall = $TAAtotal/$totaltuples*100;

my $TGAfreqstop = $TGAtotal/$totalstopcodons*100;
my $TAGfreqstop = $TAGtotal/$totalstopcodons*100;
my $TAAfreqstop = $TAAtotal/$totalstopcodons*100;


print "CODON\t\tRelFreqAll\tRelFreqStp\n";
print "TGA\t\t",$TGAfreqall,"%\t\t",$TGAfreqstop,"%\n";
print "TAG\t\t",$TAGfreqall,"%\t\t",$TAGfreqstop,"%\n";
print "TAA\t\t",$TAAfreqall,"%\t\t",$TAAfreqstop,"%\n";

exit(0);

## FUNCTIONS

sub loop_str {
   my $frame = shift @_;
   for (my $i = $frame; $i < $seqlen; $i+=3) { # i+=3 because we are reading codons
    	my $codon = uc(substr($dna_seq,$i,3));
    	$Ntuples++;
    	SWITCH: {
		$codon eq 'TGA' && ($TGA++, push (@ORFlengths,$ORFcounter), $ORFcounter = 0, last SWITCH);
		$codon eq 'TAG' && ($TAG++,push (@ORFlengths,$ORFcounter), $ORFcounter = 0, last SWITCH); 
		$codon eq 'TAA' && ($TAA++,push (@ORFlengths,$ORFcounter), $ORFcounter = 0, last SWITCH); 
		$ORFcounter+=3;# 3 nucleotides per codon
	};
   };
   print STDOUT "Forward strain ORF lengths for reading frame ", $frame, ": \n", join(",", @ORFlengths), "\n\n";
   @ORFlengths=();
   return (@ORFlengths);
   $ORFcounter = 0;  
};


sub loop_rev() {
    my $frame = shift @_;
    for (my $i = $frame; $i < $seqlen; $i+=3) { 
    	my $revcodon = uc(substr($revdna_seq,$i,3));
    	$Nrevtuples++; 
    	SWITCH: {
		$revcodon eq 'TGA' && ($TGAr++, push (@revORFlengths,$revORFcounter), $revORFcounter = 0, last SWITCH);
		$revcodon eq 'TAG' && ($TAGr++, push (@revORFlengths,$revORFcounter), $revORFcounter = 0,last SWITCH); 
		$revcodon eq 'TAA' && ($TAAr++, push (@revORFlengths,$revORFcounter), $revORFcounter = 0,last SWITCH);
		$revORFcounter+=3;
	};
    };
    print STDOUT "Reverse strain ORF lengths for reading frame ", $frame, ": \n", join(",", @revORFlengths), "\n\n";
    @revORFlengths=();
    return(@revORFlengths);
    $revORFcounter = 0;  
};


