#!/bin/perl
#17 July 2015 - Kuan-Lin Huang @ WashU - 
# find druggable mutations in a MAF

use strict;
use warnings;
 
my $usage = 'perl find_dMut.pl <input1>  
';
 
die $usage , unless @ARGV == 1;
my ( $input1 ) = @ARGV;
my $gene2dmut = parsef("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/1_unionOfGeneVariantAADrug.tsv_hugoified");

open ( IN , "<$input1" ) or die "Cannot open $input1: $!";
while ( <IN> )
{
	chomp;
	my @F = split "\t" , $_;
	my $sample = $F[0];
	my $gene = $F[7];
	my $mut = $F[16];
	if (exists($gene2dmut ->{$gene}->{$mut})){
		print "$sample\t$gene $mut\n";
	}
 
}
close IN;

sub parsef {
	my $f = shift;
	my $gene2dmut = {};
	open(FILE, $f) or die "Unable to open file $f due to $!";
	while(<FILE>) {
                chomp;
		my @F = split "\t" , $_;
		$gene2dmut->{$F[0]}->{$F[2]}=1;
        }
        close FILE;
	return $gene2dmut;
}
