#!/bin/perl
#27 July 2015 - Kuan-Lin Huang @ WashU -
# outlier_panOmics.pl: from an outlier table with one column for one omic level, one row for one sample
# get the same gene's expression on other levels 
 
use strict;
use warnings;
use Time::Piece;

my $date = localtime->strftime('%Y-%m-%d');
my $outfile_all = "results/".$date."_all_outliers.txt";
my $outfile_filter = "results/".$date."_filtered_outliers.txt";

my $usage = 'perl summarize_outlier.pl  
';
 
die $usage , unless @ARGV == 0;

# dependencies
my $base_path = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/figures/";
my $fDate = "2015-10-09/2015-10-09";
#my $drug_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/reference_files/gene_drug_list/UnionOfVariantDrug_babyHugo_genedrug_TARGET_db_v3_merged_hugoified_curated_noKns_noTSG.list.txt";

my @cancers = ("BRCA", "CRC", "OV PNNL", "OV JHU");

my $outliers = {}; # ex: datatype -> gene -> outlier; 
#my $gene2drug = {}; # ex. gene -> drug;
my $samples = {}; # ex: datatype -> sample -> has_outlier;
my $data2n = {};
print "Cancer\tLevel\tGene\tNum_outlier\tNum_samples\tOutlier_precentage\n";

foreach my $cancer (@cancers){
# for BRCA
my $CNV_f = "_KH_".$cancer."\ druggable\ CNV_outlier_score_table.txt";
my $RNA_f = "_KH_".$cancer."\ druggable\ RNA_outlier_score_table.txt";
my $Pro_f = "_KH_".$cancer."\ druggable\ proteome_outlier_score_table.txt";
my $Pho_f = "_KH_".$cancer."\ druggable\ phosphoproteome_outlier_score_table.txt";

if ($cancer eq "OV PNNL" || $cancer eq "OV JHU"){
	$CNV_f = "_KH_OV\ druggable\ CNV_outlier_score_table.txt";
	$RNA_f = "_KH_OV\ druggable\ RNA_outlier_score_table.txt";
}
# global variables

$outliers = {}; # ex: datatype -> gene -> outlier; 
#my $gene2drug = {}; # ex. gene -> drug;
$samples = {}; # ex: datatype -> sample -> has_outlier;
$data2n = {}; # datatype -> sample size

# read and store outlier results
read_score($base_path.$fDate.$CNV_f, "CNV"); #print STDERR $base_path.$fDate.$CNV_f."\n";
read_score($base_path.$fDate.$RNA_f, "RNA");
read_score($base_path.$fDate.$Pro_f, "PRO");

# bonus levels
if ($cancer eq "BRCA" || $cancer eq "OV PNNL"){
	read_score($base_path.$fDate.$Pho_f, "PHO");
}
if ($cancer eq "OV JHU"){
	my $Gly_f = "_KH_".$cancer."\ druggable\ glycoproteome_outlier_score_table.txt";
	read_score($base_path.$fDate.$Gly_f, "GLY");
}

#read_drug($drug_f);

#open ( my $all_out, ">$outfile_all" ) or die "Cannot open $outfile_all: $!";
# print all data levels except phospho

#print $cancer."\n";
my %outliers = %$outliers;
my %data2n = %$data2n;
my %samples = %$samples;
my $dtype_outlier_num;
foreach my $dtype (sort keys %outliers){
	$dtype_outlier_num = 0;
	foreach my $sample_wOutliers (keys %{$samples{$dtype}}){
		$dtype_outlier_num++; #print "Test\t".$sample."\t".$samples{$dtype}{$sample}."\n";}
	}
	my $dtype_num = $data2n{$dtype};
	my $outlier_perA = 100*$dtype_outlier_num/$dtype_num;
	print $cancer."\t".$dtype."\tSum\t".$dtype_outlier_num."\t".$dtype_num."\t".$outlier_perA."\n";
	foreach my $gene (sort keys %{$outliers{$dtype}}){		
		my $outlier_num = $outliers{$dtype} -> {$gene};
		my $outlier_per = 100*$outlier_num/$dtype_num;
		print $cancer."\t".$dtype."\t".$gene."\t".$outlier_num."\t".$dtype_num."\t".$outlier_per."\n";
	}
}

sub read_score {
	my $dFile = shift;
	my $type = shift;
	open ( IN , "<$dFile" ) or die "Cannot open $dFile: $!";
	my @G = split "\t", <IN>; 
	my $dsize = $#G-1;
	$data2n -> {$type} = $dsize;
	while (<IN>)
	{
		chomp;
		my @F = split "\t", $_;
		my $gene = $F[0];
		my $outlier_num = 0;

		for (my $i  = 1; $i < $#F; $i++){
			if ($F[$i] eq "NA"){next;}
			elsif ($F[$i] >= 1.5){ $outlier_num++; $samples -> {$type} -> {$G[$i]} = 1}
		} 
		$outliers -> {$type} -> {$gene} = $outlier_num;
	} 
	close (IN);
}

}

#close $all_out;

# open ( my $f_out, ">$outfile_filter" ) or die "Cannot open $outfile_filter: $!";
# # print all data levels except phospho
# print $f_out "sample\tCNV_outlier\tmRNA_outlier\tLFQ_protein_outlier\tPro_protein_outlier\tPhosphosite_outlier\n";
# foreach my $sample (@samples){
# 	print $f_out "$sample\t";
# 	foreach my $data (@data){
# 		if ($data eq "CNV" || $data eq "RNA"){
# 			if (exists($outliers -> {$sample} -> {$data})){
# 				my @outliers = @{$outliers -> {$sample} -> {$data}};
# 				OUTER:
# 				my $outlierP={};
# 				foreach my $outlier (@outliers){
# 					foreach my $data (@data){
# 						if ($data eq "CNV" || $data eq "RNA") {next;}
# 						if (exists($gene2score -> {$sample} -> {$data} -> {$outlier})){
# 							my $score = $gene2score -> {$sample} -> {$data} -> {$outlier};
# 							if (exists($outlierP->{$outlier})){next;}
# 							if ($score >= 1.0) {print $f_out "$outlier"."; "; $outlierP->{$outlier}=1;} #last OUTER;} # print and break out of the outlier loop
# 						}
# 					}
# 				}
# 			} #else {print $f_out "NA";} #data type not available
# 		} elsif ($data eq "PHO"){
# 			if (exists($outliers -> {$sample} -> {$data})){
			
# 				my @outliers = @{$outliers -> {$sample} -> {$data}};
# 				if ($#outliers > 3){@outliers = @outliers[0..2];}
# 				foreach my $outlier (@outliers){
# 					my $h_val = 0;
# 					if (exists($outlier_num_h -> {"PHO_H"} -> {$outlier})){$h_val = $outlier_num_h -> {"PHO_H"} -> {$outlier};} 
# 					print $f_out "$outlier($h_val)"."; ";
# 				}
# 			} #else {print $f_out "NA";} #data type not available
# 		} else {
# 			if (exists($outliers -> {$sample} -> {$data})){
			
# 				my @outliers = @{$outliers -> {$sample} -> {$data}};
# 				foreach my $outlier (@outliers){
# 					my $h_val = 0;
# 					if (exists($outlier_num_h -> {"PRO_H"} -> {$outlier})){$h_val = $outlier_num_h -> {"PRO_H"} -> {$outlier};} 
# 					print $f_out "$outlier($h_val)"."; ";
# 				}
# 			} #else {print $f_out "NA";} #data type not available
# 		}
# 		print $f_out "\t";
# 	}
# 	print $f_out "\n";
# }
# close $f_out;

# functions
# sub read_outlier {
# 	my $dFile = shift;
# 	my $type = shift;
# 	open ( IN , "<$dFile" ) or die "Cannot open $dFile: $!";
# 	<IN>; # skip header
# 	while (<IN>)
# 	{
# 		chomp;
# 		my @F = split "\t", $_;
# 		$data2n -> {$type} = $#F;
# 		# first line: read in sample and gene and its rank
# 		my $sample = $F[0];
# 		my $rank2gene = {};
# 		for (my $i  = 1; $i < $#F; $i++){
# 			$rank2gene -> {$i} = $F[$i];
# 		} 

# 		# second line: retrieve outlier boolean for each gene
# 		my $line = <IN>;
# 		my @L = split "\t", $line;
# 		for (my $i  = 1; $i < $#L; $i++){
# 			my $gene = $rank2gene -> {$i};
# 			if ($L[$i] eq "TRUE"){
# 				if (exists($outliers -> {$sample} -> {$type})){ push $outliers -> {$sample} -> {$type}, $gene;} else {$outliers -> {$sample} -> {$type} = [$gene];}
# 			}
# 		} 
# 	} 
# 	close (IN);
# }

# sub read_score {
# 	my $dFile = shift;
# 	my $type = shift;
# 	open ( IN , "<$dFile" ) or die "Cannot open $dFile: $!";
# 	my @F = split "\t", <IN>; 
# 	my $dsize = $#F-1;
# 	$data2n -> {$type} = $dsize;
# 	while (<IN>)
# 	{
# 		chomp;
# 		my @F = split "\t", $_;
# 		my $gene = $F[0];
# 		my $outlier_num = 0;

# 		for (my $i  = 1; $i < $#F; $i++){
# 			if ($F[$i] eq "NA"){next;}
# 			elsif ($F[$i] >= 1.5){ $outlier_num++;}
# 		} 
# 		$outliers -> {$type} -> {$gene} = $outlier_num;
# 	} 
# 	close (IN);
# }

# sub read_drug {
# 	my $file = shift; 
# 	open ( IN , "<$file" ) or die "Cannot open $file: $!";
# 	<IN>; # skip header
# 	while (<IN>)
# 	{
# 		chomp;
# 		my @F = split "\t", $_;
# 		my $gene = $F[0];
# 		my $drug = $F[3];
# 		$gene2drug -> {$gene} -> {$drug} = 1;
# 	} 
# 	close (IN);
# }
