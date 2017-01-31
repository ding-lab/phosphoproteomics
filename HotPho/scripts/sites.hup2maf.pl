#!/usr/bin/perl
#18 November 2016 - Adam Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

use TGI::Mutpro::Preprocess::Uniprot;

my $usage = 'perl phospho2maf.pl <input.hup> <output.maf>
';

die $usage , unless @ARGV == 2;

my ( $input , $output ) = @ARGV;

my $IN = FileHandle->new( $input , "r" );
if ( not defined $IN ) { die "ADSERROR: Could not open/write $input\n"; }

my $OUT = FileHandle->new( $output , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $output\n"; }

$OUT->print( "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\ttranscript_name\tamino_acid_change\n" );

while ( my $line = $IN->getline() ) {
	chomp( $line );
	my ( $hugo , $id ) = (split( /\t/ , $line ))[0,1];
	print $hugo."\n";
	my $uniprot = new TGI::Mutpro::Preprocess::Uniprot( $id );
	my $sites = &getSites( $uniprot );
	#my ( $gene , $species ) = split( /\_/ , &getGeneID( $uniprot ) );
	my $trans = &getCanonicalTranscript( $uniprot );
	my $chromosome = &getChromosome( $uniprot , $sites );
	&printAsMAF( $OUT , $sites , $hugo , $trans , $chromosome );
}
$OUT->close();


sub printAsMAF {
	my $OUT = shift;
	my $sites = shift;
	my $gene = shift;
	my $trans = shift;
	my $chromosome = shift;
	foreach my $site ( sort {$a<=>$b} keys %{$sites} ) {
		my $category = $sites->{$site}->{category};
		my $type = $sites->{$site}->{type};
		my $aa = "X";
		# if ( $type =~ /serine/ ) {
		# 	$aa = "S";
		# } elsif ( $type =~ /tyrosine/ ) {
		# 	$aa = "Y";
		# } elsif ( $type =~ /threonine/ ) {
		# 	$aa = "T";
		# } elsif ( $type =~ /aspar/ ) {
		# 	$aa = "D";
		# } elsif ( $type =~ /histidine/ ) {
		# 	$aa = "H";
		# }

		my $start = int( 1000000000*rand() );
		my $stop = $start;
		$OUT->print( join( "\t" , ( $gene , $chromosome , $start , $stop, "Missense_Mutation" , "A" , "T" , "T" , $category.":".$type , $trans , "p.".$aa.$site."Z" ) )."\n" );
	}
}

sub getChromosome {
	my $this = shift;
	my $sites = shift;
	foreach my $site ( sort {$a<=>$b} keys %{$sites} ) {
		my $proteomes = $this->annotations( "Proteomes" );
		my $chromosome = 1;
		foreach my $line ( @{$proteomes} ) {
			next unless ( $line =~ /Chromosome/ );
			chomp( $line );
			$chromosome = $line;
			$chromosome =~ s/.*Chromosome\s+([A-Z0-9]+)\./$1/g;
			return $chromosome;
		}
	}
}

sub getCanonicalTranscript {
    my $this = shift;
    foreach my $line ( split /\n/ , $this->entireRecord() ) {
        if ( $line =~ /^DR\s+Ensembl;\s+(ENST\d+);.*\[(\w+)-(\d+)\]/ ) {
            if ( $3 == 1 ) {
                return $1;
            }
        }
    }
    return "";
}

sub getGeneID {
    my $this = shift;
    foreach my $line ( split /\n/ , $this->entireRecord() ) {
        chomp( $line );
        if ( $line =~ m/^ID\s+(\w+)\s+\w+;\s+\d+\sAA\./ ) {
            return $1;
        }
    }
    return "";
}

sub getSites { # change this to get active sites and others
    my $this = shift;
    my $sites = {};
	print "#Found: \n";
    foreach my $feature ( split /\n/ , $this->entireRecord() ) {
        if ( $feature =~ /MOD_RES|ACT_SITE|BINDING|SITE/ ) {
            my ( $category, $position , $type ) = $feature =~ m/FT\s+(\w+)\s+(\d+)\s+\d+\s+(.*)/; 
            #if ($feature =~ /MOD_RES/) {next unless ( $type =~ /Phospho/ );} # get only phosphosites
            my $detail = "";
            $sites->{$position}->{category} = $category;
            $sites->{$position}->{type} = $type;
            if ( $type =~ m/(.*); (.*)\./ ) {
                $sites->{$position}->{type} = $1;
                $sites->{$position}->{detail} = $2;
            } elsif ( $type =~ m/(.*)\.$/ ) {
                $sites->{$position}->{type} = $1;
                $sites->{$position}->{detail} = "";
            } elsif ( $type =~ m/(.*)\.(.*)/ ) {
                $sites->{$position}->{type} = $1;
                $sites->{$position}->{detail} = $2;
            }
			print "#";
			&printPos( $sites , $position );
        } 
    }
    return $sites;
}

sub getInfo {
	my $info3d = shift;
	my $pdb = shift;
	foreach my $info ( split /\|/ , $info3d ) {
		if ( $info =~ /$pdb/ ) {
			return ( split / / , $info3d );
		}
	}
	return 0; 
}

sub printPos {
	my $sites = shift;
	my $position = shift;
	print join( "\t" , ( $position ,
						 $sites->{$position}->{type} ,
						 $sites->{$position}->{detail}
						 )
		  )."\n";
}
