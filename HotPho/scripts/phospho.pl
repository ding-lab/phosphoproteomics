#!/usr/bin/perl
#18 November 2016 - Adam Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

use TGI::Mutpro::Preprocess::Uniprot;

my $usage = 'perl phospho2maf.pl #uniprot-id# <output.maf>
';

die $usage , unless @ARGV == 2;

my $id = shift;
my $output = shift;
my $uniprot = new TGI::Mutpro::Preprocess::Uniprot( $id );
my $sites = &getPhosphosites( $uniprot );
my ( $gene , $species ) = split( /\_/ , &getGeneID( $uniprot ) );
my $trans = &getCanonicalTranscript( $uniprot );
my $chromosome = &getChromosome( $uniprot );

my $OUT = FileHandle->new( $output , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $output\n"; }

foreach my $site ( sort {$a<=>$b} keys %{$sites} ) {
	my $type = $sites->{$site}->{type};
	my $aa = "X";
	if ( $type =~ /serine/ ) {
		$aa = "S";
	} elsif ( $type =~ /tyrosine/ ) {
		$aa = "Y";
	} elsif ( $type =~ /threonine/ ) {
		$aa = "T";
	} elsif ( $type =~ /aspar/ ) {
		$aa = "D";
	} elsif ( $type =~ /histidine/ ) {
		$aa = "H";
	}

	$OUT->print( join( "\t" , ( $gene , $chromosome , 1 , 1 , "Missense_Mutation" , "A" , "T" , "T" , $type , $trans , "p.".$aa.$site."C" ) )."\n" );
}
$OUT->close();

sub getChromosome {
	my $this = shift;
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

sub getPhosphosites {
    my $this = shift;
    return &getModifiedResidues( $this , "Phospho" );
}

sub getModifiedResidues {
    my $this = shift;
    my $description = shift;
    my $sites = {};
	print "#Found: \n";
    foreach my $feature ( split /\n/ , $this->entireRecord() ) {
        if ( $feature =~ /MOD_RES/ ) {
            my ( $position , $type ) = $feature =~ m/FT\s+MOD_RES\s+(\d+)\s+\d+\s+(.*)/;
            next unless ( $type =~ /$description/ );
            my $detail = "";
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
