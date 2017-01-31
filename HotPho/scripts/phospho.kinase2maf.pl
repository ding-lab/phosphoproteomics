#!/usr/bin/perl
#18 November 2016 - Adam Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

use TGI::Mutpro::Preprocess::Uniprot;

my $usage = 'perl phospho.kinase2maf.pl <input.kinase> <input.hup> <output.maf>
';

die $usage , unless @ARGV == 3;

my ( $kinase , $hup , $output ) = @ARGV;

my $MS = FileHandle->new( $kinase , "r" );
if ( not defined $MS ) { die "ADSERROR: Could not open/write $kinase\n"; }

my $HUP = FileHandle->new( $hup , "r" );
if ( not defined $HUP ) { die "ADSERROR: Could not open/write $hup\n"; }

my $OUT = FileHandle->new( $output , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $output\n"; }

my $hugo2uniprot = {};
while ( my $line = $HUP->getline() ) {
	chomp( $line );
	my ( $hugo , $id ) = (split( /\t/ , $line ))[0,1];
	$hugo2uniprot->{$hugo} = $id;
}
$HUP->close();

$OUT->print( join( "\t" , ( "Hugo_Symbol" , "Chromosome" , "Start_Position" , 
			 "End_Position" , "Variant_Classification" , "Reference_Allele" , 
			 "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , 
			 "Tumor_Sample_Barcode" , "transcript_name" , "amino_acid_change" , 
			 "UniProt_ID" ) )."\n" );


my $msSites = {};
while ( my $line = $MS->getline() ) {
	next if ( $line =~ /KINASE/ );
	chomp( $line );
	my ( $id , $gene , $phosphosite ) = (split( /\t/ , $line ) )[6,7,9];
	print join( "\t" , ( $id , $gene , $phosphosite ) )."\n";
	if ( $phosphosite =~ /\D\d+\D\d+/ ) {
		my @sites;
		my $site = "";
		foreach my $char ( split // , $phosphosite ) {
			if ( $char =~ /\D/ ) {
				if ( $site ne "" ) {
					push @sites , $site;
				}
				$site = $char;
			} else {
				$site .= $char;
			}
		}
		if ( $site ne "" ) {
			push @sites , $site;
		}
		foreach my $phosphosite ( @sites ) {
			$msSites->{$id}->{$phosphosite} = 1;
		}
	} else {
		$msSites->{$id}->{$phosphosite} = 1;
	}
}
$MS->close();

foreach my $uniprotID ( sort keys %{$msSites} ) {
	my $uniprot = new TGI::Mutpro::Preprocess::Uniprot( $uniprotID );
	my $hugo = &getGeneID( $uniprot );
	my $transcriptID = &getCanonicalTranscript( $uniprot );
	my $chromosome = &getChromosome( $uniprot );
	my @sites = sort keys %{$msSites->{$uniprotID}};
	&printAsMAF( $OUT , \@sites , $hugo , $transcriptID , $chromosome , $uniprot->uniprotId() );
}
$OUT->close();

sub printAsMAF {
	my $OUT = shift;
	my $sites = shift;
	my $gene = shift;
	my $trans = shift;
	my $chromosome = shift;
	my $uniprotID = shift;
	foreach my $phosphoSite ( @{$sites} ) {
		my ( $aa , $site ) = $phosphoSite =~ m/(\w)(\d+)/;
		my $type = "Phospho";
		if ( $aa =~ m/[sS]/ ) {
			$aa = "S";
			$type .= "serine";
		} elsif ( $aa =~ /[yY]/ ) {
			$aa = "Y";
			$type .= "tyrosine";
		} elsif ( $aa =~ /[tT]/ ) {
			$aa = "T";
			$type .= "threonine";
		} elsif ( $aa =~ /[dD]/ ) {
			$aa = "D";
			$type .= "aspartate";
		} elsif ( $aa =~ /[hH]/ ) {
			$aa = "H";
			$type .= "histidine";
		}

		my $start = int( 1000000000*rand() );
		my $stop = $start;
		$OUT->print( join( "\t" , ( $gene , $chromosome , $start , $stop , 
									"Missense_Mutation" , "A" , "T" , "T" , 
									$type , $trans , "p.".$aa.$site."C" , 
									$uniprotID ) )."\n" );
	}
}

sub getChromosome {
	my $this = shift;
	my $proteomes = $this->annotations( "Proteomes" );
	my $chromosome = 1;
	foreach my $line ( @{$proteomes} ) {
		next unless ( $line =~ /Chromosome/ );
		chomp( $line );
		$chromosome = $line;
		$chromosome =~ s/.*Chromosome\s+([A-Z0-9]+)\./$1/g;
		return $chromosome;
	}
	return $chromosome;
}

sub getCanonicalTranscript {
    my $this = shift;
	my $specific = "1";
	if ( $this->uniprotId() =~ /\-/ ) {
		$this->uniprotId() =~ /\w\d+\-(\d+)/;
		$specific = $1;
	}
	print $this->uniprotId()."\t".$specific."\n";
    foreach my $line ( @{$this->annotations( "Ensembl" )} ) {
		print $line."\n";
        if ( $line =~ /(ENST\d+);.*\[(\w+)-(\d+)\]/ ) {
            if ( $3 == $specific ) {
                return $1;
            }
        }
    }
    return "ENST";
}

sub getGeneID {
    my $this = shift;
    foreach my $line ( split /\n/ , $this->entireRecord() ) {
        chomp( $line );
        if ( $line =~ m/^GN\s+Name=(\w+);/ ) {
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
