#!/usr/bin/perl
#18 November 2016 - Adam Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

use TGI::Mutpro::Preprocess::Uniprot;

my $usage = 'perl phospho.ms2maf.pl <input.ms> <input.hup> <output.maf>
';

die $usage , unless @ARGV == 3;

my ( $ms , $hup , $output ) = @ARGV;

my $MS = FileHandle->new( $ms , "r" );
if ( not defined $MS ) { die "ADSERROR: Could not open/write $ms\n"; }

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
			 "RefSeq_Protein_ID" , "Specific_UniProt_ID" , 
			 "UniProt_ID" ) )."\n" );


my $msSites = {};
my $header = $MS->getline();
while ( my $line = $MS->getline() ) {
	chomp( $line );
	my ( $siteIdentifier ) = (split( /\t/ , $line ) )[0];
	my ( $hugo , $refseqProteinID , $phosphosite ) = split( /\:/ , $siteIdentifier );
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
			$msSites->{$hugo}->{$refseqProteinID}->{$phosphosite} = 1;
		}
	} else {
		$msSites->{$hugo}->{$refseqProteinID}->{$phosphosite} = 1;
	}
}
$MS->close();

foreach my $hugo ( sort keys %{$msSites} ) {
	my $uniprotID = $hugo2uniprot->{$hugo};
	print $uniprotID."\t";
	my $uniprot = new TGI::Mutpro::Preprocess::Uniprot( $uniprotID );
	my ( $protein2uniprot , $uniprot2transcript ) = &refseqProteinID2EnsemblTranscriptID( $uniprot );
	print( (scalar keys %{$protein2uniprot})."\n" );
	foreach my $refseqProteinID ( sort keys %{$msSites->{$hugo}} ) {
		my $specificUniProtID = $protein2uniprot->{$refseqProteinID};
		my $transcriptID = $uniprot2transcript->{$specificUniProtID};
		my $chromosome = &getChromosome( $uniprot );
		my @sites = sort keys %{$msSites->{$hugo}->{$refseqProteinID}};
		&printAsMAF( $OUT , \@sites , $hugo , $transcriptID , $chromosome , $refseqProteinID , $specificUniProtID , $uniprotID );
	}
}
$OUT->close();

sub refseqProteinID2EnsemblTranscriptID {
	my $this = shift;
	my $protein2uniprot = {};
	my $uniprot2transcript = {};
	my $id = $this->uniprotId();
	my $refseqLines = $this->annotations( "RefSeq" );
	my $refseqEntries = scalar @{$refseqLines};
	my $ensemblLines = $this->annotations( "Ensembl" );
	my $ensemblEntries = scalar @{$ensemblLines};
	my $overrideRefSeq = 0;
	my $overrideEnsembl = 0;
	if ( $refseqEntries > 1 && $ensemblEntries == 1 ) {
		$overrideRefSeq = 1;
	}
	if ( $refseqEntries == 1 && $ensemblEntries > 1 ) {
		$overrideEnsembl = 1;
	}
	my ( $proteinID , $uniprotID , $transcriptID ) = ( "" , $id , "" );
	foreach my $refseqLine ( @{$refseqLines} ) {
		print $refseqLine." => ";
		if ( $refseqLine =~ /$id/ ) {
			( $proteinID , $uniprotID ) = $refseqLine =~ m/(.+);\s+.+\.\s+\[(.+)\]/;
		} else {
			$refseqLine =~ m/(.+);.+\./;
			$proteinID = $1;
			$uniprotID = $id;
			if ( $overrideRefSeq ) {
				$uniprotID .= "-1";
			}
		}
		print join( "\t" , ( $proteinID , $uniprotID ) )."\n";
		$protein2uniprot->{$proteinID} = $uniprotID;
	}
	foreach my $ensemblLine ( @{$ensemblLines} ) {
		print $ensemblLine." => ";
		if ( $ensemblLine =~ /$id/ ) {
			( $transcriptID , $uniprotID ) = $ensemblLine =~ m/(.+);\s+.+;\s+.+\.\s+\[(.+)\]/;
		} else {
			$ensemblLine =~ m/(ENST\d+);.*\./;
			$transcriptID = $1;
			$uniprotID = $id;
			if ( $overrideEnsembl ) {
				$uniprotID .= "-1";
			}
		}
		print join( "\t" , ( $transcriptID , $uniprotID ) )."\n";
		$uniprot2transcript->{$uniprotID} = $transcriptID;
	}
	return ( $protein2uniprot , $uniprot2transcript );
}

sub printAsMAF {
	my $OUT = shift;
	my $sites = shift;
	my $gene = shift;
	my $trans = shift;
	my $chromosome = shift;
	my $refseqProteinID = shift;
	my $specificUniProtID = shift;
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
		$OUT->print( join( "\t" , ( $gene , $chromosome , $start , $stop, "Missense_Mutation" , "A" , "T" , "T" , $type , $trans , "p.".$aa.$site."C" , $refseqProteinID , $specificUniProtID , $uniprotID ) )."\n" );
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
