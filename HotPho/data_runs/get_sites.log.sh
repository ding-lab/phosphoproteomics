bsubl -oo run_phospho.hup2maf.pl.log 'perl ../scripts/phospho.hup2maf.pl /gscmnt/gc2706/dinglab/medseq/HotSpot3D/AvgDistance/Prepr*/hugo.uniprot.pdb.transcript.csv PDB_phosphosite.maf'
bsubl -oo run_sites.hup2maf.pl.log 'perl ../scripts/sites.hup2maf.pl /gscmnt/gc2706/dinglab/medseq/HotSpot3D/AvgDistance/Prepr*/hugo.uniprot.pdb.transcript.csv PDB_site.maf'

cat PDB_site.maf /gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/pan3can_shared_data/maf_w_enst/cptac.brca.phosphosites.mafish /gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/pan3can_shared_data/maf_w_enst/cptac.ov.phosphosites.mafish > PDB_cptac.brca.ov_sites.maf
