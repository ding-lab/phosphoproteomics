### run hotspot3d ###
ln -s /gscmnt/gc2521/dinglab/Structure_Guided/Preprocessing_AvgDist/

bsub -oo search_PDBsites.log 'hotspot3d search --prep-dir Preprocessing_AvgDist --maf-file ../data_runs/PDB_site.maf --output-prefix PDB_sites --skip-silent --missense-only --3d-distance-cutoff 10'
bsub -oo search_PDB.brca.ov_sites.log 'hotspot3d search --prep-dir Preprocessing_AvgDist --maf-file ../data_runs/PDB_cptac.brca.ov_sites.maf --output-prefix PDB_cptac.brca.ov_sites --skip-silent --missense-only --3d-distance-cutoff 10'

bsubl -oo post.log 'hotspot3d post --maf-file=../data_runs/PDB_cptac.brca.ov_sites.maf --input-prefix PDB_cptac.brca.ov_sites'

bsubl -oo cluster.log 'hotspot3d cluster --collapsed-file=PDB_cptac.brca.ov_sites.pairwise.singleprotein.collapsed --pairwise-file=PDB_cptac.brca.ov_sites.pairwise --maf-file=../data_runs/PDB_site.maf'
hotspot3d summary --clusters-file=PDB_site.maf.PDB_cptac.brca.ov_sites.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters


### visualization ###
# get stats of which sites have most residues
perl /Users/khuang/Box\ Sync/PhD/cloned_repository/hotspot3d/scripts/clusterPDBPresence.pl PDB_cptac.brca.ov_sites.pairwise PDB_site.maf.PDB_cptac.brca.ov_sites.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters "PDB_cptac.brca.ov_sites"
perl /Users/khuang/Box\ Sync/PhD/cloned_repository/hotspot3d/scripts/genePDBPresence.pl PDB_cptac.brca.ov_sites.pairwise "PDB_cptac.brca.ov_sites"

mkdir pml_scripts
hotspot3d visual --clusters-file PDB_site.maf.PDB_cptac.brca.ov_sites.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters --pairwise-file PDB_cptac.brca.ov_sites.pairwise --output-file pml_scripts/MAP2K6.3ENM.hotspot.pml --pdb 3ENM --script-only
hotspot3d visual --clusters-file PDB_site.maf.PDB_cptac.brca.ov_sites.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters --pairwise-file PDB_cptac.brca.ov_sites.pairwise --output-file pml_scripts/MAP2K6.3VN9.hotspot.pml --pdb 3VN9 --script-only
