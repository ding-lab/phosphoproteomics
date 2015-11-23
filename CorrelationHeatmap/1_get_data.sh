SRC="linus300.gsc.wustl.edu:/gscuser/mwyczalk/projects/Virus/Virus_2013.9a/analysis/Music/UnifiedVirus2"
DEST="origdata"
mkdir -p $DEST

scp $SRC/?_\*-\*/dat/glm\*  $DEST
