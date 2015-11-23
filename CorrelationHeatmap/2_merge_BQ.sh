# Combine clinical correlation data

# B file:
#     1  y
#     2  y_type
#     3  x
#     4  degrees_freedom
#     5  deviance
#     6  residual_degrees_freedom
#     7  residual_deviance
#     8  p-value
#     9  coefficient
#    10  covariants
#    11  memo
#    12  FDR
#
# Q file:
#     1  y
#     2  y_type
#     3  x
#     4  degrees_freedom
#     5  deviance
#     6  residual_degrees_freedom
#     7  residual_deviance
#*    8  F_statistic
#     9  p-value
#    10  coefficient
#    11  covariants
#    12  memo
#    13  FDR

# combine B and Q files, stripping header off second one
# Note that headers of B, Q differ
# Usage BFN QFN fout
function combine_BQ {
BFN=$1
QFN=$2
FOUT=$3

cp $BFN $FOUT
cut -f 1-7,9- $QFN | grep -v "p-value" >> $FOUT
echo Written to $FOUT
}

# normalize Q files when B not present
function process_Q {
QFN=$1
FOUT=$2

cut -f 1-7,9- $QFN  > $FOUT
echo Written to $FOUT
}

mkdir -p dat

# combine B and Q data
combine_BQ origdata/glm.Integration-Clinical.B.dat origdata/glm.Integration-Clinical.Q.dat dat/glm.Integration-Clinical.dat
combine_BQ origdata/glm.MAF-Clinical.B.dat origdata/glm.MAF-Clinical.Q.dat dat/glm.MAF-Clinical.dat
combine_BQ origdata/glm.Virus-Clinical.B.dat origdata/glm.Virus-Clinical.Q.dat dat/glm.Virus-Clinical.dat

# process Q and copy to ./dat
process_Q origdata/glm.Integration-RPPA.Q.dat dat/glm.Integration-RPPA.dat
process_Q origdata/glm.Integration-RSEM.Q.dat dat/glm.Integration-RSEM.dat
process_Q origdata/glm.MAF-RPPA.Q.dat dat/glm.MAF-RPPA.dat
process_Q origdata/glm.MAF-RSEM.Q.dat dat/glm.MAF-RSEM.dat
process_Q origdata/glm.Virus-RPPA.Q.dat dat/glm.Virus-RPPA.dat
process_Q origdata/glm.Virus-RSEM.Q.dat dat/glm.Virus-RSEM.dat
process_Q origdata/glm.Methylation-RSEM.Q.dat dat/glm.Methylation-RSEM.dat
process_Q origdata/glm.Virus-Methylation.Q.dat dat/glm.Virus-Methylation.dat

# simply copy B to dat, normalizing name
cp origdata/glm.MAF-Integration.B.dat dat/glm.MAF-Integration.dat
cp origdata/glm.MAF-Virus.B.dat dat/glm.MAF-Virus.dat

