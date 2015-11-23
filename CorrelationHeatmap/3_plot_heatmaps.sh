mkdir -p plots


# Usage: Rscript PvalHeatmapPlotter.R [-v] [-t title] [-h height] [-w width] [-H hm.height] [-W hm.width] [-P Pval.cutoff] [-f] [-n] [-d] [-o Pval.list.out] [-F filter.fn] 
#    PvalHeatmapPlotter.R glm.data.fn out.pdf
#
# Create a heatmap of Pvalues from MuSiC Clinical Correlation
# rows and columns are clutered according to similarity
# Pvalues indicate positive or negative correlation
# -t title
# -h, -w: height, width of plot in inches
# -H, -W: relative height and width of heatmap.  Between 0 and 1.
# -P Pval.cutoff: remove columns and rows where all entries have Pval > Pval.cutoff
# -f: perform filter to get rid of trivial correlations
# -d: do not plot dendrograms 
# -n: do not order heatmap (implies -d)
# -l: do not show legend
# -o out.fn: output filtered variant/trait pairs sorted by Pvalue 
# -c Pval.crop: Pvalue significance cropping. E.g., -c 1e-10 crops Pvalue of 1e-12 to 1e-10 for plotting purposes.

PVAL="-P 0.01"

ARGS="$PVAL -h 8 -w 18"

BIN="Rscript src/PvalHeatmapPlotter.R"

$BIN -c 1e-30 -l $PVAL -t "Virus vs. Methylation" -w 18 -h 10 -H 0.8 dat/glm.Virus-Methylation.dat plots/Virus-Methylation.pdf
$BIN -c 1e-30 -n $PVAL -t "Virus vs. Methylation" -w 18 -h 10 -H 0.8 dat/glm.Virus-Methylation.dat plots/Virus-Methylation-noclust.pdf

$BIN -n $PVAL -c 1e-30 -t "Methylation vs. RSEM" -w 15 -h 13 dat/glm.Methylation-RSEM.dat plots/Methylation-RSEM-noclust.pdf
$BIN -l $PVAL -c 1e-30 -t "Methylation vs. RSEM" -w 18 -h 18 dat/glm.Methylation-RSEM.dat plots/Methylation-RSEM.pdf

$BIN -n $PVAL -t "Integration vs. Clinical" -w 8 -h 6 dat/glm.Integration-Clinical.dat plots/Integration-Clinical-noclust.pdf
$BIN -l $PVAL -t "Integration vs. Clinical" -w 8 -h 6 dat/glm.Integration-Clinical.dat plots/Integration-Clinical.pdf

$BIN -n $PVAL -t "MAF vs. Clinical" -w 12 -h 12 dat/glm.MAF-Clinical.dat plots/MAF-Clinical-noclust.pdf
$BIN -l $PVAL -t "MAF vs. Clinical" -w 12 -h 12 dat/glm.MAF-Clinical.dat plots/MAF-Clinical.pdf

$BIN -n $PVAL -t "Virus vs. Clinical" -w 8 -h 6 -H 0.8 dat/glm.Virus-Clinical.dat plots/Virus-Clinical-noclust.pdf
$BIN -l $PVAL -t "Virus vs. Clinical" -w 8 -h 6 -H 0.8 dat/glm.Virus-Clinical.dat plots/Virus-Clinical.pdf

$BIN -n $PVAL -t "Integration vs. RPPA" -w 8 -h 6 -H 0.8 dat/glm.Integration-RPPA.dat plots/Integration-RPPA-noclust.pdf
$BIN -l $PVAL -t "Integration vs. RPPA" -w 8 -h 6 -H 0.8 dat/glm.Integration-RPPA.dat plots/Integration-RPPA.pdf

$BIN -n $PVAL -t "Integration vs. RSEM" -w 8 -h 6 -H 0.8 dat/glm.Integration-RSEM.dat plots/Integration-RSEM-noclust.pdf
$BIN -l $PVAL -t "Integration vs. RSEM" -w 8 -h 6 -H 0.8 dat/glm.Integration-RSEM.dat plots/Integration-RSEM.pdf

$BIN -n $PVAL -t "MAF vs. RPPA" -w 18 -h 10 -H 0.8 dat/glm.MAF-RPPA.dat plots/MAF-RPPA-noclust.pdf
$BIN -l $PVAL -t "MAF vs. RPPA" -w 18 -h 10 -H 0.8 dat/glm.MAF-RPPA.dat plots/MAF-RPPA.pdf

$BIN -n $PVAL -t "MAF vs. RSEM" -w 18 -h 12 dat/glm.MAF-RSEM.dat plots/MAF-RSEM-noclust.pdf
$BIN -l $PVAL -t "MAF vs. RSEM" -w 18 -h 12 dat/glm.MAF-RSEM.dat plots/MAF-RSEM.pdf

$BIN -n $PVAL -t "Virus vs. RPPA" -w 18 -h 10 -H 0.8 dat/glm.Virus-RPPA.dat plots/Virus-RPPA-noclust.pdf
$BIN -l $PVAL -t "Virus vs. RPPA" -w 18 -h 10 -H 0.8 dat/glm.Virus-RPPA.dat plots/Virus-RPPA.pdf

$BIN -n $PVAL -t "Virus vs. RSEM" -w 18 -h 8 -H 0.8 dat/glm.Virus-RSEM.dat plots/Virus-RSEM-noclust.pdf
$BIN -l $PVAL -t "Virus vs. RSEM" -w 18 -h 8 -H 0.8 dat/glm.Virus-RSEM.dat plots/Virus-RSEM.pdf

$BIN -n $PVAL -t "MAF vs. Virus" -w 8 -h 6 -H 0.8 dat/glm.MAF-Virus.dat plots/MAF-Virus-noclust.pdf
$BIN -l $PVAL -t "MAF vs. Virus" -w 8 -h 6 -H 0.8 dat/glm.MAF-Virus.dat plots/MAF-Virus.pdf

$BIN -n -P 0.1 -t "MAF vs. Integration" -w 8 -h 6 -H 0.8 dat/glm.MAF-Integration.dat plots/MAF-Integration-noclust.pdf
$BIN -l -P 0.1 -t "MAF vs. Integration" -w 8 -h 6 -H 0.8 dat/glm.MAF-Integration.dat plots/MAF-Integration.pdf

