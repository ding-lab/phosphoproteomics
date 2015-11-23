# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript PvalHeatmapPlotter.R [-v] [-t title] [-h height] [-w width] [-H hm.height] [-W hm.width] [-P Pval.cutoff] 
#           [-n] [-d] [-c Pval.crop] [-o Pval.list.out] [-F filter.fn] [-l] [-R order.fn] glm.data.fn out.pdf
#
# Create a heatmap of Pvalues from MuSiC Clinical Correlation
# rows and columns are clutered according to similarity
# Pvalues indicate positive or negative correlation
#
# -t title
# -h, -w: height, width of plot in inches
# -H, -W: relative height and width of heatmap.  Between 0 and 1.
# -P Pval.cutoff: remove columns and rows where all entries have Pval > Pval.cutoff
# -d: do not plot dendrograms 
# -n: do not order heatmap (implies -d)
# -l: do not show legend
# -o out.fn: output filtered variant/trait pairs sorted by Pvalue 
# -c Pval.crop: Pvalue significance cropping. E.g., -c 1e-10 crops Pvalue of 1e-12 to 1e-10 for plotting purposes.
# -R order.fn: define name of file which stores row and column order.  If -n is not set, will write this file
#    after row/column order has been defined by similarity clustering.  If -n is set, will read from this file
#    in order to define row/column order.  This is useful when want same ordering for two different plots.
#
# Can filter out arbitrary variant/trait pairs with the -F flag.  For
# convenience in constructing this file the [-o out.fn] option can be used,
# which outputs the most significant correlations which remain after filtering
# the idea is that -o and -F can be used to iteratively construct a list of
# uninteresting variant/trait pairs It is expected that [-F filter.fn] is in
# same format as the generated output file [-o out.fn]: tab separated with
# headers, with variant in first column and trait in second.  Wildcards "*" are
# supported for variants or traits.  Third column (signed.lnP) is optional and
# unused.

# Dendrogram and correlation based on: http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
# Package "ggdendro" needs to be installed: install.packages("ggdendro", repos="http://cran.wustl.edu/")

# support for self-self correlations (where variant, trait are the same) has been removed in this version.

options("width"=180) # useful for debugging

library("reshape2")
library("plyr")
suppressPackageStartupMessages(library("ggplot2"))
library("RColorBrewer")
library("ggdendro")
library("grid")  
library("gridExtra")  


# Return the command line argument associated with a given flag (i.e., -o foo),
# or the default value if argument not specified.
# Note that this will break if an argument is not supplied after the flag.
get_val_arg = function(args, flag, default) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = args[ix+1] } else { val = default }
    return(val)
}

# Return boolean specifying whether given flag appears in command line (i.e., -o),
get_bool_arg = function(args, flag) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = TRUE } else { val = FALSE }
    return(val)
}

# Usage: 
#   args = parse_args()
#   print(args$disease.filter)
parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    title = get_val_arg(args, "-t", "")
    height = as.numeric(get_val_arg(args, "-h", "16"))
    width = as.numeric(get_val_arg(args, "-w", "16"))
    hm.height = as.numeric(get_val_arg(args, "-H", "0.9"))
    hm.width = as.numeric(get_val_arg(args, "-W", "0.9"))
    Pval.cutoff = as.numeric(get_val_arg(args, "-P", NA))
    Pval.list.fn = get_val_arg(args, "-o", NA)
    filter.fn = get_val_arg(args, "-F", NA)
    no.cluster = get_bool_arg(args, "-n")
    no.legend = get_bool_arg(args, "-l")
    no.dendrogram = get_bool_arg(args, "-d")
    Pval.crop = as.numeric(get_val_arg(args, "-c", NA))
    order.fn = get_val_arg(args, "-R", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    plot.fn = args[length(args)];               args = args[-length(args)]
    data.fn = args[length(args)];               args = args[-length(args)]

    # additional arguments logic
    if (no.cluster) no.dendrogram = TRUE

    val = list( "plot.fn" = plot.fn, "data.fn" = data.fn,
                'verbose'=verbose, 'title'=title, 
                'height'=height, 'width'=width,
                'hm.height'=hm.height, 'hm.width'=hm.width, 'Pval.cutoff'=Pval.cutoff,
                'no.cluster'=no.cluster, 'no.dendrogram'=no.dendrogram,
                'Pval.list.fn' = Pval.list.fn, 'filter.fn' = filter.fn, 'Pval.crop'=Pval.crop,
                'no.legend' = no.legend, 'order.fn'=order.fn)
    if (val$verbose) { print(val) }

    return (val)
}

# Read correlation data, optionally apply cutoff, and calculate signed lnP
get.correlation.data = function(data.fn, Pval.crop=NULL) {
    data = read.table(data.fn, row.names=NULL, header=TRUE, sep = "\t", comment.char="#")
    #                                   y y_type        x degrees_freedom  deviance residual_degrees_freedom residual_deviance   p.value coefficient covariants memo       FDR
    # 1 CLIN.Dx_merged.T.is.Adenosquamous      B   ARID1A               1 0.5670171                      190          38.31867 0.4514467   -15.78758         NA   NA 0.8308658
    data = data[,c("x", "y", "p.value", "coefficient")]
    names(data)[names(data)=="x"] = "variant"
    names(data)[names(data)=="y"] = "trait"

    data$p.value = as.numeric(data$p.value)

    # crop the Pvalues, so that their minimum (most significant) value is no less than Pval.crop
    # This has to be done row-wise
    if (!is.na(Pval.crop)) {
        data$min.value = Pval.crop
        data$p.value = apply(data[,c("p.value", "min.value")], 1, max)
        data$min.value = NULL
    }
    data$lnP = -log10(data$p.value)
    # positive correlations -> lnP > 0 -> red color
    data$signed.lnP = ifelse(data[,"coefficient"]>0, data$lnP, -data$lnP)
    return(data)
}

# read in filter.fn file where first column is variant and second column is trait. All variant/trait pairs are set to NA.
# warn user about pairs not in dataset.
# It is expected that filter.fn is in same format as the generated output file out.fn: tab separated with headers,
# with variant in first column and trait in second.  Third column (signed.lnP) is optional and unused.
apply.file.filter = function(filter.fn) {
    filter.list = read.table(filter.fn, row.names=NULL, header=TRUE, sep = "\t", comment.char="#")
    if (!(all(c("variant", "trait") %in% names(filter.list)))) {
        stop(paste(filter.fn, "does not have columns \"variant\" and \"trait\""))
    }
    # check to make sure that all variants and traits are in data
    filter.variant.match = filter.list$variant %in% c(as.character(data$variant), "*")
    filter.trait.match = filter.list$trait %in% c(as.character(data$trait), "*")
    if (any(!filter.variant.match)) {
        unknown.variants = paste(filter.list$variant[!filter.variant.match])
        stop(paste("unknown variants in ", filter.fn, ":", unknown.variants))
    }
    if (any(!filter.trait.match)) {
        unknown.traits = paste(filter.list$trait[!filter.trait.match])
        stop(paste("unknown traits in ", filter.fn, ":", unknown.traits))
    }

    for (i in 1:nrow(filter.list)) {
        v=as.character(filter.list[i,"variant"])
        t=as.character(filter.list[i,"trait"])
#        cat(paste("Removing ", v, ",", t,"\n"))
        if ( (t == "*") & (v == "*")) stop("variant and trait both *.  Quitting.")
        if ( t == "*") {
            data[data$variant == v ,"signed.lnP"] = NA
        } else if (v == "*") {
            data[data$trait == t,"signed.lnP"] = NA
        } else {
            data[data$variant == v & data$trait == t,"signed.lnP"] = NA
        }
    }
    return(data)
}

# remove data based on values:
# * remove rows where all values are NA
# * remove columns, rows where all values are below a given Pval.cutoff
apply.value.filter = function(data.wide, Pval.cutoff=NA) {
    # remove genes which have all NA (e.g., "OTHER")
    # from http://stackoverflow.com/questions/4862178/remove-rows-with-nas-in-data-frame
    # 1 indicates rows, 2 indicates columns
    row.all.na = apply(data.wide, 1, function(x){all(is.na(x))})
    if (any(row.all.na)) {
        cat("NOTE: Removing rows because all NA: \n")
        print(row.names(data.wide)[row.all.na])
        data.wide = data.wide[!row.all.na,]
    }

    # Remove columns, rows with all |lnP| < ln(Pval.cutoff)
    if (!is.na(Pval.cutoff)) {
        lnP.cutoff = -log10(Pval.cutoff)
        row.all.below.cutoff = apply(data.wide, 1, function(x){all( na.omit(abs(x)) < lnP.cutoff )})
        col.all.below.cutoff = apply(data.wide, 2, function(x){all( na.omit(abs(x)) < lnP.cutoff )})
        n.row.removed = sum(row.all.below.cutoff, na.rm=TRUE)
        n.col.removed = sum(col.all.below.cutoff, na.rm=TRUE)
        cat( paste("Removing", n.row.removed, "rows and", n.col.removed, "columns with all Pval <", Pval.cutoff,"\n") )
        data.wide = data.wide[!row.all.below.cutoff,]
        data.wide = data.wide[,!col.all.below.cutoff]
    }
    return(data.wide)
}

write.filter.file = function(data.wide, Pval.list.fn) {
    data.wide$variant = row.names(data.wide)
    data.filtered = melt(data.wide, id.vars="variant", variable.name="trait", value.name="signed.lnP")
    # sort by magnitude of lnP, and ignore NA entries.  
    data.filtered = na.omit(data.filtered[order(-abs(data.filtered$signed.lnP)), ])
    write.table(data.filtered, Pval.list.fn, sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("    Saved ordered list of P-values to %s\n", Pval.list.fn))
}

write.order.file = function(ordered.data.melt, order.fn) {
    # construct data frame consisting of levels of variants and traits in ordered.data.melt data frame
    variant.df = data.frame(name=levels(ordered.data.melt$variant), type="variant", order=seq_along(levels(ordered.data.melt$variant)))
    trait.df = data.frame(name=levels(ordered.data.melt$trait), type="trait", order=seq_along(levels(ordered.data.melt$trait)))
    print(order.fn)
    write.table(rbind(variant.df, trait.df), order.fn, sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("    Saved list of ordered variants and traits to %s\n", order.fn))
}

# evaluate row (if is.row=TRUE) or column similarity and find relative positions (data.order)
# construct data used to construct dendrograms
# return list (data.order, dendrogram.data)
get.dendrogram.order = function(dm, is.row) {
# core code from: http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
    if (is.row)
        dd = as.dendrogram(hclust(dist(dm)))
    else
        dd = as.dendrogram(hclust(dist(t(dm))))
    data.order = order.dendrogram(dd)
    dendrogram.data = dendro_data(dd)
    return(list('data.order'=data.order, 'dendrogram.data'=dendrogram.data))
}

get.file.order = function(order.fn, dm.cols, dm.rows) {
    order.dat = read.table(order.fn, row.names=NULL, header=TRUE, sep = "\t", comment.char="#", stringsAsFactors=FALSE)
#       name    type order
# 1  Any.HHV variant     1

    order.cols = order.dat[order.dat$type == 'trait','name']
    # check to make sure that requested column names match actual column names
    if (length(order.cols) != length(dm.cols))
        stop(paste("Incorrect traits count in", order.fn))
    missing.cols = !(order.cols %in% dm.cols)
    if (any(missing.cols))
        stop(paste("Unknown traits in", order.fn, order.cols[missing.cols]))

    # get new order for dm.cols by merging with requested order and reordering.  There may be a better way to do this.
    merged.cols = merge(data.frame(name=dm.cols), order.dat[order.dat$type == 'trait',])
    merged.cols = merged.cols[order(merged.cols$order),]
    col.ord = as.numeric(row.names(merged.cols))

    # now do the same thing with rows
    order.rows = order.dat[order.dat$type == 'variant','name']
    if (length(order.rows) != length(dm.rows))
        stop(paste("Incorrect traits count in", order.fn))
    missing.rows = !(order.rows %in% dm.rows)
    if (any(missing.rows))
        stop(paste("Unknown traits in", order.fn, order.rows[missing.rows]))
    merged.rows = merge(data.frame(name=dm.rows), order.dat[order.dat$type == 'variant',])
    merged.rows = merged.rows[order(merged.rows$order),]
    row.ord = as.numeric(row.names(merged.rows))

    return(list('row.ord'=row.ord, 'col.ord'=col.ord))
}

# define column and row order in one of three ways:
# 1. by similiarity.  Also creates dendrogram information.
# 2. from a file which defines the order.
# 3. default order.
#
# If do.dendrogram.order is defined, the first method is chosen; if order.fn is not null, the second
# is chosen; otherwise, the third.
get.ordered.data = function(data.wide, do.dendrogram.order, order.fn=NULL) {
    dm = as.matrix(data.wide)
    dendrogram.data.cols = NULL
    dendrogram.data.rows = NULL

#print(dm[,1:10])
#print(colnames(dm))
#print(rownames(dm))
    if (do.dendrogram.order) {
        cols = get.dendrogram.order(dm, FALSE) 
        col.ord = cols$data.order 
        dendrogram.data.cols = cols$dendrogram.data
        rows = get.dendrogram.order(dm, TRUE) 
        row.ord = rows$data.order 
        dendrogram.data.rows = rows$dendrogram.data
    } else if (!is.null(order.fn)) {
        order.data = get.file.order(order.fn, colnames(dm), rownames(dm))
        col.ord = order.data$col.ord
        row.ord = order.data$row.ord
    } else {
        # retain existing row and column order if not clustering
        col.ord = seq_along(colnames(dm))
        row.ord = seq_along(rownames(dm))
    }
    # re-order heatmap names according to dendrogram (or trivial order if not reordering)
    dm = dm[row.ord, col.ord] 
    df = as.data.frame(dm)  
    df$variant = rownames(df)
    df$variant = with(df, factor(variant, levels=variant, ordered=TRUE))
    ordered.data.melt = melt(df, id.vars="variant", variable.name="trait", value.name="signed.lnP")
    return(list('ordered.data.melt'=ordered.data.melt, 'dendrogram.data.cols'=dendrogram.data.cols, 'dendrogram.data.rows'=dendrogram.data.rows))
}

get.dendrogram.plot = function(dendrogram.data.cols, dendrogram.data.rows) {
    theme_none = theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_text(colour=NA),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(), 
      axis.ticks = element_blank())

    dendrogram.top = ggplot(segment(dendrogram.data.cols)) + 
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="gray50") + 
      theme_none + theme(axis.title.x=element_blank(), 
      plot.margin=unit(c(0,0,-0.15,0), "npc"))  # top, right, bottom, and left

    dendrogram.right = ggplot(segment(dendrogram.data.rows)) + 
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="gray50") + 
      coord_flip() + theme_none + theme(plot.margin=unit(c(0,0,0,-0.15), "npc"))  # top, right, bottom, and left
    return(list('dendrogram.top'=dendrogram.top, 'dendrogram.right'=dendrogram.right))
}

get.gradient.fill = function(ordered.data.melt, quantitative=TRUE) {
    # either of these works:
    if (quantitative) 
        scale = scale_fill_gradient2( low = "#377EB8", mid = "white", high = "#E41A1C", name="log10(P)")
    else {
        # this from /Users/mwyczalk/Documents/Manuscripts/PanGermline.Resubmission.May2015/Figure4/Figure4E/src/CorrelationHeatmapPlotter.R
        # want to set up heatmap so it ranges -A:A (i.e., 0 is midpoint)
        max.lnP = max(abs(ordered.data.melt$signed.lnP))
        RdBu.11 = brewer.pal(11, "RdBu")
        RdBu.11[6]="#FFFFFF"
        getPalette = colorRampPalette(RdBu.11)
        scale = scale_fill_gradientn(colours=rev(getPalette(100)), limits=c(-max.lnP,max.lnP), breaks=c(-0.8*max.lnP, 0.8*max.lnP), labels=c("Negative", "Positive"), name="Correlation")
    }
    return(scale)
}

get.heatmap.plot = function(ordered.data.melt, quantitative.scale=TRUE) {
    p.heatmap = ggplot(ordered.data.melt, aes(x=trait, y=variant)) + geom_tile(aes(fill=signed.lnP)) 
    p.heatmap = p.heatmap + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=6, hjust=1), axis.text.y = element_text(size=5))
    p.heatmap = p.heatmap + theme(axis.title = element_blank())
    p.heatmap = p.heatmap + theme(plot.margin=unit(c(0,0,0,0), "npc"), panel.border = element_blank(), panel.background=element_blank(), axis.ticks = element_blank())
    p.heatmap = p.heatmap + get.gradient.fill(ordered.data.melt, quantitative=TRUE) 
    return(p.heatmap)
}

arrange.heatmap.dendrogram = function(p.heatmap, dendrogram.top, dendrogram.right, title, hm.width, hm.height) {
    # Align the dendrograms with the heatmap.  See /Users/mwyczalk/Data/Virus/Virus_2013.9a/CombinedPlot/BreakpointSurveyor/S_Draw/src/BreakpointDrawer.R
    # for details
    p.blank = rectGrob(gp = gpar(col = "white"))

    heatmap.gb = ggplot_build(p.heatmap)
    dend.top.gb = ggplot_build(dendrogram.top)
    dend.right.gb = ggplot_build(dendrogram.right)

    dend.top.gb$panel$ranges[[1]]$x.range = heatmap.gb$panel$ranges[[1]]$x.range
    dend.right.gb$panel$ranges[[1]]$y.range = heatmap.gb$panel$ranges[[1]]$y.range

    heatmap.gt = ggplot_gtable(heatmap.gb)
    dend.top.gt = ggplot_gtable(dend.top.gb)
    dend.right.gt = ggplot_gtable(dend.right.gb)

    dend.top.gt$widths = heatmap.gt$widths
    dend.right.gt$heights = heatmap.gt$heights

    main.grob = arrangeGrob(dend.top.gt, p.blank, heatmap.gt, dend.right.gt, 
        widths=c(hm.width,1-hm.width), heights=c(1-hm.height,hm.height), ncol=2, nrow=2)

    if (!is.null(title)) {
        title.ggp = textGrob(title)
        grid.arrange(title.ggp, main.grob, heights = c(0.025, 0.975), ncol=1, nrow=2, newpage=FALSE)
    } else {
        print(main.grob)
    }

}

#options("width"=270) 
args = parse_args()

data = get.correlation.data(args$data.fn, args$Pval.crop)

if (!is.na(args$filter.fn)) {
    data = apply.file.filter(args$filter.fn)
}

data.wide = dcast(data, variant ~ trait, value.var = "signed.lnP")
rownames(data.wide) = data.wide$variant
data.wide$variant = NULL

data.wide = apply.value.filter(data.wide, args$Pval.cutoff)

# output a list of variant/trait pairs along with their Pvalues if the -o flag is set
# this is used for creating the filter file read previously in apply.file.filter()
if (!is.na(args$Pval.list.fn)) {
    write.filter.file(data.wide, args$Pval.list.fn)
}

# no.cluster defined implies no ordering of data by similarity
# args$no.dendrogram implies dendrogram not displayed
order.data = get.ordered.data(data.wide, !args$no.cluster, args$order.fn)
ordered.data.melt = order.data$ordered.data.melt

# Write an order.fn file is 1) order.fn is defined and 2) no.cluster is not defined
# this is used in cases where one plot needs to define the ordering for another
if (!args$no.cluster & !is.null(args$order.fn)) {
    write.order.file(ordered.data.melt, args$order.fn)
}


### Plotting begins ###
p.heatmap = get.heatmap.plot(ordered.data.melt, quantitative.scale=TRUE)

if (!args$no.dendrogram) 
    dendrograms = get.dendrogram.plot( order.data$dendrogram.data.cols, order.data$dendrogram.data.rows)

if (!is.null(args$title) & args$no.dendrogram)
    p.heatmap = p.heatmap + ggtitle(args$title)

if (args$no.legend) 
    p.heatmap = p.heatmap + theme(legend.position="none")

pdf(file=args$plot.fn, width=args$width, height=args$height, useDingbats=FALSE)

if (args$no.dendrogram) {
    print(p.heatmap)
} else {
    arrange.heatmap.dendrogram(p.heatmap, dendrograms$dendrogram.top, dendrograms$dendrogram.right, args$title, args$hm.width, args$hm.height)
}

cat(sprintf("Heatmap Saved to %s\n", args$plot.fn))
