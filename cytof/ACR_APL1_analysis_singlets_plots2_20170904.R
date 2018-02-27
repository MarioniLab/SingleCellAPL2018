## Rscript to make plots of differential abundance analysis on APL1 data

library(ncdfFlow)
library(cydar)
library(flowCore)

############################## functions #################################

plotCellGeneric <- function (x, y, generic_value, grange = NULL, length.out = 100, pch = 16, 
                             virid_option='D',
                             ...) 
{
  require(viridis)
  if (is.null(grange)) {
    grange <- range(generic_value)
  }
  else {
    generic_value[generic_value > grange[2]] <- grange[2]
    generic_value[generic_value < grange[1]] <- grange[1]
  }
  all.cols <- viridis(n = length.out, option=virid_option)
  mdpts <- seq(from = grange[1], to = grange[2], length.out = length.out)
  actual.threshold <- mdpts + (mdpts[2] - mdpts[1])/2
  ix <- pmin(length.out, findInterval(generic_value, actual.threshold) + 
               1)
  cur.cols <- all.cols[ix]
  if (pch %in% 21:25) {
    plot(x, y, pch = pch, bg = cur.cols, ...)
  }
  else {
    plot(x, y, pch = pch, col = cur.cols, ...)
  }
  names(all.cols) <- mdpts
  return(invisible(all.cols))
}

############################## plotting script #################################


setwd('~/Documents/Post_doc/Sync/Experiments/CyTOF/APL 1/20161209-12 APL1 CyTOF/analysis2/')
load('ACR_APL1_singletfilt2_inprogress_20170904.RData')

set.seed(100)
library(Rtsne)
t.out <- Rtsne(intensities(cd2), perplexity=50)

## save the tsne coordinates for later use if want to make more plots
saveRDS(t.out, file='ACR_APL1_tSNE_singletfilt2_allspheres_20170904.rds')

library(edgeR)
pdf('~/Documents/Post_doc/Sync/Experiments/CyTOF/APL 1/20161209-12 APL1 CyTOF/analysis2/APL_v_NP68_singletfilt2_20170904.pdf', height=9, width=3)
layout(cbind(c(1:4),5), widths=c(8, 1.5))
par(mar=c(2.1, 2.1, 2.1, 1.1), oma=c(3,3,1,1), mgp=c(2.2,1,0))
adjc <- cpm(y, log=TRUE, prior=3)
for (peptide in c('N4', 'T4', 'Q4H7', 'G4')) {
  logfc <- res.all$table[is.sig.all,paste0("logFC.", peptide)]
  col <- plotCellLogFC(t.out$Y[is.sig.all,1], t.out$Y[is.sig.all,2], logfc, max.logFC=5,
                       cex.axis=1.5, cex.lab=1.5, 
                       cex.main=1.5, cex=1, main=peptide)
}
par(mar=c(0,0,0,0))
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(col))
rect(-0.2, start.loc, 0.5, start.loc+diff(start.loc)[1], col=col, border=col)
text(0.1, -0.5, pos=1, as.character(round(as.numeric(names(col)[1]), 2)), cex=1.2)
text(0.1, 0.5,  pos=3, as.character(round(as.numeric(names(col)[length(col)]), 2)), cex=1.2)
text(-0.7, 0, pos=1, srt=90, "Log-FC", cex=1.5)
mtext('t-SNE1', side=1, line=1, outer=TRUE, cex=1.2)
mtext('t-SNE2', side=2, line=1, out=TRUE, cex=1.2)
dev.off()

## Now look at the intensities so we can see how these hyperspheres are defined in terms of markers.


pdf('~/Documents/Post_doc/Sync/Experiments/CyTOF/APL 1/20161209-12 APL1 CyTOF/analysis2/Intensity_plots_sig_singletfilt2_20170904.pdf', height=20, width=18)
rownames(channels) <- channels[,1]
reranges <- intensityRanges(cd2)
lmat <- cbind(matrix(seq_len(5*4), ncol=4, nrow=5), 21)
par(oma=c(4.1,4.1,1.1,1.1))
layout(lmat, widths=c(rep(1, 4), 0.2))
for (i in order(colnames(intensities(cd2)))) {
  par(mar=c(2.1, 2.1, 2.1, 2.1))
  col <- plotCellIntensity(t.out$Y[is.sig.all,1], t.out$Y[is.sig.all,2], intensities(cd2)[is.sig.all,i],
                           irange=reranges[,i], main=as.character(channels[colnames(intensities(cd2)),2])[i],
                           cex=1.5, cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
}
for (j in seq_len(20-ncol(intensities(cd2)))) { plot.new() }
par(mar=c(0,0,0,0))
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(col))
interval <- diff(start.loc)[1]
rect(-0.5, start.loc, 0.5, start.loc+interval, col=col, border=col)
text(0, -0.5, pos=1, "Low", cex=1.5)
text(0, 0.5+interval,  pos=3, "High", cex=1.5)
text(-0.9, 0, pos=1, srt=90, "Marker intensity", cex=1.5)
mtext(text='t-SNE1', side=1, line=2, outer=TRUE, cex=1.5)
mtext(text='t-SNE2', side=2, line=1, outer=TRUE, cex=1.5)
dev.off()

## plot intensities of specific genes
reranges2 <- reranges[,c('Dy163Di', 'Gd155Di', 'Yb173Di', 'Sm154Di', 'Dy161Di', 'Eu151Di')]
int_filt <- intensities(cd2)[,c('Dy163Di', 'Gd155Di', 'Yb173Di', 'Sm154Di', 'Dy161Di', 'Eu151Di')]
pdf('~/Documents/Post_doc/Sync/Experiments/CyTOF/APL 1/20161209-12 APL1 CyTOF/analysis2/Intensity_plots_sig_filtered_singletfilt2_20170904.pdf', height=7.2, width=5.5)
lmat <- cbind(c(1,3,5), c(2,4,6), 7)
layout(lmat, widths=c(rep(8, 2), 2))
par(mar=c(2.1, 2.1, 2.1, 1.1), oma=c(3,3,1,1), mgp=c(2.2,1,0))
for (i in 1:ncol(int_filt)) {
  col <- plotCellIntensity(t.out$Y[is.sig.all,1], t.out$Y[is.sig.all,2], int_filt[is.sig.all,i],
                           irange=reranges2[,i], 
                           main=strsplit(as.character(channels[colnames(int_filt),2])[i], '_')[[1]][2], 
                           cex=1, cex.axis=1.5, 
                           cex.main=1.5, xlab="", ylab="")
}
par(mar=c(0,0,0,0))
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(col))
interval <- diff(start.loc)[1]
rect(-0.3, start.loc, 0.3, start.loc+interval, col=col, border=col)
text(0, -0.5, pos=1, "Low", cex=1.2)
text(0, 0.5+interval,  pos=3, "High", cex=1.2)
text(-0.8, 0, pos=1, srt=90, "Marker intensity", cex=1.5)
mtext('t-SNE1', side=1, line=1, outer=TRUE, cex=1.2)
mtext('t-SNE2', side=2, line=1, out=TRUE, cex=1.2)
dev.off()

## and now plot the relative proportion of each hypersphere occupied by cells of each sample 
## (after normalizing for cell number in each sample)

props <- t(apply(assay(cd2), 1, function(x) {x/cd2$totals * 100}))
peptide <- samples[colnames(cd2), 'peptide']
props_avg <- data.frame(N4=apply(props[,peptide %in% 'N4'], 1, mean),
                        T4=apply(props[,peptide %in% 'T4'], 1, mean),
                        Q4H7=apply(props[,peptide %in% 'Q4H7'], 1, mean),
                        G4=apply(props[,peptide %in% 'G4'], 1, mean),
                        NP68=apply(props[,peptide %in% 'NP68'], 1, mean))

rel_props <- data.frame(t(apply(props_avg, 1, function(x) {x/sum(x)})))


pdf('~/Documents/Post_doc/Sync/Experiments/CyTOF/APL 1/20161209-12 APL1 CyTOF/analysis2/Relative_proportions_HSs_sig_singletfilt2_horizontal_20170904.pdf', 
    height=2.8, width=12)
layout(matrix(c(1:6), nrow=1, ncol=6), widths=c(8, 8, 8, 8, 8, 1.5))
par(mar=c(2.1, 2.1, 2.1, 1.1), oma=c(3,3,1,1), mgp=c(2.2,1,0))
max_prop <- max(rel_props)
for(i in 1:(ncol(rel_props))){
  col <- plotCellGeneric(t.out$Y[is.sig.all,1], t.out$Y[is.sig.all,2], generic_value = rel_props[is.sig.all,i], 
                         grange=c(0,max_prop), virid_option = 'B', 
                         cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex=1, 
                         main=colnames(rel_props)[i])
}
par(mar=c(0,0,0,0))
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(col))
rect(-0.2, start.loc, 0.5, start.loc+diff(start.loc)[1], col=col, border=col)
text(0, -0.5, pos=1, as.character(round(as.numeric(names(col)[1]), 2)), cex=1.2)
text(0, 0.5,  pos=3, as.character(round(as.numeric(names(col)[length(col)]), 2)), cex=1.2)
text(-0.7, 0, pos=1, srt=90, "Relative Fraction", cex=1.5)
mtext('t-SNE1', side=1, line=1, outer=TRUE, cex=1.2)
mtext('t-SNE2', side=2, line=1, out=TRUE, cex=1.2)
dev.off()

