## path to R???

anno.files <- c("/lustre/jmlab/resources/annotation/original/ERCC92.gtf", "/lustre/jmlab/resources/annotation/processed/mm10.gtf")
bam.files <- list.files("bam", full=TRUE, pattern="*r_1.bam$")
stat.file <- "all_qual.tsv"
ispet <- FALSE

# Files of interest.
bam.files
anno.files
stat.file

if (!exists("ispet")) { 
    ispet <- FALSE
} 
if (!exists("strandspec")) {
    strandspec <- 0
}

ispet
strandspec

if (length(anno.files)==1L) {
    file.symlink(anno.files, "temp.gtf")
} else {
    system(paste(c("cat", anno.files, "> temp.gtf"), collapse=" "))
}

# Running featureCounts.

## For some reason, sample SLX-12612.i724_i513.HHH7TBBXX.s_7.r_1.bam is failing every time we try to count it. But for some reason when I run up to it and then again past it, it all seems fine. I need to switch to the new cluster... perhaps a memory issue?

bad_index <- which(bam.files %in% 'bam/SLX-12612.i724_i513.HHH7TBBXX.s_7.r_1.bam')

require(Rsubread)
outa <- featureCounts(bam.files[1:bad_index-1], annot.ext="temp.gtf", isGTFAnnotationFile=TRUE, minMQS=10, nthreads=4, isPairedEnd=ispet, strandSpecific=strandspec, ignoreDup=FALSE)

# Saving counts to file, with gene names.
colnames(outa$counts) <- sub("\\.bam$", "", basename(bam.files[1:bad_index-1]))
final <- data.frame(GeneID=rownames(outa$counts), Length=outa$annotation$Length, outa$counts)
write.table(file="genic_counts_leaveDupa.tsv", final, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
save(outa, file="counts_leaveDupa.RData")

rm(outa, final)

outb <- featureCounts(bam.files[(bad_index):length(bam.files)], annot.ext="temp.gtf", isGTFAnnotationFile=TRUE, minMQS=10, nthreads=4, isPairedEnd=ispet, strandSpecific=strandspec, ignoreDup=FALSE)

# Saving counts to file, with gene names.
colnames(outb$counts) <- sub("\\.bam$", "", basename(bam.files[(bad_index):length(bam.files)]))
final <- data.frame(GeneID=rownames(outb$counts), Length=outb$annotation$Length, outb$counts)
write.table(file="genic_counts_leaveDupb.tsv", final, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
save(outb, file="counts_leaveDupb.RData")

rm(final)

## concatenate
load('counts_leaveDupa.RData')
out <- outa
out$counts <- cbind(out$counts, outb$counts)
out$stat <- cbind(out$stat, outb$stat[,2:ncol(outb$stat)])
save(out, file="counts_leaveDup.RData")

# Augmenting the stats. ## put this in analysis script instead
##my.stats <- read.table(stat.file, header=TRUE)
##m <- match(my.stats$Sample, colnames(out$counts))
##my.stats$Genic <- as.integer(out$stat[1,-1][m])
##write.table(file="my_qual.tsv", my.stats, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")



# Saving the session information.
unlink("temp.gtf")

sink("sessionInfo.txt")
sessionInfo()
sink()

