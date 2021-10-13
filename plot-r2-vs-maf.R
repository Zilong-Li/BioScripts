#!/usr/bin/env Rscript

getmaf <- function(x){
    ifelse(x > 0.5, 1-x, x)
}

plotInfoMaf <- function(dat, ytext, title) {
    colnames(dat) = c("bin", "info", "counts")
    barplot(dat$counts, main=title)
    par(new=T)
    plot(dat$bin * 0.01 - 0.01/2, dat$info, type='o', lwd=2, ylim=c(0.4, 1), axes=F, xlab="", ylab="", col="red", bty="n")
    mtext("MAF Bin", side=1, line=2.5)
    mtext("Counts", side=2, line=2.5)
    mtext(ytext, side=4, line=2.5)
    axis(side=1, at=bins)
    axis(side=4, at=seq(0, 1, 0.2), col="red")
}

args <- commandArgs(trailingOnly=T)
finfo <- args[1]
fout <- args[2]
imgsize <- 6
png(paste0(fout, ".png"), h=imgsize, w=imgsize*2, units="in", res=300)
header <- c("ID", "R2", "MAF")
par(family='arial', mar=c(4, 4, 4, 4), mfrow=c(1,2))

df <- data.table::fread(finfo, select=header)
bins = seq(0, 0.5, 0.01)
df$bin = cut(as.numeric(df$MAF), breaks=bins, labels=F)
# plot INFO vs MAF
dat = as.data.frame(matrix(unlist(parallel::mclapply(1:50, function(x){d = subset(df, bin==x); m = mean(as.numeric(d$R2)); n = length(d$R2); o = c(x, m, n)}, mc.cores=20)), ncol=3, byrow=T))
plotInfoMaf(dat, "R2", paste("R2 vs MAF \n#Sites: ", sum(dat$V3)))
# plot R2 >0.95 vs MAF
dat = as.data.frame(matrix(unlist(parallel::mclapply(1:50, function(x){d = subset(df, bin==x & R2>0.95); m = mean(as.numeric(d$R2)); n = length(d$R2); o = c(x, m, n)}, mc.cores=20)), ncol=3, byrow=T))
plotInfoMaf(dat, "R2", paste("R2(>0.95) vs MAF \n#Sites: ", sum(dat$V3)))

dev.off()

