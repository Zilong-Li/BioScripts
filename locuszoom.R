#!/usr/bin/env Rscript


locusZoomNoLD <- function(pval, pos, chr, refGenes = NULL, w = 1, plog = FALSE, main = "LocusZoom") {
  xforrs <- 0.03
  regsz <- 1.2
  width <- 22
  rsqplus <- 0.045
  rightylabxplus <- 0.05
  xmargin <- 0.005
  cexaxis <- 1.3
  cexlab <- 1.5
  blue <- "dodgerblue4"

  pval[is.na(pval)] <- 1

  N <- length(pval)

  rsqwithn <- rep(0, N - 1)
  ## palette(brewer.pal(10, "RdBu"))[1:10]
  cols <- c(10, 8, 6, 4, 2, 2)
  getcol <- function(x) {
    cuts <- c(.2, .4, .6, .8, 1)
    cols[6 - sum(as.numeric(x) < cuts)]
  }

  min.pos <- min(pos) / 1e6
  max.pos <- max(pos) / 1e6

  adj <- 0
  nf <- layout(matrix(c(1, 2), 2, 1, byrow = TRUE), heights = c(4.0, 2.5), widths = c(9))
  par(mar = c(0, 5.5, 2, 1))

  bgcols <- sapply(rsqwithn, getcol)
  bgcols[w] <- 1
  if (plog) {
    ylim <- c(0, max(10.5, max(-log10(pval))))
    pval <- -log10(pval)
  } else {
    ylim <- c(0, max(pval)*1.05)
  }
  plot(pos / 1e6, pval, col = "black", bg = bgcols, xaxs = "i", xlim = c(min.pos - adj, max.pos + adj), xaxt = "n", xlab = "", ylab = "", ylim = ylim, pch = 21, cex = 2, las = 1, cex.axis = cexaxis, main = main)
  ## axis(2, ylim = c(0, 15), col = "black", label = FALSE, cex.axis = cexaxis)
  ## mtext(2, text = expression(paste("-l", og[10], "[", italic(P), "]")), line = 3, cex = cexlab)
  mtext(2, text = "SNP Loadings of PC1", line = 4, cex = cexlab)
  xstart <- max.pos - (max.pos - min.pos) * 0.15
  scalewidth <- (max.pos - min.pos) * 0.0125
  xx <- c(xstart, xstart, xstart + scalewidth, xstart + scalewidth)
  #  legend("topleft",legend=c(rs),pt.bg=1,col="black",cex=1.6,pch=21)

  ybot <- ylim[2] * 0.65
  ytop <- ylim[2] * 0.95
  ysz <- ytop - ybot
  cuts <- c(.2, .4, .6, .8, 1)
  txtcuts <- c("0.2", "0.4", "0.6", "0.8", "1.0")
  txtposs <- c(0.2, 0.4, 0.6, 0.8, 1.00)
  scalexplus <- (max.pos - min.pos) * 0.06
  ## if (FALSE) {
  ##   polygon(xx + scalexplus, c(ybot, ybot + cuts[1] * ysz, ybot + cuts[1] * ysz, ybot), col = cols[1])

  ##   for (i in 2:length(cols)) {
  ##     polygon(xx + scalexplus, c(ybot + cuts[(i - 1)] * ysz, ybot + cuts[i] * ysz, ybot + cuts[i] * ysz, ybot + ysz * cuts[(i - 1)]), col = cols[i])
  ##   }
  ##   scalenumplus <- (max.pos - min.pos) * 0.09
  ##   # text((xx[1]+0.4+xx[3])/2+rsqplus,ytop+0.35,expression(italic(r)^2),cex=cexaxis)
  ##   text(xx[3] + scalenumplus, ytop + 0.35, expression(italic(r)^2), cex = cexaxis)

  ##   text(rep(xx[3] + scalenumplus, 2), txtposs[3:4] * ysz + ybot - 0.07, txtcuts[3:4], cex = cexlab - 0.1)
  ##   text(rep(xx[3] + scalenumplus, 1), txtposs[2] * ysz + ybot - 0.07, txtcuts[2], cex = cexlab - 0.1)
  ##   text(rep(xx[3] + scalenumplus, 1), txtposs[1] * ysz + ybot - 0.07, txtcuts[1], cex = cexlab - 0.1)
  ## }

  dat <- read.table(refGenes, as.is = T, head = T, comment.char = "")
  xx2 <- dat[dat[, "chrom"] == paste("chr", chr, sep = "") & dat[, "cdsStart"] < max.pos * 1e6 & dat[, "cdsEnd"] > min.pos * 1e6, ]

  start <- xx2$txStart # casia. this column is needed in the refgene file
  end <- xx2$txEnd # casia. this column is needed in the refgene file
  nams <- xx2$name2 # casia. this column is needed in the refgene file
  cnts <- xx2$exonCount # casia. this column is needed in the refgene file but you can have all as NA

  abline(v = pos, lty = 2)
  ## par(mar = c(5.2, 6.2, -0.1, 6.3) + 0.1)
  par(mar = c(5, 5.5, 1, 1))

  plot(c(0, 0), c(0, 0), type = "n", xlim = c(min.pos - adj, max.pos + adj), ylim = c(-0.8, 0.1), xlab = "", xaxs = "i", yaxt = "n", ylab = "", main = "", cex.lab = 2.6, cex.axis = cexaxis, tck = -0.05)
  mtext(1, text = paste("Position on chromosome ", chr, " (Mb)", sep = ""), line = 3, cex = cexlab)
  abline(v = pos, lty = 2)

  ord <- order(start)
  start <- start[ord]
  end <- end[ord]
  exoncnts <- cnts[ord]
  nams <- nams[ord]
  keep <- !duplicated(nams)
  start <- start[keep]
  end <- end[keep]
  exoncnts <- cnts[keep]
  nams <- nams[keep]
  ord <- ord[keep]
  he <- rep(c(0, -0.18, -0.36, -0.54, -0.72), 100)[1:length(nams)] - 0.05
  if (length(start) > 0) {
    segments(start / 1e6, he, end / 1e6, he)
    keep <- !duplicated(nams)
    sapply(1:sum(keep), function(x) {
      text((end[keep][x] + start[keep][x]) / 2e6, he[keep][x] + 0.08, bquote(italic(.(nams[keep][x]))), cex = cexlab - 0.6)
    })
    estart <- as.numeric(unlist(sapply(xx2$exonStarts[ord], function(y) {
      strsplit(y, ",")[[1]]
    }))) / 1e6 # Casia outcomment
    eend <- as.numeric(unlist(sapply(xx2$exonEnds[ord], function(y) {
      strsplit(y, ",")[[1]]
    }))) / 1e6 # Casia outcomment
    rect(estart, rep(he, xx2$exonCount[ord]) - 0.01, eend, rep(he, xx2$exonCount[ord]) + 0.01, col = "black") # Casia outcomment
  }
}

refGenes <- "~/Downloads/refGene.hg38.txt.gz"
locusZoomNoLD(t$pc, t$pos, 5, refGenes = refGenes)


d <- read.table("t.loadings.chr5.pc1.csv", h = F)
colnames(d) <- c("chr", "pos", "pc")
t <- subset(d, pos > 33851693 & pos < 33991693)

head(sort(t$pc, de = T))

dev.off()

pos <- t$pos
min.pos <- min(pos) / 1e6
max.pos <- max(pos) / 1e6
head(xx2)
