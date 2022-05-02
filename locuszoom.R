#!/usr/bin/env Rscript


colpal <- rev(c("#b2182b", "#ef8a62", "#fddbc7", "#d1e5f0", "#67a9cf", "#2166ac"))

locusZoomNoLD <- function(pval, w = 1, pos, chr, main = "LocusZoom") {
  geneticMap <- "/home/albrecht/geneticMap/hg19/genetic_map_GRCh37_chr" # Casia. You can remove genetic map
  refGenes <- "/home/albrecht/geneticMap/hg19/refGeneHG19.gz"

  xforrs <- 0.03
  regsz <- 1.2
  width <- 22

  rsqplus <- 0.045
  rightylabxplus <- 0.05
  xmargin <- 0.005
  cexaxis <- 1.5
  cexlab <- 1.5
  blue <- "dodgerblue4"

  pval[is.na(pval)] <- 1

  N <- length(pval)

  rsqwithn <- rep(0, N - 1)
  #  palette(brewer.pal(10, "RdBu"))[1:10]
  cols <- c(10, 8, 6, 4, 2, 2)
  getcol <- function(x) {
    cuts <- c(.2, .4, .6, .8, 1)
    cols[6 - sum(as.numeric(x) < cuts)]
  }

  min.pos <- min(pos) / 1e6
  max.pos <- max(pos) / 1e6
  rec <- read.table(paste(geneticMap, chr, ".txt", sep = ""), header = T)
  relrec <- rec[(rec[, 2] / 1e6) > min.pos & (rec[, 2] / 1e6) < max.pos, ]

  adj <- 0
  nf <- layout(matrix(c(1, 2), 2, 1, byrow = TRUE), heights = c(5.5, 2.65), widths = c(8))
  par(mar = c(0.2, 6.2, 4, 6.3) + 0.1)
  plot(relrec[, "Position.bp."] / 1e6, relrec[, "Rate.cM.Mb."], type = "l", ylim = c(0, 105), xaxs = "i", xlim = c(min.pos - adj, max.pos + adj), axes = F, col = blue, xlab = "", ylab = "", lwd = 3.5, cex.axis = cexaxis)
  axis(4, labels = seq(0, 100, 20), at = seq(0, 100, 20), col = blue, las = 1, col.axis = blue, cex.axis = cexaxis)
  #  text(par("usr")[2]+rightylabxplus,45,xpd=TRUE,srt=-90,labels=c(expression(paste("Recombination rate (cM ", Mb^-1,")"))),col=blue,cex=cexlab)
  mtext(c(expression(paste("Recombination rate (cM ", Mb^-1, ")"))), 4, 4, col = blue, cex = cexlab)

  par(new = T)
  bgcols <- sapply(rsqwithn, getcol)
  bgcols[w] <- 1
  ylim <- c(0, max(10.5, max(-log10(pval))))
  plot(pos / 1e6, -log10(pval), col = "black", bg = bgcols, xaxs = "i", xlim = c(min.pos - adj, max.pos + adj), xaxt = "n", xlab = "", ylab = "", ylim = ylim, pch = 21, cex = 2, las = 1, cex.axis = cexaxis, main = main)
  axis(2, ylim = c(0, 15), col = "black", label = FALSE, cex.axis = cexaxis)
  mtext(2, text = expression(paste("-l", og[10], "[", italic(P), "]")), line = 3, cex = cexlab)
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
  if (FALSE) {
    polygon(xx + scalexplus, c(ybot, ybot + cuts[1] * ysz, ybot + cuts[1] * ysz, ybot), col = cols[1])

    for (i in 2:length(cols)) {
      polygon(xx + scalexplus, c(ybot + cuts[(i - 1)] * ysz, ybot + cuts[i] * ysz, ybot + cuts[i] * ysz, ybot + ysz * cuts[(i - 1)]), col = cols[i])
    }
    scalenumplus <- (max.pos - min.pos) * 0.09
    # text((xx[1]+0.4+xx[3])/2+rsqplus,ytop+0.35,expression(italic(r)^2),cex=cexaxis)
    text(xx[3] + scalenumplus, ytop + 0.35, expression(italic(r)^2), cex = cexaxis)

    text(rep(xx[3] + scalenumplus, 2), txtposs[3:4] * ysz + ybot - 0.07, txtcuts[3:4], cex = cexlab - 0.1)
    text(rep(xx[3] + scalenumplus, 1), txtposs[2] * ysz + ybot - 0.07, txtcuts[2], cex = cexlab - 0.1)
    text(rep(xx[3] + scalenumplus, 1), txtposs[1] * ysz + ybot - 0.07, txtcuts[1], cex = cexlab - 0.1)
  }

  dat <- read.table(refGenes, as.is = T, head = T, comment.char = "")
  xx2 <- dat[dat[, "chrom"] == paste("chr", chr, sep = "") & dat[, "cdsStart"] < max.pos * 1e6 & dat[, "cdsEnd"] > min.pos * 1e6, ]

  start <- xx2$txStart # casia. this column is needed in the refgene file
  end <- xx2$txEnd # casia. this column is needed in the refgene file
  nams <- xx2$name2 # casia. this column is needed in the refgene file
  cnts <- xx2$exonCount # casia. this column is needed in the refgene file but you can have all as NA

  abline(v = pos, lty = 2)
  par(mar = c(5.2, 6.2, -0.1, 6.3) + 0.1)


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


manPlot <- function(x, chr, thinLarge = TRUE, pass, cap = 1e-30, collar = c("darkblue", "#67a9cf", "orange", "grey"), main = "") {
  keep <- !is.na(x)
  x <- x[keep]
  chr <- chr[keep]
  N <- length(x)
  bonf <- -log(0.05 / N, base = 10)

  if (thinLarge & N > 1.1e5) {
    num <- 10
    if (N > 1e6) {
      num <- 25
    }
    if (N > 5e6) {
      num <- 50
    }
    n1 <- round(N / num)
    k <- c(1:n1, 1:n1 * num)
    k <- sort(order(x)[k])
    x <- x[k]

    chr <- chr[k]
    if (!missing(pass)) {
      pass <- pass[k]
    }
  }


  x <- ifelse(x < cap, cap, x)
  col <- chr %% 2 + 1
  pch <- rep(16, length(x))
  if (!missing(pass)) {
    col[!pass[keep]] <- 4
    pch[!pass[keep]] <- 1
  }
  #     par( mar = c( 5.1, 5.1, 4.1, 1.1 ) )
  maxy <- max((-log(x, base = 10)))
  plot(-log(x, base = 10), col = collar[col], ylab = expression(-log[10](italic(P))), xlab = "Chromosomes", main = main, cex = 1, lwd = 2, pch = pch, axes = F, cex.lab = 1, cex.main = 2, ylim = c(0, maxy + 0.05 * maxy))
  box()
  axis(2, las = 1, cex.axis = 1.8)
  abline(h = bonf, lty = 2, col = 1)
  ## legend( 0, t-0.5, "Bonferroni correction", lty = 2, bty = "n", col = 1, cex = 2 )
  ## legend( "topright", c( "known", "novel" ), pch = c( 4, 1 ), bty = "n", cex = 2 )


  if (maxy > -log10(cap) * .99) {
    text(0, -log10(cap) * 1.03, "capped", adj = 0)
    abline(h = -log10(cap), lty = 2, col = "grey")
  }
}

qqp <-
  function(x, ci = TRUE, add = FALSE, ylab = "Observed log10(p-value)", xlab = "Expected log10(p-value)", maxLogP, thinLarge = TRUE, col = 1, ...) {
    if (length(col) > 1) {
      col <- col[!is.na(x)]
    } # NB ida changed this from col<-col[!is.a(x)]
    x <- x[!is.na(x)]
    if (!missing(maxLogP)) {
      x[x < 10^-maxLogP] <- 10^-maxLogP
    }
    N <- length(x)
    ord <- order(x)
    x <- x[ord]
    if (length(col) > 1) {
      col <- col[ord]
    }
    ## lambda<-round(median(x)/qchisq(0.5,1),2)

    if (thinLarge & N > 1.1e5) {
      n1 <- round(N / 100)
      keep <- c(1:n1, 1:n1 * 100)
      e <- -log((1:N - 0.5)[keep] / N, 10)
      x <- x[keep]
    } else {
      e <- -log((1:N - 0.5) / N, 10)
    }

    if (add) {
      points(e, -log(x, 10), col = col, ...)
    } else {
      plot(e, -log(x, 10), ylab = ylab, xlab = xlab, col = col, ...)
      abline(0, 1, col = 2, lwd = 2)
    }
    # legend("topleft",paste("lambda=",lambda))

    if (ci) {
      # https://www.tandfonline.com/doi/abs/10.1080/00949658008810388https://www.tandfonline.com/doi/abs/10.1080/00949658008810388
      # concentration bands
      if (thinLarge & N > 1.1e5) {
        c97.5 <- qbeta(0.975, keep, N - (keep) + 1)
        c02.5 <- qbeta(0.025, keep, N - (keep) + 1)
      } else {
        c97.5 <- qbeta(0.975, 1:N, N - (1:N) + 1)
        c02.5 <- qbeta(0.025, 1:N, N - (1:N) + 1)
      }
      lines(e, -log(c97.5, 10))
      lines(e, -log(c02.5, 10))
    }
  }

qqPlot <- function(x, cap = 1e-30, main = "") {
  keep <- !is.na(x)
  x <- x[keep]
  x <- ifelse(x < cap, cap, x)
  maxy <- max((-log(x, base = 10)))

  qqp(x,
    pch = 16, col = "darkblue", main = main, las = 1, cex.lab = 1, cex.main = 2, ylim = c(0, maxy + 0.05 * maxy),
    xlab = expression(Expected ~ ~ -log[10](italic(P))), ylab = expression(Observed ~ ~ -log[10](italic(P))), cex.axis = 1.8
  )

  if (maxy > -log10(cap) * .99) {
    text(0, -log10(cap) * 1.03, "capped", adj = 0)
    abline(h = -log10(cap), lty = 2, col = "grey")
  }
}



plink.logistic <- function(x, pval) {
  p <- read.table(x, colC = c("integer", "character", "integer", "NULL", "character", "integer", "numeric", "NULL", "numeric"), head = T)
  p <- p[p$TEST == "ADD", ]
  if (!missing(pval)) {
    print(subset(p, P < pval))
  }

  manPlot(p$P, p$CHR, main = basename(x))
  qqPlot(p$P, main = basename(x))
}

plink.genome <- function(x) {
  d <- read.table(x, hea = T, colC = c("character", "character", "character", "character", "character", rep("numeric", 14 - 5)))


  plot(d$Z1, d$Z2, ylab = "k2", xlab = "k1", main = " pairwise relatedness", col = "darkblue", lwd = 2, ylim = c(0, 1), xlim = c(0, 1))

  text(0, 1, "MZ", font = 2)
  text(1, 0, "PO", font = 2)
  text(0.5, 0.25, "FS", font = 2)
  text(0.5, 0, "HS", font = 2)
  text(0.25, 0, "C1", font = 2)
  text(1 / 16, 0, "C2", font = 2)
  text(0, 0, "UR", font = 2)

  legend("topright", c("MZ: Monozygotic Twins", "PO: Parent offspring", "Full siblings", "HS: halfsibling", "C1: first cousins", "UR: unrelated"))



  if (any(d$Z1 > 0.2)) {
    o <- subset(d, Z1 > 0.2)
    ord <- order(o$PI_HAT, decreasing = T)
    return(o[ord, , drop = F])
  }
}

plink.mds <- function(x, fam) {
  d <- read.table(x, head = T)
  col <- "darkblue"

  ccpal <- colpal[c(1, 6, 2)]
  if (!missing(fam)) {
    f <- read.table(fam, head = F, as.is = T)
    status <- as.factor(f[, 6])
    col <- ccpal[status]
  }

  main <- paste("MDS", basename(x))
  plot(d$C1, d$C2, xlab = "dim 1", ylab = "dim 2", main = main, lwd = 2, col = col)
  if (!missing(fam)) {
    legend("top", legend = levels(status), fill = ccpal[1:length(levels(status))], hor = T, title = "Disease Status")
  }
}

plotPlink <- function(x, ...) {
  type <- sapply(strsplit(x, ".", fixed = T), function(y) y[length(y)])

  if (type == "logistic") {
    try(plink.logistic(x, ...))
  }
  if (type == "genome") {
    try(plink.genome(x))
  }
  if (type == "mds") {
    try(plink.mds(x, ...))
  }
}
