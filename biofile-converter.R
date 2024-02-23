library(genio)

## convert phase file to plink bed files
## phase file format
#  0
#  ninds
#  nsites
#  P bp1 bp2 bp3 ...
#  SSSSSSSS
#  mat: 2 x ninds x nsites


convert.phase.plink <- function(phasefiles) {
  for(chr in names(phasefiles)){
    phase <- readLines(phasefiles[chr])
    nind <- as.integer(phase[2])
    nsites <- as.integer(phase[3])
    pos <- as.integer(unlist(strsplit(phase[4], split = " "))[-1])
    mat <- phase[6:length(phase)]

    gt <- as.integer(unlist(strsplit(mat[1], split="")))
    stopifnot(length(pos)==length(gt))

    bim <- make_bim( n = nsites )
    bim$chr <- chr
    bim$pos <- pos
    bim$id <- paste0(chr, "_", pos)

    fam <- make_fam( n = nind )
    G <- sapply(1:length(mat), function(i){
      gt <- as.integer(unlist(strsplit(mat[i], split="")))
    })

    G <- G[, seq(1, ncol(G), 2)] + G[, seq(2, ncol(G), 2)]
    fileplink <- paste0("example.", chr)
    o <- write_plink(fileplink, G, bim, fam, verbose = FALSE)
  }
}

phasefiles <- c("chr1"= "data/example.phase", "chr2"="data/example.phase")

convert.phase.plink(phasefiles)
