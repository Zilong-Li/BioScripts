#!/usr/bin/env Rscript
## Usage: Rscript imputation-concordance.R truth-vcf list-of-imputed-vcf outdir bins

args <- commandArgs(trailingOnly = TRUE)
print(paste0(args, collapse = "', '"))
print(system("hostname"))

truth_vcf <- args[1]
list_vcf <- args[2]
outdir <- args[3]

bins <- sort(unique(c(
  c(0, 0.05 / 1),
  seq(0.1, 0.5, length.out = 5)
)))


library(data.table)

## input is matrix
r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE, per_snp = FALSE) {
  if (flip) {
    w <- af > 0.5
    af[w] <- 1 - af[w]
    truthG[w, ] <- 2 - truthG[w, ]
    testDS[w, ] <- 2 - testDS[w, ]
  }
  if (!is.null(which_snps)) {
    af <- af[which_snps]
    truthG <- truthG[which_snps, ]
    testDS <- testDS[which_snps, ]
  }
  truthG <- as.matrix(truthG)
  testDS <- as.matrix(testDS)
  x <- cut(af, breaks = breaks)
  if (ncol(truthG) > 1 && per_snp) {
    # for multiple sample, calculate r2 per snp then average them
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        simple = mean(sapply(w, function(ww) {
          cor(truthG[ww, ], testDS[ww, ], use = "pairwise.complete")**2
        }), na.rm = TRUE),
        norm = mean(sapply(w, function(ww) {
          cor(truthG[ww, ] - 2 * af[ww], testDS[ww, ] - 2 * af[ww], use = "pairwise.complete")**2
        }), na.rm = TRUE)
      )
    })
  } else {
    # calculate r2 based on a big matrix nsamples x nsnps
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        simple = cor(as.vector(truthG[w, ]), as.vector(testDS[w, ]), use = "pairwise.complete")**2,
        norm = cor(as.vector(truthG[w, ] - 2 * af[w]), as.vector(testDS[w, ] - 2 * af[w]), use = "pairwise.complete")**2
      )
    })
  }
  # fill with NA for AF bins without SNPs
  cors_per_af <- t(sapply(cors_per_af, function(a) {
    if (is.null(a[1])) {
      return(c(n = NA, nA = NA, simple = NA, norm = NA))
    }
    a
  }))
  return(cors_per_af)
}

## (gt0, gt1, ds) x nsamples
## system("bcftools query -f \"%CHROM:%POS[\t%GT]\n\" tmp/0000.vcf |sed -E 's/\\/|\\|/\t/g' > tmp/0000.truth.gt")
## system("bcftools query -f \"%CHROM:%POS[\t%GT\t%DS]\n\" tmp/0003.vcf |sed -E 's/\\/|\\|/\t/g' > tmp/0003.stitch")
parse.imputed.gts2 <- function(fn) {
  d1 <- fread(fn, data.table = F)
  id <- d1[, 1]
  d1 <- as.matrix(suppressWarnings(sapply(d1[, -1], as.numeric)))
  rownames(d1) <- id
  return(d1)
}

## get AF and truth genotypes
af_file <- file.path(outdir, "af.tsv")
system(paste("bcftools +fill-tags", truth_vcf,"-- -t AC,AN,AF | bcftools query -f \"%CHROM:%POS\t%AF\n\" >", af_file))
truth_file <- file.path(outdir, "truth.tsv")
system(paste("bcftools query -f \"%CHROM:%POS[\t%GT]\n\"", truth_vcf,"| sed -E 's/\\/|\\|/\t/g' >", truth_file))

truth <- fread(truth_file, data.table=F)
truth.id <- truth[,1]
truth <- sapply(seq(1, dim(truth)[2] - 1, 2), function(i) {
                    rowSums(truth[, (i + 1):(i + 2)])
        }) # matrix: nsnps x nsamples
rownames(truth) <- truth.id

d.af <- read.table(af_file)
af <- as.numeric(d.af[, 2])
names(af) <- d.af[, 1]
rm(d.af)
af <- ifelse(af > 0.5, 1-af, af) # make MAF

## get list of imputed genotypes
ff <- read.table(list_vcf, h = F)
colnames(ff) <- c("id", "vcf")
imputed <- lapply(1:nrow(ff), function(i) {
  imputed_vcf <- ff[i,2]
  id <- ff[i,1]
  out <- file.path(outdir, paste0(id, ".tsv"))
  system(paste("bcftools query -f \"%CHROM:%POS[\t%GT\t%DS]\n\" ", imputed_vcf,"| sed -E 's/\\/|\\|/\t/g' >", out))
  parse.imputed.gts2(out)
})
names(imputed) <- ff$id

acc <- lapply(imputed, function(test) {
  iid <- intersect(rownames(truth), rownames(test))
  test <- test[,seq(3, ncol(test), by = 3 )] # get dosage
  o <- r2_by_freq(bins, af, truth, test, which_snps=iid, flip=F, per_snp=F)
  o[,c("n","simple")]
})

acc <- do.call("cbind",acc)
rm <- seq(1, ncol(acc), 2)[-1]
acc <- acc[, -rm]
colnames(acc) <- c("nsites", ff$id)
print(acc)
