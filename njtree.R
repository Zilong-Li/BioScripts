#!/home/zilong/miniconda3/envs/R4/bin/Rscript
library(ape)
library(tidyverse)
library(ggtree)
library(wesanderson)
library(RColorBrewer)
library(ggpubr)
library(ggsci)

info <- read.table("data/ibsmat.pops", header = T)
mat <- as.matrix(read.table("data/ibsmat.txt", header = F))
colnames(mat) <- info$ID
rownames(mat) <- info$ID

# unrooted tree
t1 <- mat %>% nj()

# pick outgroup 557 Bohor_reedbuck
t2 <- mat %>%
  nj() %>%
  root(outgroup = "557", resolve.root = T)

lpps <- sapply(unique(info$Pop), function(x) {
  info[info$Pop == x, ]$ID
})


mycol <- wes_palette("Darjeeling1", length(lpps), type = "continuous")
mycol <- colorRampPalette(brewer.pal(8, "Dark2"))(length(lpps))

# allows up up to 20 different populations, kelly palette from package catecolors https://github.com/btupper/catecolors
cols <- c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A", "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16")
mycol <- cols[1:length(lpps)]
names(mycol) <- names(lpps)

p1 <- t1 %>%
  ggtree(layout = "circular", branch.length = "none") %>%
  groupOTU(lpps, "Pops") + aes(color = Pops) + scale_color_manual(values = mycol) + theme(legend.position = "right") + geom_tiplab(size = 3, show.legend = FALSE)
## print(p1)

p2 <- t2 %>%
  ggtree(layout = "circular", branch.length = "none") %>%
  groupOTU(lpps, "Pops") + aes(color = Pops) + scale_color_manual(values = mycol) + theme(legend.position = "right") + geom_tiplab(size = 3, show.legend = FALSE)
## print(p2)

p <- ggarrange(p1, p2, nrow = 1, ncol = 2, labels = "AUTO")
## print(p)

ggsave("njtree-circular.png", p2, h = 7, w = 7)
ggsave("njtree-circular.pdf", p, h = 7, w = 14)
