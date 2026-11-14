suppressMessages({library(ggtree);library(ggplot2);library(ape)})
tr <- read.tree("data/mag_tree.treefile")
if("Fillmore" %in% tr$tip.label) tr <- root(tr, "Fillmore", resolve.root=TRUE, edgelabel=TRUE)
mag <- c("Mag1_15620328","Mag1_10513216")
LI <- c("9a5c","CVC0251","CVC0256","J1a12","XRB","50066","50080","11399","B111","U24D","Hib4","Pr8x","6c","COF0324","RAAR16_Co13","RAAR15_Co33","Fb7","ALM4","AM2-Angelina","OLI17A","32","3124","XF39","XF43","XF64.2","XF71_BR","XF-73-BR","XF143","XF217-2","XF221","XF45-A3","Reference","RAAR17_CiUb7","50088")
LII <- c("De_Donno","CoDiRO","APL2","APL9","APL12","APL18","APL31","APL32","APL37","APL42","APL50","APL57","APL64","APL69","APL71","APL72","APL75","APL76","MentOliv","MentPoly","COF0407","OLS0478","OLS0479","CFBP8072","CFBP8498","PD7211","PD7202","LSV47.09","LSV47.33","M5","M7","Salento-1","Salento-2")
g <- setNames(rep("other", length(tr$tip.label)), tr$tip.label)
g[intersect(LI,tr$tip.label)]  <- "Lineage I"
g[intersect(LII,tr$tip.label)] <- "Lineage II"
if("Fillmore" %in% tr$tip.label) g["Fillmore"] <- "outgroup"
g[mag] <- "Mag (metagenome)"
cols <- c("Lineage I"="#1B9E77","Lineage II"="#1F78B4","Mag (metagenome)"="#D62728","outgroup"="grey55","other"="grey45","0"="grey45")
trg <- groupOTU(tr, split(names(g), g))
dl <- data.frame(label=tr$tip.label, g=g[tr$tip.label], ismag=tr$tip.label %in% mag)

p <- ggtree(trg, aes(color=group), size=0.45) %<+% dl +
  scale_color_manual(values=cols, name=NULL, breaks=c("Lineage I","Lineage II","Mag (metagenome)","outgroup")) +
  geom_tiplab(aes(label=label, fontface=ifelse(ismag,"bold","plain"),
                  size=ifelse(ismag,3.2,2.1)), show.legend=FALSE) +
  scale_size_identity() +
  geom_tippoint(aes(subset=ismag), color="#D62728", size=2.4) +
  geom_treescale(width=0.01, fontsize=3) +
  theme(legend.position=c(0.12,0.85), legend.text=element_text(size=9)) +
  ggtitle("Placement of two plant-metagenome subsp. pauca assemblies (Mag1)") +
  theme(plot.title=element_text(size=11,face="bold")) + coord_cartesian(clip="off") +
  theme(plot.margin=margin(6,80,6,6))

ggsave("figures/FigureS3_mag.pdf", p, width=9, height=11, device=cairo_pdf)
ggsave("figures/FigureS3_mag.png", p, width=9, height=11, dpi=150, type="cairo")
cat("hecho\n")
