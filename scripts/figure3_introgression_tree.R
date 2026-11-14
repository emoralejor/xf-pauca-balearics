suppressMessages({library(ggtree);library(treeio);library(ggplot2);library(ape)})
TREE <- "data/bloques32.treefile"
tr <- read.tree(TREE)

ibiza     <- c("IXF31SC2","IAS-XYL1513-4","IAS-XYL1518","IVIA11529","XYL1961")
pauca     <- c("Hib4","CoDiRO","COF0407","OLS0479","9a5c","U24D","J1a12","CVC0251","XRB")
multiplex <- "Fillmore"
want <- c(ibiza,pauca,multiplex)
miss <- setdiff(want,tr$tip.label); if(length(miss)){cat("FALTAN:\n");print(miss);quit(save="no")}
fast <- setdiff(tr$tip.label, want)
cat("fastidiosa n =",length(fast),"\n")
tr <- root(tr, outgroup=multiplex, resolve.root=TRUE, edgelabel=TRUE)

disp <- tr$tip.label
disp[disp=="XYL2093_18"]<-"XYL2093/18"; disp[disp=="XYL2017_18"]<-"XYL2017/18"
info <- setNames(rep(NA_character_,length(tr$tip.label)), tr$tip.label)
info["CoDiRO"]<-"Apulia"; info["9a5c"]<-"reference"
info["XYL2093_18"]<-"Mallorca, ST1"; info["XYL1732"]<-"Mallorca, ST1"; info["XYL2017_18"]<-"Mallorca, ST1"
info["M23"]<-"California"; info["Temecula1"]<-"California"
grp <- setNames(rep("fast",length(tr$tip.label)), tr$tip.label)
grp[pauca]<-"pauca"; grp[multiplex]<-"mult"; grp[ibiza]<-"ibiza"
dlab <- data.frame(label=tr$tip.label, disp=disp, info=info[tr$tip.label],
                   grp=grp[tr$tip.label], face=ifelse(tr$tip.label %in% ibiza,"bold","plain"),
                   stringsAsFactors=FALSE)

trg <- groupOTU(tr, split(names(grp),grp))
orange<-"#E8820C"; green<-"#1B9E77"; blue<-"#1F78B4"; grey<-"grey45"
cols <- c("0"="grey30", fast=blue, pauca=green, mult=grey, ibiza=orange)
node_ib <- getMRCA(tr,ibiza); node_pa <- getMRCA(tr,pauca); node_fa <- getMRCA(tr,c(fast,ibiza))
node_circle <- tr$edge[tr$edge[,2]==node_ib,1]
fill_n <- which(tr$tip.label==multiplex)
xr <- max(ggtree(tr)$data$x[ggtree(tr)$data$isTip], na.rm=TRUE)

p <- ggtree(trg, aes(color=group), size=0.55) %<+% dlab +
  scale_color_manual(values=cols, guide="none") +
  geom_tippoint(size=1.5) +
  geom_tiplab(aes(label=disp, fontface=face), offset=xr*0.008, size=2.9) +
  geom_text2(aes(subset=(!isTip & label!="" & !is.na(label) & branch.length>=1e-5),
                 label=label, fontface=ifelse(suppressWarnings(as.numeric(label))>=95,"bold","plain")),
             size=2.4, color="grey30", hjust=1.2, vjust=-0.4) +
  geom_point2(aes(subset=(node==node_circle)), shape=21, fill="white", color=orange, size=3.2, stroke=1.3) +
  geom_treescale(width=0.02, x=0, y=-0.5, fontsize=3.2) +
  theme(plot.margin=margin(5,95,5,5)) + coord_cartesian(clip="off")

## --- lideres de puntos manuales: empiezan TRAS la ultima letra del nombre ---
d <- p$data; maxx <- max(d$x[d$isTip], na.rm=TRUE)
CW <- 0.00085                 # dial: ancho medio por caracter
host_x <- maxx + xr*0.12
tips <- d[d$isTip, c("label","x","y")]
tips$disp <- dlab$disp[match(tips$label, dlab$label)]
tips$info <- dlab$info[match(tips$label, dlab$label)]
tips$x0 <- tips$x + xr*0.008 + nchar(tips$disp)*CW + xr*0.015
tips$x0 <- pmin(tips$x0, host_x - xr*0.02)
p <- p +
  geom_segment(data=tips, aes(x=x0, xend=host_x-xr*0.02, y=y, yend=y),
               linetype="dotted", linewidth=0.2, color="grey70", inherit.aes=FALSE) +
  geom_text(data=tips[!is.na(tips$info),], aes(x=host_x, y=y, label=info),
            inherit.aes=FALSE, hjust=0, size=2.4, color="grey55")

## --- corchetes de (sub)especie ---
o1 <- xr*0.30; o2 <- xr*0.60
cl <- function(n,lab,off,bc,tc) geom_cladelab(node=n,label=lab,offset=off,barcolor=bc,textcolor=tc,
                                              fontsize=3.1,barsize=2.2,angle=0,hjust=0,align=TRUE,offset.text=xr*0.02)
p <- p +
  cl(fill_n,"subsp. multiplex", o2, grey,  "black") +
  cl(node_pa,"subsp. pauca",    o2, green, "black") +
  cl(node_fa,"subsp. fastidiosa",o2, blue,  "black") +
  cl(node_ib,"Ibiza,\nsubsp. pauca\nST80", o1, orange, orange)

ggsave("figures/Figure3.pdf", p, width=9.5, height=8, device=cairo_pdf)
ggsave("figures/Figure3.png", p, width=9.5, height=8, dpi=150, type="cairo")
cat("hecho\n")
