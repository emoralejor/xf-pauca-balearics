suppressMessages({library(ggtree);library(treeio);library(ggplot2);library(ape);library(patchwork)})

## ====================== ARBOL (introgresion, 16.3 kb) ======================
TREE <- "data/bloques32.treefile"
tr <- read.tree(TREE)

ibiza     <- c("IXF31SC2","IAS-XYL1513-4","IAS-XYL1518","IVIA11529","XYL1961")
pauca     <- c("Hib4","CoDiRO","COF0407","OLS0479","9a5c","U24D","J1a12","CVC0251","XRB")
multiplex <- "Fillmore"
want <- c(ibiza,pauca,multiplex)
miss <- setdiff(want,tr$tip.label); if(length(miss)){cat("FALTAN:\n");print(miss);quit(save="no")}
fast <- setdiff(tr$tip.label, want)
tr <- root(tr, outgroup=multiplex, resolve.root=TRUE, edgelabel=TRUE)

disp <- tr$tip.label
disp[disp=="XYL2093_18"]<-"XYL2093/18"; disp[disp=="XYL2017_18"]<-"XYL2017/18"
grp <- setNames(rep("fast",length(tr$tip.label)), tr$tip.label)
grp[pauca]<-"pauca"; grp[multiplex]<-"mult"; grp[ibiza]<-"ibiza"
dlab <- data.frame(label=tr$tip.label, disp=disp,
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
  geom_treescale(width=0.02, x=0, y=-0.5, fontsize=3.2)

## corchetes de (sub)especie
o_ib  <- xr*0.05
o_sub <- xr*0.18
clbar <- function(n,off,bc) geom_cladelab(node=n,label="",offset=off,barcolor=bc,barsize=2.2,align=TRUE)
p <- p +
  clbar(fill_n, o_sub, grey) +
  clbar(node_pa, o_sub, green) +
  clbar(node_fa, o_sub, blue) +
  geom_cladelab(node=node_ib,label="Ibiza ST80",offset=o_ib,barcolor=orange,textcolor=orange,
                fontsize=3.0,barsize=2.2,angle=90,hjust=0.5,align=TRUE,offset.text=xr*0.015)

## textos de subespecie: "subsp." redonda + epiteto en cursiva
mx  <- max(p$data$x[p$data$isTip], na.rm=TRUE)
td2 <- p$data[p$data$isTip, ]
yfa <- mean(range(td2$y[td2$label %in% c(fast,ibiza)]))
ypa <- mean(range(td2$y[td2$label %in% pauca]))
ymu <- td2$y[td2$label==multiplex]
xt  <- mx + o_sub + xr*0.06
p <- p +
  annotate("text", x=xt, y=yfa, label="plain('subsp. ')*italic('fastidiosa')", parse=TRUE, hjust=0, size=3.0) +
  annotate("text", x=xt, y=ypa, label="plain('subsp. ')*italic('pauca')",      parse=TRUE, hjust=0, size=3.0) +
  annotate("text", x=xt, y=ymu, label="plain('subsp. ')*italic('multiplex')",  parse=TRUE, hjust=0, size=3.0)

## rango Y comun para alinear la pista
yall <- p$data$y
yr <- c(min(yall, -1) - 0.5, max(yall) + 0.8)
p <- p + coord_cartesian(ylim=yr, clip="off") +
  theme(plot.margin=margin(5, 112, 5, 5, unit="pt"))
td <- p$data[p$data$isTip, c("label","y")]

## ====================== PISTA fastGEAR (recent) ======================
D    <- "data/fastgear"
GLEN <- 2679306; GMB <- GLEN/1e6
lin <- read.table(file.path(D,"lineage_information.txt"), header=TRUE, stringsAsFactors=FALSE)
rec <- read.table(file.path(D,"recombinations_recent.txt"), skip=2, header=FALSE, stringsAsFactors=FALSE,
                  col.names=c("Start","End","DonorLineage","RecipientStrain","logBF","StrainName"), fill=TRUE)
rec <- rec[!is.na(rec$Start),]
Lname <- c("1"="multiplex","2"="pauca (South American)","3"="pauca ST80/ST53","4"="fastidiosa","5"="unknown")
Lcol  <- c("1"="grey55","2"="#1B9E77","3"="#1F78B4","4"="#D62728","5"="grey30")

rec$name <- lin$Name[match(rec$RecipientStrain, lin$StrainIndex)]
an <- intersect(unique(lin$Name), tr$tip.label)
cat(">> cepas fastGEAR en el arbol:", length(an), "/", length(unique(lin$Name)), "\n")
nomatch <- setdiff(unique(lin$Name), tr$tip.label)
if(length(nomatch)){ cat(">> [AVISO] nombres sin tip:\n"); print(nomatch) }

bg  <- data.frame(label=an, stringsAsFactors=FALSE)
bg$y <- td$y[match(bg$label, td$label)]
seg <- data.frame(label=rec$name, start=rec$Start/1e6, end=rec$End/1e6,
                  don=factor(rec$DonorLineage, levels=names(Lname)), stringsAsFactors=FALSE)
seg <- seg[seg$label %in% an, ]
seg$y <- td$y[match(seg$label, td$label)]

q <- ggplot() +
  geom_rect(data=bg,  aes(xmin=0, xmax=GMB, ymin=y-0.42, ymax=y+0.42), fill="grey94") +
  geom_rect(data=seg, aes(xmin=start, xmax=end, ymin=y-0.42, ymax=y+0.42, fill=don)) +
  scale_fill_manual(values=Lcol, breaks=names(Lname), labels=Lname, name="Donor lineage", drop=FALSE) +
  scale_x_continuous(name="Position in 9a5c chromosome (Mb)", limits=c(0,GMB),
                     expand=expansion(mult=c(0.01,0.01))) +
  coord_cartesian(ylim=yr, expand=FALSE, clip="off") +
  theme_minimal(base_size=9) +
  theme(panel.grid=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="bottom", plot.margin=margin(5,10,5,0))

## ====================== COMBINAR ======================
fig <- p + q + plot_layout(widths=c(0.82, 1.8))
PDF <- "figures/Figure3.pdf"
PNG <- "figures/Figure3.png"
ggsave(PDF, fig, width=16, height=8, device=cairo_pdf); cat("PDF:", PDF, "\n")

ok <- FALSE
try({ ggsave(PNG, fig, width=16, height=8, dpi=300, device=ragg::agg_png); ok <- TRUE }, silent=TRUE)
if(!ok) try({ ggsave(PNG, fig, width=16, height=8, dpi=300, type="cairo"); ok <- TRUE }, silent=TRUE)
if(!ok) try({ ggsave(PNG, fig, width=16, height=8, dpi=300); ok <- TRUE }, silent=TRUE)
if(!ok && requireNamespace("magick", quietly=TRUE))
  try({ magick::image_write(magick::image_read_pdf(PDF, density=300), PNG); ok <- TRUE }, silent=TRUE)
if(ok && file.exists(PNG)) cat("PNG:", PNG, "\n") else cat("[ERROR] PNG no guardado\n")
cat("hecho\n")
