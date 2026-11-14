suppressMessages({library(ggtree);library(treeio);library(ggplot2);library(ape)})
tr <- read.tree("data/pauca69_ML.treefile")
tr <- root(tr, outgroup="Fillmore", resolve.root=TRUE, edgelabel=TRUE)
ufb <- sapply(strsplit(tr$node.label,"/"), function(z){v<-suppressWarnings(as.numeric(z)); if(length(v)>=2) v[2] else NA_real_})
N<-length(tr$tip.label); phy2<-reorder(tr,"postorder"); cnt<-integer(N+tr$Nnode); cnt[1:N]<-1
for(i in 1:nrow(phy2$edge)) cnt[phy2$edge[i,1]]<-cnt[phy2$edge[i,1]]+cnt[phy2$edge[i,2]]
dsup <- data.frame(node=(N+1):(N+tr$Nnode), ufb=ufb, ntips=cnt[(N+1):(N+tr$Nnode)])
fi<-which(tr$tip.label=="Fillmore"); e<-which(tr$edge[,2]==fi); tr$edge.length[e]<-0.06

meta <- read.delim("data/meta.tsv", stringsAsFactors=FALSE)
norm <- function(s) gsub("[^a-z0-9]","", tolower(s))
k <- norm(tr$tip.label); k[tr$tip.label=="Reference"]<-"9a5c"
dinfo <- data.frame(label=tr$tip.label, info=meta$info[match(k,meta$key)], stringsAsFactors=FALSE)

ib   <- c("XYL1961","IXF31SC2","IAS-XYL1513-4","IAS-XYL1518","IVIA11529")
linI <- c("Hib4","RAAR16_Co13","Pr8x","6c","COF0324","11399","XRB","50080","50088","J1a12",
          "B111","50066","CVC0251","CVC0256","U24D","RAAR17_CiUb7","Reference","ALM4","Fb7",
          "AM2-Angelina","OLI17A","32","RAAR15_Co33","XF43","3124","XF64.2","XF45-A3",
          "XF71_BR","XF-73-BR","XF143","XF39","XF217-2","XF221")
st53 <- c("APL37","APL31","De_Donno","APL42","APL50","APL12","APL2","MentPoly","MentOliv",
          "COF0407","OLS0479","OLS0478","APL32","APL9","CoDiRO","APL75","APL18","APL76",
          "APL71","APL57","M7","M5","APL64","APL72","APL69","PD7202","LSV47.33")
rel  <- c("CFBP8072","CFBP8498","PD7211")
node_ib<-getMRCA(tr,ib); node_I<-getMRCA(tr,linI); node_ST53<-getMRCA(tr,st53)
node_pauca<-getMRCA(tr, setdiff(tr$tip.label,"Fillmore"))
st53_xmax <- max(ggtree(tr)$data$x[ggtree(tr)$data$label %in% st53], na.rm=TRUE)
trg<-groupClade(tr,node_ib)
orange<-"#E8820C"; green<-"#1B9E77"; blue<-"#1F78B4"
name_sz<-3.0; host_sz<-2.9; num_sz<-2.9; star_sz<-4.3; TS<-3.9

## arbol + nombres de cepa (SIN la columna de huesped todavia)
p <- ggtree(trg, aes(color=group), size=0.45) %<+% dsup +
  scale_color_manual(values=c("black",orange), guide="none") +
  geom_tiplab(size=name_sz, offset=0.0015) +
  geom_text2(aes(subset=(!isTip & !is.na(ufb) & ntips>=4 & branch.length>=3e-4 & ufb>=99), label="✺"),
             color="black", size=star_sz, vjust=0.35) +
  geom_text2(aes(subset=(!isTip & !is.na(ufb) & branch.length>=1e-5 &
                         ((ntips>=4 & ufb>=90 & ufb<99) | node==node_pauca)), label=ufb),
             color="grey35", size=num_sz, hjust=1.1, vjust=-0.4) +
  geom_treescale(width=0.01, x=0, y=-1, fontsize=3.6) +
  theme(plot.margin=margin(5,55,5,5)) + coord_cartesian(clip="off")

p <- scaleClade(p, node=node_ST53, scale=0.6)
p <- collapse(p, node=node_ST53, mode="mixed", fill=blue, alpha=0.55, color=blue)
y_st53 <- p$data$y[p$data$node==node_ST53]
p <- p + annotate("text", x=st53_xmax+0.004, y=y_st53, label="ST53 (n=27)", color=blue, size=3.8, hjust=0, vjust=0.5)

## ---- columna de huesped con lideres de puntos que empiezan TRAS la ultima letra ----
d <- p$data
maxx <- max(d$x[d$isTip], na.rm=TRUE)
CW      <- 0.00135          # dial: ancho medio por caracter (nombre de cepa)
host_x  <- maxx + 0.020     # dial: inicio de la columna de huesped
tips <- d[d$isTip, c("label","x","y")]
tips$info <- dinfo$info[match(tips$label, dinfo$label)]
tips <- tips[!is.na(tips$info), ]
tips$x0 <- tips$x + 0.0015 + nchar(tips$label)*CW + 0.0028     # primer punto = fin del nombre + hueco
tips$x0 <- pmin(tips$x0, host_x-0.004)
p <- p +
  geom_segment(data=tips, aes(x=x0, xend=host_x-0.004, y=y, yend=y),
               linetype="dotted", linewidth=0.22, color="grey60", inherit.aes=FALSE) +
  geom_text(data=tips, aes(x=host_x, y=y, label=info),
            inherit.aes=FALSE, hjust=0, size=host_sz, color="grey45", fontface="italic")

## ---- barras de linaje (mas a la izquierda, mas gruesas, titulos centrados) ----
host_end <- host_x + max(nchar(tips$info))*0.00140           # fin del texto de huesped
barx <- host_end + 0.006; LW <- 2.6
h <- 0.6*(length(st53)-1)/2
yI  <- range(d$y[d$label %in% linI])
yIB <- range(d$y[d$label %in% ib])
yII <- c(min(d$y[d$label %in% rel])-0.3, yI[1]-0.8)          # sube hasta casi tocar la verde
seg <- function(a,b,col) annotate("segment", x=barx, xend=barx, y=a, yend=b, color=col, linewidth=LW)
lab <- function(a,b,txt,col) annotate("text", x=barx+0.004, y=(a+b)/2, label=txt, color=col, angle=90, size=TS, fontface="bold")
p <- p + seg(yI[1],yI[2],green)  + lab(yI[1],yI[2],"Lineage I: 9a5c and Brazilian strains",green) +
         seg(yII[1],yII[2],blue) + lab(yII[1],yII[2],"Lineage II",blue) +
         seg(yIB[1],yIB[2],orange)+ lab(yIB[1],yIB[2],"Ibiza, ST80",orange)

xf<-d$x[d$label=="Fillmore"]; yf<-d$y[d$label=="Fillmore"]; xpar<-d$x[d$node==d$parent[d$label=="Fillmore"]]
p <- p + annotate("text", x=(xf+xpar)/2, y=yf, label="//", size=5, color="black")

## ---- inset: mapa de cafetales / registros de X. fastidiosa en cafe ----
for(pk in c("cowplot","png")) if(!requireNamespace(pk,quietly=TRUE)) install.packages(pk,repos="https://cloud.r-project.org")
suppressMessages(library(cowplot))
MAP <- "data/map_coffee_americas.png"   # <-- AJUSTA la ruta si esta en otro sitio
g <- grid::rasterGrob(png::readPNG(path.expand(MAP)), interpolate=TRUE)
# posicion/tamano del inset en coordenadas 0-1 del lienzo (dial): esquina inferior-izquierda
INS_X<-0.045; INS_Y<-0.05; INS_W<-0.36; INS_H<-0.36
final <- ggdraw(p) + draw_grob(g, x=INS_X, y=INS_Y, width=INS_W, height=INS_H)

ggsave("figures/Figure1.pdf", final, width=11.5, height=10, device=cairo_pdf)
ggsave("figures/Figure1.png", final, width=11.5, height=10, dpi=150, type="cairo")
cat("hecho\n")
