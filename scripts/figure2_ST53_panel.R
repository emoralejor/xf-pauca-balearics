suppressMessages({library(ape);library(ggplot2);library(grid)})
for(pk in c("pheatmap","patchwork","ggplotify","maps")) if(!requireNamespace(pk,quietly=TRUE)) install.packages(pk,repos="https://cloud.r-project.org")
library(pheatmap); library(patchwork); library(ggplotify)
ALN <- "data/core78.filtered_polymorphic_sites.fasta"
labmap <- c(
  "SRR30248731"="50080","SRR3134552"="CVC0251","SRR3134553"="CVC0256","SRR28211529"="50066",
  "SRR3134554"="COF0324","SRR3134557"="COF0407","SRR3134556"="OLS0478","SRR8501430"="CFBP8072",
  "SRR8187197"="APL42","SRR1707402"="CoDiRO","SRS22193438"="De_Donno","SRR32619863"="LSV47.09",
  "SRR32619866"="MentPoly","SRR3134555"="OLS0479","SRR8263259"="PD7202","SRR9849809"="Salento-1",
  "SRR9849810"="Salento-2","SRR8263205"="PD7211","SRR8187184"="APL12","SRR8187164"="APL18",
  "SRR8187140"="APL2","SRR8187151"="APL31","SRR8187152"="APL32","SRR8187194"="APL37",
  "SRR8187208"="APL50","SRR8187209"="APL71","SRR8187148"="APL9","SRR32619864"="LSV47.33",
  "SRR8187158"="APL57","SRR8187180"="APL64","SRR8187172"="APL69","SRR8187211"="APL72",
  "SRR8187168"="APL75","SRR8187169"="APL76","SRR32619865"="MentOliv","SRR22122615"="ALM4",
  "SRR22122611"="AM2-Angelina","SRR22122613"="OLI17A","SRR15620328"="Mag1_15620328","SRR10513216"="Mag1_10513216",
  "SRS22070466"="50088","SRS8523271"="XF71_BR","SRS8523272"="XF64.2","SRS5491332"="RAAR17_CiUb7",
  "SRS5491331"="RAAR16_Co13","SRS5491330"="RAAR15_Co33","SRR11931328"="XYL2017/18","SRR8082354"="XYL1732",
  "SRR11931329"="XYL2093/18","SRR27912192"="XYL1752/17","SRR11931333"="XYL1968/18","SRR10246938"="XF1094",
  "SRR10246946"="XF1110","SRR10246936"="XF70","SRR10246940"="XF1090","SRR10246941"="XF75",
  "SRR10246944"="XF72","SRR10246945"="XF71","SRR10246935"="XF1105","SRR10246937"="XF68",
  "SRR10246939"="XF1093","SRR10246943"="XF73","SRR10246942"="XF74","SRS8523265"="XF143",
  "SRS8523264"="XF217-2","SRS8523270"="XF221","SRS8523267"="XF43","SRS8523269"="XF45-A3",
  "SRS8523263"="XF-69-BR","SRS8523266"="XF-73-BR","SRS8523268"="XF39","SRR8187176"="APL30",
  "SRR8187141"="APL1","SRR8187147"="APL10","SRR8187185"="APL11","SRR8187187"="APL13",
  "SRR8187186"="APL14","SRR8187189"="APL15","SRR8187188"="APL17","SRR8187191"="APL19",
  "SRR8187190"="APL20","SRR8187193"="APL21","SRR8187165"="APL22","SRR8187192"="APL23",
  "SRR8187154"="APL24","SRR8187155"="APL25","SRR8187156"="APL26","SRR8187157"="APL27",
  "SRR8187150"="APL28","SRR8187139"="APL3","SRR8187153"="APL33","SRR8187162"="APL34",
  "SRR8187163"="APL35","SRR8187195"="APL36","SRR8187160"="APL38","SRR8187207"="APL39",
  "SRR8187138"="APL4","SRR8187199"="APL40","SRR8187198"="APL41","SRR8187196"="APL43",
  "SRR8187210"="APL46","SRR8187145"="APL5","SRR8187177"="APL54","SRR8187166"="APL55",
  "SRR8187161"="APL56","SRR8187144"="APL6","SRR8187159"="APL60","SRR8187146"="APL61",
  "SRR8187149"="APL62","SRR8187200"="APL63","SRR8187179"="APL65","SRR8187178"="APL66",
  "SRR8187206"="APL67","SRR8187183"="APL68","SRR8187143"="APL7","SRR8187201"="APL70",
  "SRR8187182"="APL73","SRR8187167"="APL74","SRR8187142"="APL8","SRR11638835"="Fillmore",
  "Reference"="De_Donno","0297_M5_S15_L001"="M5","0298_M7_S16_L001"="M7")

ren <- function(v){i<-match(v,names(labmap)); v[!is.na(i)]<-labmap[i[!is.na(i)]]; v}
grp <- function(x){
  if(x %in% c("M5","M7")) return("Mallorca")
  if(grepl("^APL|De.?Donno|CoDiRO|Salento", x, ignore.case=TRUE)) return("Apulia")
  if(grepl("^Ment", x)) return("France")
  if(grepl("^OLS|^COF", x)) return("Costa Rica")
  if(grepl("^LSV|PD7202", x)) return("Lineage B")
  "other"}
gcol <- c(Apulia="#1F78B4", Mallorca="#D62728", France="#2CA02C", "Costa Rica"="#7E3F9E", "Lineage B"="grey45")
m <- as.character(as.matrix(read.dna(ALN, format="fasta"))); rownames(m) <- ren(rownames(m))
m <- m[rownames(m)!="PD7211",,drop=FALSE]; labs<-rownames(m); n<-length(labs)
valid<-function(v) v %in% c("a","c","g","t"); D<-matrix(0,n,n,dimnames=list(labs,labs))
for(i in 1:(n-1)) for(j in (i+1):n){a<-m[i,];b<-m[j,];ok<-valid(a)&valid(b);dd<-sum(ok&a!=b);D[i,j]<-dd;D[j,i]<-dd}
ann <- data.frame(Lineage=vapply(labs,grp,character(1)), row.names=labs); annc <- list(Lineage=gcol)
brk <- unique(seq(0, sqrt(max(D)), length.out=101)^2)
pal <- colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(length(brk)-1)
d <- as.dist(D)
ph <- pheatmap(D, clustering_distance_rows=d, clustering_distance_cols=d, clustering_method="average",
               color=pal, breaks=brk, annotation_row=ann, annotation_col=ann, annotation_colors=annc,
               fontsize_row=5, fontsize_col=5, treeheight_row=45, treeheight_col=45, border_color=NA, silent=TRUE,
               main="Pairwise non-recombinant SNP differences (77 ST53 isolates)")
A <- as.ggplot(ph$gtable)
anc<-c("GCC","AAC","AGC","ACC","ATG"); aaA<-c("Ala","Asn","Ser","Thr","Met")
der<-c("GCC","AAC","GGC","ACC","ATG"); aaD<-c("Ala","Asn","Gly","Thr","Met")
dB<-rbind(data.frame(x=1:5,y=2,cod=anc,aa=aaA,hl=c(F,F,T,F,F)),data.frame(x=1:5,y=1,cod=der,aa=aaD,hl=c(F,F,T,F,F)))
B <- ggplot(dB,aes(x,y))+geom_tile(aes(fill=hl),colour="grey40",width=0.92,height=0.6)+
  geom_text(aes(label=cod),size=3,fontface="bold")+geom_text(aes(y=y-0.42,label=aa),size=2.4,colour="grey30")+
  annotate("text",x=3,y=1.5,label="A → G",size=2.8,fontface="bold",colour="#B2182B")+
  annotate("text",x=0.3,y=2,label="Ancestral",hjust=1,size=2.6)+annotate("text",x=0.3,y=1,label="Derived",hjust=1,size=2.6)+
  scale_fill_manual(values=c(`FALSE`="grey95",`TRUE`="#FDDBC7"),guide="none")+
  coord_cartesian(xlim=c(-1.1,5.6),ylim=c(0.3,2.7),clip="off")+
  labs(title="ALS2, codon 565 (pos 1,251,537): Ser → Gly",
       subtitle="ancestral (Ser): Mallorca, Costa Rica, European ST53  —  derived (Gly): Apulia only")+
  theme_void(base_size=8)+theme(plot.title=element_text(size=8.5,face="bold"),plot.subtitle=element_text(size=6,colour="grey30"),plot.margin=margin(2,2,2,2))
dC<-data.frame(x=c(1,2),lab=c("All other ST53\n(incl. 67 Apulian)","M5, M7"),st=c("present","absent"))
C <- ggplot(dC,aes(x,1))+geom_tile(aes(fill=st),colour="grey40",width=0.8,height=0.5)+
  geom_text(aes(label=lab),size=3.4)+scale_fill_manual(values=c(present="#92C5DE",absent="white"),name=NULL)+
  coord_cartesian(ylim=c(0.55,1.45))+labs(title="Plasmid pXF-De_Donno")+
  theme_void(base_size=8)+theme(plot.title=element_text(size=9.5,face="bold"),legend.position="bottom",
        legend.text=element_text(size=8),legend.key.size=unit(0.35,"cm"),plot.margin=margin(2,2,2,2))
## D: mapa — introduccion transatlantica a un reservorio europeo difuso, luego difusion a los focos
wmap <- map_data("world")
pts <- data.frame(site=c("Costa Rica","Apulia","Menton","Mallorca"),
                  lon=c(-84.1,17.9,7.5,2.9), lat=c(9.9,40.4,43.8,39.6),
                  g=c("Costa Rica","Apulia","France","Mallorca"),
                  nx=c(0,1.0,-0.5,-2.8), ny=c(-4,-2.8,-2.8,-3.2))
bl <- c(5,52)   # blob difuso ~ Paises Bajos
arrT <- data.frame(x=-84.1,y=9.9,xend=bl[1],yend=bl[2])
arrE <- data.frame(x=bl[1],y=bl[2],xend=pts$lon[2:4],yend=pts$lat[2:4])
blob <- data.frame(s=c(34,27,20,14,9), a=c(0.05,0.06,0.08,0.10,0.12))
Dmap <- ggplot()+
  geom_polygon(data=wmap, aes(long,lat,group=group), fill="grey92", color="grey75", linewidth=0.12)+
  lapply(seq_len(nrow(blob)), function(i) annotate("point", x=bl[1], y=bl[2], size=blob$s[i], color="#4292C6", alpha=blob$a[i]))+
  geom_curve(data=arrT, aes(x=x,y=y,xend=xend,yend=yend), curvature=0.14, linewidth=0.55, color="grey35",
             arrow=arrow(length=unit(0.12,"cm"),type="closed"))+
  geom_curve(data=arrE, aes(x=x,y=y,xend=xend,yend=yend), curvature=-0.25, linewidth=0.32, color="grey45",
             linetype="22", arrow=arrow(length=unit(0.08,"cm"),type="closed"), alpha=0.9)+
  geom_point(data=pts, aes(lon,lat,color=g), size=2.2)+
  annotate("text", x=bl[1], y=bl[2]+6, label="putative European reservoir\n(uncertain entry / intra-EU trade)",
           size=1.85, fontface="italic", color="grey40", lineheight=0.9)+
  geom_text(data=pts, aes(lon+nx,lat+ny,label=site,color=g), size=2.2, fontface="bold", show.legend=FALSE)+
  scale_color_manual(values=gcol, guide="none")+
  coord_quickmap(xlim=c(-95,22), ylim=c(3,60))+
  labs(title="ST53: Central American origin and route into Europe",
       subtitle="solid: single transatlantic introduction  —  dashed: circulation / establishment within Europe")+
  theme_void(base_size=8)+theme(plot.title=element_text(size=8.5,face="bold"),
        plot.subtitle=element_text(size=5.3,colour="grey35"),
        panel.border=element_rect(fill=NA,color="grey70",linewidth=0.3), plot.margin=margin(2,16,2,2))

left <- (B / C) + plot_layout(heights=c(1.35,0.65))
bottom <- (left | Dmap) + plot_layout(widths=c(1.5,0.82))
fig <- A / bottom + plot_layout(heights=c(3,1.25)) + plot_annotation(tag_levels="A")
ggsave("figures/Figure2.pdf", fig, width=12, height=11, device=cairo_pdf)
ggsave("figures/Figure2.png", fig, width=12, height=11, dpi=150, type="cairo")
cat("hecho; n =",n,"\n")
