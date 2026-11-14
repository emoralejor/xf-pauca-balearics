suppressMessages({library(ggtree);library(treeio);library(ggplot2);library(ape)})
TREE <- "data/core78_ML.treefile"
tr <- read.tree(TREE)
labmap <- c("0297"="M5","0298"="M7")
for(i in seq_along(labmap)){ j<-which(tr$tip.label==names(labmap)[i]); if(length(j)) tr$tip.label[j]<-labmap[i] }
mallorca<-c("M5","M7"); menton<-c("MentOliv","MentPoly"); cr<-c("OLS0479","OLS0478"); cof<-"COF0407"
linB<-c("LSV47.09","LSV47.33","PD7202"); OUT<-"PD7211"
want<-c(mallorca,menton,cr,cof,linB,OUT)
miss<-setdiff(want,tr$tip.label); if(length(miss)){cat("FALTAN:\n");print(miss);cat("\nTIPS:\n");print(sort(tr$tip.label));quit(save="no")}
tr<-root(tr,outgroup=OUT,resolve.root=TRUE,edgelabel=TRUE)
apulia<-setdiff(tr$tip.label,want); cat("Apulia n =",length(apulia),"\n")

## --- colapsar Apulia: podar a 1 tip representativo (el triangulo se dibuja a mano) ---
rep_tip<-apulia[1]; tr2<-drop.tip(tr,setdiff(apulia,rep_tip))
APU<-"Apulia (67)"; tr2$tip.label[tr2$tip.label==rep_tip]<-APU

## --- acortar rama del outgroup respecto a la profundidad del ingroup ---
dall<-node.depth.edgelength(tr2); ing<-max(dall[which(tr2$tip.label!=OUT)])
fe<-which(tr2$edge[,2]==which(tr2$tip.label==OUT)); real_out<-tr2$edge.length[fe]
tr2$edge.length[fe]<-ing*0.5

## --- formatear soporte aLRT/UFboot (redondear aLRT) ---
fmt<-function(s){s<-as.character(s); if(is.na(s)||s=="")return(""); z<-strsplit(s,"/")[[1]]; if(length(z)<2)return(s); a<-suppressWarnings(as.numeric(z[1])); paste0(ifelse(is.na(a),z[1],round(a)),"/",z[2])}
tr2$node.label<-vapply(tr2$node.label,fmt,character(1))

orange<-"#E8820C"; blue<-"#1F78B4"; green<-"#1B9E77"; grey<-"grey45"
grp<-rep("other",length(tr2$tip.label)); names(grp)<-tr2$tip.label
grp[mallorca]<-"mallorca"; grp[c(cr,cof)]<-"cr"; grp[linB]<-"linB"; grp[OUT]<-"out"; grp[APU]<-"apu"
trg<-groupOTU(tr2,split(names(grp),grp))
rootn<-length(tr2$tip.label)+1

p<-ggtree(trg,aes(color=group),size=0.6)+
  scale_color_manual(values=c(other="grey30",mallorca=orange,cr=green,linB=grey,out=grey,apu=blue),guide="none")+
  geom_tippoint(aes(subset=(label!=APU),color=group),size=1.8)+
  geom_tiplab(aes(subset=(label!=APU)),size=3.0,offset=ing*0.015,color="black")+
  geom_text2(aes(subset=(!isTip & node!=rootn & label!="" & !is.na(label)),label=label),size=2.4,color="grey45",hjust=1.15,vjust=-0.45)+
  geom_treescale(width=1e-4,x=0,y=-0.4,fontsize=3.2)+
  theme(plot.margin=margin(5,55,5,5))+coord_cartesian(clip="off")

## --- triangulo Apulia dibujado a mano ---
d<-p$data; ar<-d[d$label==APU,]; xpar<-d$x[d$node==ar$parent]
W<-ing*0.7; hh<-0.85
tri<-data.frame(x=c(xpar,xpar+W,xpar+W),y=c(ar$y,ar$y-hh,ar$y+hh))
p<-p+geom_polygon(data=tri,aes(x=x,y=y),fill=blue,color=blue,alpha=0.55,inherit.aes=FALSE)+
    annotate("text",x=xpar+W+ing*0.04,y=ar$y,label=APU,color=blue,fontface="bold",size=3.1,hjust=0)

## --- corchetes de clado ---
node_ml<-getMRCA(tr2,mallorca); node_cr<-getMRCA(tr2,cr); node_B<-getMRCA(tr2,linB)
node_A<-getMRCA(tr2,c(APU,mallorca,menton,cr,cof)); apu_n<-which(tr2$tip.label==APU)
o1<-ing*0.4; o2<-ing*1.1
cl<-function(n,lab,off,bc,tc,fs,bs,ang,hj) geom_cladelab(node=n,label=lab,offset=off,barcolor=bc,textcolor=tc,fontsize=fs,barsize=bs,angle=ang,hjust=hj,align=TRUE,offset.text=ing*0.08)
p<-p+
  cl(apu_n,"",o1,blue,blue,3.0,2.2,0,0)+
  cl(node_ml,"Mallorca\nST53, 2025",o1,orange,orange,3.0,2.2,0,0)+
  cl(node_cr,"Costa Rica",o1,green,green,3.0,2.2,0,0)+
  cl(node_A,"Lineage A",o2,"grey20","black",3.5,2.8,90,0.5)+
  cl(node_B,"Lineage B",o2,"grey20","black",3.5,2.8,90,0.5)

ggsave("figures/FigureS2.pdf",p,width=7.6,height=6.5,device=cairo_pdf)
ggsave("figures/FigureS2.png",p,width=7.6,height=6.5,dpi=150,type="cairo")
cat("hecho\n")
