library(ggtree); library(treeio); library(ggplot2)

tr <- read.beast("data/apulia_2013_mcc.tree")
ph <- as.phylo(tr)
apu <- c("APL2","APL9","APL12","APL18","APL31","APL32","APL37","APL42","APL50","APL57",
         "APL64","APL69","APL71","APL72","APL75","APL76","CoDiRO","De_Donno","Salento-1","Salento-2")
mal <- c("0297_M5","0298_M7"); fra <- c("MentOliv","MentPoly"); cr <- c("COF0407","OLS0478","OLS0479")
c_mal<-"#c0392b"; c_fra<-"#2e7d32"; c_cr<-"#8e44ad"; italia<-"#0072BC"

node_apu <- ape::getMRCA(ph, apu)
d0 <- ggtree(tr, mrsd="2025-06-30")$data
apu_x_max <- max(d0$x[d0$label %in% apu], na.rm=TRUE)
x_am  <- d0$x[d0$node==ape::getMRCA(ph, c(apu, mal))]
pn <- suppressWarnings(as.numeric(d0$posterior))
x_066 <- d0$x[d0$node==d0$node[!d0$isTip & !is.na(pn) & abs(pn-0.66)<0.02][1]]

TS <- 3.4  # tamaño base de texto (≈30% mayor que antes)

p <- ggtree(tr, mrsd="2025-06-30", size=0.4) +
  theme_tree2() +
  geom_range(range="height_0.95_HPD", center="height", color="grey65", alpha=0.6, size=2) +
  geom_text2(aes(subset=(!isTip & !is.na(as.numeric(posterior)) &
                         as.numeric(posterior)>=0.5 & as.numeric(posterior)<0.99),
                 x=branch, label=sprintf("%.2f", as.numeric(posterior))),
             size=3.0, vjust=-0.4, color="grey30") +
  geom_text2(aes(subset=(!isTip & !is.na(as.numeric(posterior)) & as.numeric(posterior)>=0.99),
                 label="\u2739"), size=4.2, color="black", vjust=0.35) +
  # puntas coloreadas por grupo (capas de color fijo, robustas en PDF)
  geom_tippoint(aes(subset=(label %in% mal)), color=c_mal, size=2.4) +
  geom_tippoint(aes(subset=(label %in% fra)), color=c_fra, size=2.4) +
  geom_tippoint(aes(subset=(label %in% cr )), color=c_cr,  size=2.4) +
  geom_tiplab(aes(subset=(label %in% mal)), color=c_mal, size=TS, offset=0.4) +
  geom_tiplab(aes(subset=(label %in% fra)), color=c_fra, size=TS, offset=0.4) +
  geom_tiplab(aes(subset=(label %in% cr )), color=c_cr,  size=TS, offset=0.4) +
  geom_tiplab(aes(subset=(!label %in% c(mal,fra,cr,apu))), color="black", size=TS, offset=0.4) +
  geom_strip("0298_M7","0297_M5", label="Mallorca", align=TRUE,
             offset=4, offset.text=0.8, barsize=0.7, fontsize=TS, color=c_mal) +
  scale_x_continuous(breaks=seq(1900,2020,20), expand=expansion(mult=c(0.02,0.0))) +
  coord_cartesian(xlim=c(1885,2026), clip="off") +
  theme(plot.margin=margin(5,95,5,5),
        axis.text.x=element_text(size=12, color="black", face="bold"))

p <- scaleClade(p, node=node_apu, scale=0.25)
p <- collapse(p, node=node_apu, mode="mixed", fill=italia, alpha=0.6, color=italia)

yr <- range(p$data$y, na.rm=TRUE)
p <- p +
  annotate("segment", x=x_am,  xend=x_am,  y=yr[1]-0.4, yend=yr[2]+0.4,
           linetype="dashed", color="grey55", linewidth=0.4) +
  annotate("segment", x=x_066, xend=x_066, y=yr[1]-0.4, yend=yr[2]+0.4,
           linetype="dashed", color="grey55", linewidth=0.4) +
  annotate("text", x=x_am,  y=yr[2]+1.0, label="2007", color="grey35", size=TS) +
  annotate("text", x=x_066, y=yr[2]+1.0, label="1990", color="grey35", size=TS)

y_apu <- p$data$y[p$data$node==node_apu]
p <- p + annotate("text", x=apu_x_max+1, y=y_apu, label="Apulia (n=20)",
                  color=italia, size=TS, hjust=0, vjust=0.5)

# leyenda manual (color fijo, en zona vacia arriba-izquierda)
p <- p +
  annotate("text", x=1888, y=yr[2]-0.2, label="Mallorca",   color=c_mal, hjust=0, size=TS, fontface="bold") +
  annotate("text", x=1888, y=yr[2]-1.2, label="France",     color=c_fra, hjust=0, size=TS, fontface="bold") +
  annotate("text", x=1888, y=yr[2]-2.2, label="Costa Rica", color=c_cr,  hjust=0, size=TS, fontface="bold")

ggsave("figures/Figure4.pdf", p, width=9.5, height=6.8, device=cairo_pdf)
ggsave("figures/Figure4.png", p, width=9.5, height=6.8, dpi=150, type="cairo")
cat("hecho\n")
