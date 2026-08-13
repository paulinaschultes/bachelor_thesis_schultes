# =============================================================================
# Venn diagrams of shared top regions, labelled with Allen CCF acronyms
#
#   Figure 26 - top 40 high-intensity, low-variability regions per sex
#   Figure 29 - top 20 candidate input regions per retrograde cohort
#
# Regions are written into the diagram rather than counted, so the figure
# doubles as the region list. White background, no additional packages
# required beyond dplyr and ggplot2.
#
# Outputs (written to ~/Desktop/Bachelor):
#   venn2_sex.png, venn2_retro.png
# =============================================================================

suppressMessages({library(dplyr); library(ggplot2)})
B <- "~/Desktop/Bachelor"
pink <- "#D6547F"; violet <- "#6B4C9A"
topn <- function(d, rank, acr, n) d %>% arrange(.data[[rank]]) %>% slice_head(n=n) %>% pull(.data[[acr]])

sx <- read.csv(file.path(B,"sex_top_regions_low_deviation.csv"), stringsAsFactors=FALSE)
Fm <- topn(filter(sx, sex=="Female"),"combined","region_acronym",40)
Ml <- topn(filter(sx, sex=="Male"),  "combined","region_acronym",40)
r1 <- read.csv(file.path(B,"retrograde","input_regions_RET1_restraint_female_n3.csv"), stringsAsFactors=FALSE)
r3 <- read.csv(file.path(B,"retrograde","input_regions_RET3_forcedswim_male_n2.csv"), stringsAsFactors=FALSE)
R1 <- topn(r1,"combined_rank","acronym",20); R3 <- topn(r3,"combined_rank","acronym",20)

circle <- function(x0,r,n=300){t<-seq(0,2*pi,length.out=n); data.frame(x=x0+r*cos(t),y=r*sin(t))}

# lay labels out in ncol columns, centred on (cx,0)
grid_lab <- function(v, cx, ncol, dx, dy){
  v <- sort(v); n <- length(v); if(!n) return(NULL)
  nrow <- ceiling(n/ncol)
  i <- seq_len(n)-1
  col <- i %/% nrow; row <- i %% nrow
  data.frame(lab=v,
             x = cx + (col - (ncol-1)/2)*dx,
             y = ((nrow-1)/2 - row)*dy)
}

venn <- function(a, b, na, nb, title, ncol_out, dx_out, dy, size_out, size_sh, sep=0.74, r=1.0){
  onlyA <- setdiff(a,b); onlyB <- setdiff(b,a); sh <- intersect(a,b)
  cA <- circle(-sep,r); cB <- circle(sep,r)
  lA <- grid_lab(onlyA, -sep-0.26, ncol_out, dx_out, dy)
  lB <- grid_lab(onlyB,  sep+0.26, ncol_out, dx_out, dy)
  lS <- grid_lab(sh, 0, 1, 0, dy)
  ggplot()+
    geom_polygon(data=cA,aes(x,y),fill=pink,alpha=.35,colour="grey30",linewidth=.5)+
    geom_polygon(data=cB,aes(x,y),fill=violet,alpha=.35,colour="grey30",linewidth=.5)+
    geom_text(data=lA,aes(x,y,label=lab),size=size_out,colour="grey15")+
    geom_text(data=lB,aes(x,y,label=lab),size=size_out,colour="grey15")+
    geom_text(data=lS,aes(x,y,label=lab),size=size_sh,fontface="bold",colour="grey5")+
    annotate("text",x=-sep-0.26,y=r+0.20,label=sprintf("%s (%d)",na,length(a)),size=4.2,fontface="bold",colour=pink)+
    annotate("text",x= sep+0.26,y=r+0.20,label=sprintf("%s (%d)",nb,length(b)),size=4.2,fontface="bold",colour=violet)+
    annotate("text",x=0,y=-r-0.22,label=sprintf("%d shared",length(sh)),size=4,fontface="bold",colour="grey20")+
    ggtitle(title)+coord_fixed(xlim=c(-2.3,2.3),ylim=c(-1.42,1.42))+
    theme_void()+
    theme(plot.background=element_rect(fill="white",colour=NA),
          panel.background=element_rect(fill="white",colour=NA),
          plot.title=element_text(hjust=.5,size=13,face="bold",margin=margin(b=6)),
          plot.margin=margin(8,8,8,8))
}

p1 <- venn(Fm,Ml,"Female","Male","Top 40 low-deviation regions per sex",
           ncol_out=2, dx_out=0.40, dy=0.107, size_out=2.25, size_sh=2.8)
p2 <- venn(R1,R3,"Restraint (female)","Forced swim (male)","Top 20 candidate input regions",
           ncol_out=2, dx_out=0.44, dy=0.185, size_out=2.85, size_sh=3.3)
ggsave(file.path(B,"venn2_sex.png"),p1,width=8,height=4.6,dpi=300,bg="white")
ggsave(file.path(B,"venn2_retro.png"),p2,width=8,height=4.6,dpi=300,bg="white")
cat("ok\n")
