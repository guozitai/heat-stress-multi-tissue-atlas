.libPaths(c("/home/guozitai/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
suppressMessages({library(dplyr); library(readxl); library(stringr)})
outdir <- "."
filter <- dplyr::filter
dat <- read_excel("/home/guozitai/20250220/targetgene_enrich.xlsx") %>%
  setNames(c("mirna","gene","regulation","tissue")) %>%
  mutate(tissue = case_when(str_detect(tissue, regex("fat|adipose", TRUE)) ~ "Adipose",
                            str_detect(tissue, regex("muscle", TRUE))      ~ "Muscle"),
         regulation = case_when(str_detect(regulation, regex("up", TRUE))   ~ "Up",
                                str_detect(regulation, regex("down", TRUE)) ~ "Down")) %>%
  filter(!is.na(tissue))
TOP_N <- 4; R_MIR <- 2.5; R_GENE <- 5.9; CX <- 7.0
make_polar <- function(data, cx, cy) {
  ## 源表每个 miRNA-基因对按其所属富集通路重复出现,先去重再取 top N,
  ## 否则同一个基因会被当成多个不同靶点画成多个节点
  dat_t  <- data %>% distinct(mirna, gene, .keep_all = TRUE) %>%
            group_by(mirna) %>% slice_head(n = TOP_N) %>% ungroup()
  mirnas <- unique(dat_t$mirna); n_mir <- length(mirnas)
  mir_ang <- seq(pi/2, pi/2 + 2*pi, length.out = n_mir + 1)[1:n_mir]
  mn<-list(); gn<-list(); ed<-list()
  for (i in seq_along(mirnas)) {
    m<-mirnas[i]; ang<-mir_ang[i]
    mx<-cx+R_MIR*cos(ang); my<-cy+R_MIR*sin(ang)
    g<-dat_t %>% filter(mirna==m); n_g<-nrow(g); sw<-2*pi/n_mir*0.76
    ga<-if(n_g==1) ang else seq(ang-sw/2, ang+sw/2, length.out=n_g)
    mn[[i]]<-data.frame(name=m,x=mx,y=my,ang=ang,size=12+n_g*2,stringsAsFactors=FALSE)
    for (j in seq_len(n_g)) {
      gx<-cx+R_GENE*cos(ga[j]); gy<-cy+R_GENE*sin(ga[j])
      gn[[length(gn)+1]]<-data.frame(name=g$gene[j],x=gx,y=gy,ang=ga[j],
        regulation=g$regulation[j],stringsAsFactors=FALSE)
      ed[[length(ed)+1]]<-data.frame(x0=mx,y0=my,x1=gx,y1=gy,
        regulation=g$regulation[j],stringsAsFactors=FALSE)
    }
  }
  list(mir=bind_rows(mn), gene=bind_rows(gn), edge=bind_rows(ed))
}
col_adi_mir<-"#FC8D62"; col_adi_gene<-"#FDE0D0"; col_adi_brd<-"#E07040"
col_mus_mir<-"#66C2A5"; col_mus_gene<-"#D4EDE4"; col_mus_brd<-"#4AA88A"
col_up<-"#C0392B"; col_dn<-"#2471A3"
ring<-function(cx,cy,r,col){t<-seq(0,2*pi,length.out=360);lines(cx+r*cos(t),cy+r*sin(t),col=col,lty=3,lwd=0.4)}
edges<-function(e){for(i in seq_len(nrow(e))){r<-e[i,];cc<-if(r$regulation=="Up") col_up else col_dn
  arrows(r$x0,r$y0,r$x1,r$y1,length=0.045,angle=18,code=2,col=adjustcolor(cc,0.7),lwd=0.7)}}
halo<-function(x,y,lab,cex,adj,col="black",font=1,hw=0.016){
  for(dx in c(-hw,0,hw)) for(dy in c(-hw,0,hw)) if(dx!=0||dy!=0)
    text(x+dx,y+dy,lab,cex=cex,col="white",font=font,adj=adj,xpd=NA)
  text(x,y,lab,cex=cex,col=col,font=font,adj=adj,xpd=NA)}

## ---------- 只把真正相撞的标签沿 y 推开 ----------
resolve_y <- function(b, pad=0.055, iters=600){
  n<-nrow(b)
  for(it in seq_len(iters)){
    moved<-FALSE
    for(i in seq_len(n-1)) for(j in (i+1):n){
      if(b$x1[i] > b$x0[j] && b$x0[i] < b$x1[j]){          # x 方向有交叠才可能撞
        dy <- b$y[i]-b$y[j]
        need <- (b$h[i]+b$h[j])/2 + pad
        if(abs(dy) < need){
          s  <- if(dy > 0) 1 else if(dy < 0) -1 else (if(i%%2==0) 1 else -1)
          ov <- need-abs(dy)
          b$y[i] <- b$y[i]+s*ov/2; b$y[j] <- b$y[j]-s*ov/2
          moved<-TRUE
        }
      }
    }
    if(!moved){cat("  收敛于第",it,"轮\n"); break}
  }
  b
}
build_box <- function(lab,x,y,adjx,cex,font,kind,nx,ny){
  w<-strwidth(lab,cex=cex,font=font); h<-strheight(lab,cex=cex,font=font)
  data.frame(lab=lab,x=x,y=y,y0=y,adjx=adjx,cex=cex,font=font,kind=kind,
             nx=nx,ny=ny,w=w,h=h,
             x0=ifelse(adjx==0,x,ifelse(adjx==1,x-w,x-w/2)),
             x1=ifelse(adjx==0,x+w,ifelse(adjx==1,x,x+w/2)),
             stringsAsFactors=FALSE)
}
draw <- function(){
  adi<-make_polar(filter(dat,tissue=="Adipose"),-CX,0)
  mus<-make_polar(filter(dat,tissue=="Muscle"),  CX,0)
  par(mar=c(1.2,0.8,2.0,0.8), bg="white")
  plot(NULL,xlim=c(-14.9,14.9),ylim=c(-8.7,8.7),asp=1,axes=FALSE,xlab="",ylab="")
  for(cx in c(-CX,CX)){cc<-if(cx<0) col_adi_mir else col_mus_mir
    ring(cx,0,R_MIR,adjustcolor(cc,0.25)); ring(cx,0,R_GENE,adjustcolor(cc,0.12))}
  abline(v=0,col="grey82",lty=2,lwd=0.7)
  edges(adi$edge); edges(mus$edge)
  points(adi$gene$x,adi$gene$y,pch=21,cex=0.85,bg=col_adi_gene,col=col_adi_brd,lwd=0.6)
  points(mus$gene$x,mus$gene$y,pch=21,cex=0.85,bg=col_mus_gene,col=col_mus_brd,lwd=0.6)
  points(adi$mir$x,adi$mir$y,pch=21,cex=adi$mir$size/11,bg=col_adi_mir,col="white",lwd=1.5)
  points(mus$mir$x,mus$mir$y,pch=21,cex=mus$mir$size/11,bg=col_mus_mir,col="white",lwd=1.5)

  boxes<-list()
  for(p in list(list(g=adi$gene,cx=-CX),list(g=mus$gene,cx=CX))){
    ## 去重已在 make_polar 里按 (miRNA, 基因) 对完成,所以每个节点都是一个不同的
    ## 预测靶点 -> 每个节点都标注。被两条 miRNA 共同靶向的基因会出现两次,这是实情。
    g<-p$g
    for(i in seq_len(nrow(g))){
      a<-g$ang[i]; tx<-p$cx+(R_GENE+0.34)*cos(a); ty<-(R_GENE+0.34)*sin(a)
      boxes[[length(boxes)+1]]<-build_box(g$name[i],tx,ty,if(cos(a)>=0) 0 else 1,
        0.40,1,"gene",g$x[i],g$y[i])
    }
  }
  for(p in list(list(m=adi$mir,cx=-CX),list(m=mus$mir,cx=CX))){
    m<-p$m; lbl<-sub("bta-","",m$name)
    for(i in seq_len(nrow(m))){
      a<-m$ang[i]
      if(abs(cos(a))<0.30){tx<-m$x[i]; ty<-m$y[i]+sign(sin(a))*0.42; adjx<-0.5}
      else {tx<-p$cx+(R_MIR+0.30)*cos(a); ty<-(R_MIR+0.30)*sin(a); adjx<-if(cos(a)>=0) 0 else 1}
      boxes[[length(boxes)+1]]<-build_box(lbl[i],tx,ty,adjx,0.50,2,"mir",m$x[i],m$y[i])
    }
  }
  b<-bind_rows(boxes)
  cat("标签总数",nrow(b),"\n")
  b2<-resolve_y(b)
  bad<-0
  for(i in seq_len(nrow(b2)-1)) for(j in (i+1):nrow(b2)){
    if(b2$x1[i]>b2$x0[j] && b2$x0[i]<b2$x1[j] &&
       abs(b2$y[i]-b2$y[j]) < (b2$h[i]+b2$h[j])/2){
      bad<-bad+1; cat("  !! 仍重叠:",b2$lab[i],"x",b2$lab[j],"
")}}
  cat("残留重叠对数 =",bad,"
")
  cat("最大位移 =",round(max(abs(b2$y-b2$y0)),3),"
")
  shifted<-which(abs(b2$y-b2$y0)>1e-6)
  cat("被推开的标签:",length(shifted),"->",paste(b2$lab[shifted],collapse=", "),"\n")
  for(i in seq_len(nrow(b2))){
    r<-b2[i,]
    if(abs(r$y-r$y0)>0.14){          # 挪得明显的才补一根细引出线
      ax<-if(r$adjx==0) r$x0-0.06 else if(r$adjx==1) r$x1+0.06 else r$x
      segments(r$nx,r$ny,ax,r$y,col="grey72",lwd=0.35,xpd=NA)
    }
    if(r$kind=="gene") text(r$x,r$y,r$lab,cex=r$cex,col="grey20",adj=c(r$adjx,0.5),xpd=NA)
    else halo(r$x,r$y,r$lab,cex=r$cex,adj=c(r$adjx,0.5),font=2)
  }
  text(-CX,8.2,"Adipose",cex=1.05,font=2,col=col_adi_brd,adj=c(0.5,0.5))
  text( CX,8.2,"Muscle", cex=1.05,font=2,col=col_mus_mir,adj=c(0.5,0.5))
  legend("bottomright",inset=c(0.002,0.002),legend=c("Up-regulated","Down-regulated"),
    col=c(col_up,col_dn),lty=1,lwd=1.5,pch=NA,cex=0.55,bty="n",title="Regulation",title.col="grey20",title.cex=0.58)
  legend("bottomleft",inset=c(0.002,0.002),legend=c("Adipose miRNA","Muscle miRNA"),
    pt.bg=c(col_adi_mir,col_mus_mir),pch=21,pt.cex=1.1,col="white",pt.lwd=0.8,cex=0.55,bty="n",
    title="Tissue",title.col="grey20",title.cex=0.58)
}
cairo_pdf(file.path(outdir,"Fig4b_network.pdf"), width=9.8, height=5.9)
draw(); dev.off()
png(file.path(outdir,"Fig4b_network.png"), width=9.8, height=5.9, units="in", res=600, type="cairo")
draw(); dev.off()
cat("done
")
