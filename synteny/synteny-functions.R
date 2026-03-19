
`%notin%` <- Negate(`%in%`)


plotSynt <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(colorRampPalette(brewer.pal(6,'BuPu'))(nrow(s2chLim)))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  plab=paste("(",nrow(synt)," orthologues)",sep="")
  labxp=paste(labx,plab,sep=" ")
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col="#6A3D9A",xlab=labxp,ylab=laby)
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',as.factor(synt$s1chr),'lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S','#810f7c','lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.6)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.6)
  
}

plotSynts <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(colorRampPalette(brewer.pal(6,'BuPu'))(nrow(s2chLim)))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  plab=paste("(",nrow(synt)," orthologues)",sep="")
  labxp=paste(labx,plab,sep=" ")
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s2chr,xlab=labxp,ylab=laby)
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S','#6A3D9A','lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',as.factor(synt$s1chr),'lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S','#810f7c','lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.6)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.6)
  
}

plotSyntCLGr <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  plab=paste("(",nrow(synt)," orthologues)",sep="")
  labxp=paste(labx,plab,sep=" ")
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=synt$clg,xlab=labxp,ylab=laby)
  #plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labx,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(h=s1chLim$max,lty=3,lwd=0.5)
  abline(v=s2chLim$max,lty=3,lwd=0.5)
  
  axis(4,at=s1chLim$mids,labels=s1chLim$s1chr,las=1,cex.axis=0.6)
  axis(3,at=s2chLim$mids,labels=s2chLim$s2chr,las=2,cex.axis=0.6)
}
plotSyntCLGrs <- function(synt,labx,laby){
  synt  %>% dplyr::select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  plab=paste("(",nrow(synt)," orthologues)",sep="")
  labxp=paste(labx,plab,sep=" ")
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  #plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labxp,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(h=s1chLim$max,lty=3,lwd=0.5)
  abline(v=s2chLim$max,lty=3,lwd=0.5)
  
  axis(4,at=s1chLim$mids,labels=s1chLim$s1chr,las=1,cex.axis=0.6)
  axis(3,at=s2chLim$mids,labels=s2chLim$s2chr,las=2,cex.axis=0.6)
}

plotSyntCLGrsi <- function(synt,labx,laby){
  synt  %>% dplyr::select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  plab=paste("(",nrow(synt)," orthologues)",sep="")
  labxp=paste(labx,plab,sep=" ")
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  #plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labxp,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.6)
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.6)
}

plotSyntCLGoc <- function(synt,labx,laby,c){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labx,ylab=laby)
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.6)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.6)
  
}


ordSynt <- function(synt){
  synt %>% dplyr::count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  return(synt.ord)
  #return(c(x_chrom,y_chrom))
}

ordSyntRw <- function(synt){
  synt %>% count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  #y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% 
    #mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    #arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  return(synt.ord)
  #return(c(x_chrom,y_chrom))
}

ordSyntCl <- function(synt){
  synt %>% count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  #x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% 
    mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) -> synt.ord
    #mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    #arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  return(synt.ord)
  #return(c(x_chrom,y_chrom))
}



ordSyntJO <- function(synt){
  synt %>% count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  #return(synt.ord)
  return(list(x_chrom,y_chrom))
}

testEnrichCLG <- function(synt){
  synt %>% 
    group_by(clg) %>% mutate(clgtot=n()) %>%
    group_by(s2chr) %>% mutate(chrtot=n())%>%
    group_by(s2chr,clgtot,chrtot,clg) %>% tally() %>% arrange(clg) -> synt.chrn
  synt.no=sum(synt.chrn$n)
  synt.chrn %>% rowwise() %>% 
    mutate(fishp=fisher.test(matrix(c(n,clgtot-n,chrtot-n,synt.no),ncol=2),alternative='greater')$p.value) %>%
    ungroup() %>% mutate(padj=p.adjust(fishp,method="bonferroni")) %>% mutate(scol=ifelse(padj<.05,'S','NS')) -> synt.test
  synt.test %>% inner_join(synt) -> synt.exp
  return(synt.test)
}

testEnrichCLG1 <- function(synt){
  synt %>% 
    group_by(clg) %>% mutate(clgtot=n()) %>%
    group_by(s1chr) %>% mutate(chrtot=n())%>%
    group_by(s1chr,clgtot,chrtot,clg) %>% tally() %>% arrange(clg) -> synt.chrn
  synt.no=sum(synt.chrn$n)
  synt.chrn %>% rowwise() %>% 
    mutate(fishp=fisher.test(matrix(c(n,clgtot-n,chrtot-n,synt.no),ncol=2),alternative='greater')$p.value) %>%
    ungroup() %>% mutate(padj=p.adjust(fishp,method="bonferroni")) %>% mutate(scol=ifelse(padj<.05,'S','NS')) -> synt.test
  #synt.test %>% inner_join(synt) -> synt.exp
  return(synt.test)
}


testEnrich <- function(synt){
  synt %>% 
    group_by(s1chr) %>% mutate(s1tot=n()) %>%
    group_by(s2chr) %>% mutate(s2tot=n())%>%
    group_by(s1chr,s2chr,s1tot,s2tot) %>% tally() -> synt.chrn
  synt.no=sum(synt.chrn$n)
  synt.chrn %>% rowwise() %>% 
    mutate(fishp=fisher.test(matrix(c(n,s1tot-n,s2tot-n,synt.no),ncol=2),alternative='greater')$p.value) %>%
    ungroup() %>% mutate(padj=p.adjust(fishp,method="bonferroni")) %>% mutate(scol=ifelse(padj<.05,'S','NS')) -> synt.test
  synt.test %>% inner_join(synt) -> synt.exp
  return(synt.exp)
} 
gchkExp <- function(synt.exp){
  synt.exp %>% arrange(s1chr,s1gi) %>% #filter(s1chr=="Scaffold_2100") %>% 
    group_by(s1chr) %>% mutate(chunk=cut_width(s1gi,width=20)) %>% ungroup() %>% filter(scol=="S") %>% 
    select(s2chr,chunk,s1chr,s1gp) %>% 
    group_by(s1chr,chunk) %>% mutate(xmin=min(s1gp),xmax=max(s1gp)) %>%
    group_by(s2chr,s1chr,chunk,xmin,xmax) %>% tally() %>% 
    group_by(chunk) %>% mutate(ymin=(cumsum(n)-n)/sum(n),ymax=cumsum(n)/(sum(n))) -> synt.gchk
  return(synt.gchk)
}

gchkExpCLG <- function(synt.exp){
  synt.exp %>% arrange(s2chr,s2gi) %>% #filter(s1chr=="Scaffold_2100") %>% 
    group_by(s2chr) %>% mutate(chunk=cut_width(s2gi,width=20)) %>% ungroup() %>% filter(scol=="S") %>% 
    select(clg,chunk,s2chr,s2gp) %>% 
    group_by(s2chr,chunk) %>% mutate(xmin=min(s2gp),xmax=max(s2gp)) %>%
    group_by(clg,s2chr,chunk,xmin,xmax) %>% tally() %>% 
    group_by(chunk) %>% mutate(ymin=(cumsum(n)-n)/sum(n),ymax=cumsum(n)/(sum(n))) -> synt.gchk
  synt.gchk$s2chr=fct_relevel(synt.gchk$s2chr,str_sort(levels(synt.gchk$s2chr), numeric = TRUE))
  return(synt.gchk)
}
