#### heatmap, box, bar, color_bar: spectral, green2orange, blue2red
draw.gsBygroup <- function(use_gs2gene=NULL,use_exp_mat=NULL,use_p=NULL,use_z_mat=NULL,use_gene=NULL,
                           pdf_file=NULL,main='',cex_gs=1,cex_gene=0.8,
                           boxmean=FALSE,scale_row=FALSE,color_bar='green2orange',width=NULL,height=NULL,label_z_col=NULL){
  use_gs2gene <- lapply(use_gs2gene,function(x)intersect(x,rownames(use_exp_mat)))
  ori_use_gs2gene <- use_gs2gene
  if(is.null(use_gene)==FALSE) use_gs2gene <- lapply(use_gs2gene,function(x)intersect(x,use_gene)) ## only display use genes
  use_gs2gene <- use_gs2gene[which(unlist(lapply(use_gs2gene,length))>3)]
  print(str(use_gs2gene))
  mat1_draw <- do.call(rbind,lapply(use_gs2gene,function(x){
    #print(str(x))
    x1 <- use_exp_mat[x,]
    if(length(x)>2){
      x2 <- hclust(as.dist(1-cor(t(x1))))
      x1[x2$order,]
    }
  }))
  all_g <- rownames(mat1_draw)
  names(use_p) <- colnames(use_exp_mat)
  tmp1 <- aggregate(t(mat1_draw),list(use_p),'mean')
  mat1_draw_p <- tmp1[,-1]; rownames(mat1_draw_p) <- tmp1[,1]
  uni_p <- unique(use_p)
  mat1_draw_p <- mat1_draw_p[uni_p,]
  mat1_draw_p <- t(mat1_draw_p)
  if(scale_row==TRUE) mat1_draw_p <- t(apply(mat1_draw_p,1,do.std))
  # 1 2 3 4
  mm <- matrix(c(rep(1,12),rep(2,length.out=length(uni_p)),rep(3,length.out=length(uni_p)),rep(4,length.out=round(1.8*ncol(use_z_mat)))),nrow=1)
  if(is.null(width)==TRUE) width=round(ncol(mm)*0.45)
  if(is.null(height)==TRUE) height=round(nrow(mat1_draw)/8)
  print(c(width,height))
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=width,height=height)
  layout(mm)
  # get color bar
  if(color_bar=='spectral'){
    cc <- brewer.pal(11,'Spectral')
    cc1 <- colorRampPalette(c(cc[11:6],'white','white','white',cc[5:1]))(100)
  }else{
    if(color_bar=='blue2red'){
      cc <- brewer.pal(12,'Paired')
      cc1 <- colorRampPalette(c(cc[6],'white',cc[2]))(100)
    }else{ # green2orange
      cc <- brewer.pal(12,'Paired')
      cc1 <- colorRampPalette(c(cc[4],'white',cc[8]))(100)
    }
  }
  bb1 <- seq(0,max(abs(mat1_draw_p)),length.out=50)
  cc2 <- cc1[c(1,25,50,75,100)]
  bb2 <- seq(-max(abs(mat1_draw_p)),max(abs(mat1_draw_p)),length.out=5)
  bb2 <- format(bb2,digits=3)
  # 1
  par(mar=c(8,0,5,0))
  plot(1,xlim=c(0,1),ylim=c(0,1),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col='white')
  yy <- seq(0,1,length.out=nrow(mat1_draw_p))
  text(1,yy,all_g,adj=1,xpd=TRUE,cex=cex_gene)
  gs_len <- unlist(lapply(use_gs2gene,length))
  gs_len_p <- cumsum(gs_len)
  gs_len_p1 <- gs_len/2+c(0,gs_len_p[1:(length(gs_len_p)-1)])
  p1 <- 1-1.05*max(strwidth(all_g)) ## right position for the gs name
  segments(x0=p1,y0=c(0,yy[1+gs_len_p[c(1:c(length(gs_len_p)-1))]]),
           x1=p1,y1=yy[gs_len_p])
  text(p1*0.98,yy[round(gs_len_p1)+1],labels=names(gs_len_p1),adj=1,xpd=TRUE,cex=cex_gs)
  # legend
  xx_left <- seq(0.4,p1,length.out=7)
  xx_mid  <- c(xx_left[1:5]+xx_left[2:6])/2
  rect(xleft=xx_left[1:5],xright=xx_left[2:6],ybottom=-2.75*(yy[2]-yy[1]),ytop=-1.75*(yy[2]-yy[1]),col=cc2,border='black',lwd=0.5,xpd=TRUE)
  text(xx_mid,-3*(yy[2]-yy[1]),adj=1,srt=60,labels=bb2,xpd=TRUE)
  if(scale_row==TRUE) text(xx_mid[3],-(yy[2]-yy[1]),'Z-value (scale by row)',adj=0.5,xpd=TRUE,cex=1) else text(xx_mid[3],-(yy[2]-yy[1]),'Raw input value',xpd=TRUE,cex=1)

  segments(x0=1,x1=p1,y0=c(0-(yy[2]-yy[1])/2,yy[gs_len_p]+(yy[2]-yy[1])/2),y1=c(0-(yy[2]-yy[1])/2,yy[gs_len_p]+(yy[2]-yy[1])/2),lwd=0.3,col='dark grey')

  # 2
  pp <- par()$usr
  par(mar=c(8,0,5,0))
  image(t(mat1_draw_p),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n',bty='n',ylim=c(pp[3],pp[4]),main=main)
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out = 1+ncol(mat1_draw_p));xxx <- xx[1:(length(xx)-1)]+(xx[2]-xx[1])/2
  text(xxx,0-c(pp[4]-pp[3])/nrow(mat1_draw_p),colnames(mat1_draw_p),srt=90,adj=1,xpd=TRUE)
  abline(h=c(0-(yy[2]-yy[1])/2,yy[gs_len_p]+(yy[2]-yy[1])/2),lwd=0.3,col='black')

  # 3
  par(mar=c(8,1.5,5,0))
  plot(1,xlim=c(0,length(uni_p)),ylim=c(pp[3],pp[4]),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col='white')
  yy <- seq(pp[3],pp[4],length.out=nrow(mat1_draw_p))
  yyy <- (yy[2]-yy[1])/2
  yyyy <- yy[gs_len_p[1:(length(gs_len_p)-1)]]+yyy
  rect(xleft = 0,xright=length(uni_p),ybottom=pp[3]-yyy,ytop=pp[4]+yyy)
  segments(x0=0,x1=length(uni_p),y0=yyyy,y1=yyyy)
  segments(x0=0:length(uni_p),x1=0:length(uni_p),y0=pp[3]-yyy,y1=pp[4]+yyy,lwd=0.3,lty=2)
  # box
  get_b <- function(x){
    x <- as.numeric(x)
    x1 <- quantile(x,probs=c(0.25,0.5,0.75,0.0,0.99))
    x1 <- c(x1[1]-1.5*(x1[3]-x1[1]),x1[1],x1[2],x1[3],x1[3]+1.5*(x1[3]-x1[1]),x1[4],x1[5],mean(x1))
    if(x1[1]<x1[6]) x1[1]<-x1[6]
    if(x1[5]>x1[7]) x1[5]<-x1[7]
    x1
  }
  yyyy1 <- c(pp[3]-yyy,yyyy);   yyyy2 <- c(yyyy,pp[4]+yyy);
  for(i in 1:length(use_gs2gene)){
    x1 <- use_exp_mat[ori_use_gs2gene[[i]],]
    if(boxmean==TRUE) x1 <- colMeans(x1)
    x20 <- get_b(x1)
    y1 <- yyyy1[i]; y2 <- yyyy2[i]; y21 <- y2-y1
    y1 <- y1+y21*0.05; y2 <- y2-y21*0.05;
    segments(x0=-0.125,x1=0,y0=c(y1,y2),y1=c(y1,y2),col='dark grey')
    #text(-0.15,c(y1,y2),format(c(x20[6],x20[7]),digits=3),adj=1,xpd=TRUE,cex=0.5)
    text(-0.15,c(y1,y2),format(c(x20[1],x20[5]),digits=3),adj=1,xpd=TRUE,cex=0.5)
    for(j in 1:length(uni_p)){
      x1 <- use_exp_mat[ori_use_gs2gene[[i]],which(use_p==uni_p[j])]
      if(boxmean==TRUE) x1 <- colMeans(x1)
      x2 <- get_b(x1)
      #y3 <- y1+(x2[c(1:5,8)]-x20[6])/(x20[7]-x20[6])*(y2-y1)
      y3 <- y1+(x2[c(1:5,8)]-x20[1])/(x20[5]-x20[1])*(y2-y1)
      rect(xleft=j-1+0.2,xright=j-0.2,ybottom=y3[2],ytop=y3[4],border='light grey',col='light grey')
      segments(x0=j-0.5,x1=j-0.5,y0=y3[1],y1=y3[5],col='grey')
      segments(x0=j-1+0.2,x1=j-0.2,y0=y3[3],y1=y3[3],col='dark grey') ## median
      points(j-0.5,y3[6],pch=16) ## mean
    }
  }
  text(c(1:length(uni_p))-0.5,pp[3]-c(pp[4]-pp[3])/nrow(mat1_draw_p),colnames(mat1_draw_p),srt=90,adj=1,xpd=TRUE)
  pp <- par()$usr
  # 4
  par(mar=c(8,0,5,2))
  plot(1,xlim=c(0,ncol(use_z_mat)),ylim=c(pp[3],pp[4]),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col='white')
  yy <- seq(pp[3],pp[4],length.out=nrow(mat1_draw_p))
  yyy <- (yy[2]-yy[1])/2
  yyyy <- yy[gs_len_p[1:(length(gs_len_p)-1)]]+yyy
  rect(xleft = 0,xright=ncol(use_z_mat),ybottom=pp[3]-yyy,ytop=pp[4]+yyy,xpd=TRUE)
  segments(x0=0,x1=ncol(use_z_mat),y0=yyyy,y1=yyyy)
  segments(x0=0:ncol(use_z_mat),x1=0:ncol(use_z_mat),y0=pp[3]-yyy,y1=pp[4]+yyy,lwd=0.3,lty=1)
  #segments(x0=1:ncol(use_z_mat)-0.5,x1=1:ncol(use_z_mat)-0.5,y0=pp[3]-yyy,y1=pp[4]+yyy,lwd=0.3,lty=2)
  cn <- colnames(use_z_mat); cn <- gsub('-','Vs.',cn)
  text(1:ncol(use_z_mat)-0.5,pp[4]+yyy*2,cn,adj=0,xpd=TRUE,cex=1.2,srt=30)
  u1 <- use_z_mat[names(use_gs2gene),]
  y1 <- yy[gs_len_p1+1]
  y2 <- yy[gs_len_p1-1]
  cc <- matrix(z2col(u1,sig_thre=qnorm(0.9)),nrow=nrow(u1))
  cc[which(cc=='white')] <- 'dark grey'
  for(i in 1:ncol(use_z_mat)){
    #print(get_z2p(u1[,i],use_star=TRUE))
    text(x=i-0.5,y=y1,labels=get_z2p(u1[,i],use_star=TRUE),cex=1,col=cc[,i])
  }
  if(is.null(label_z_col)==FALSE) points(x=label_z_col-0.5,y=pp[3]-yyy*3,cex=2,pch=17,col='red',xpd=TRUE)
  if(is.null(pdf_file)==FALSE) dev.off()
}
##
get_top_ms_tab <- function(ms_tab,top_strategy,top_num,z_col,z_thre=qnorm(1-0.05)){
  top_num <- as.numeric(top_num)
  tmp_ms_tab <- ms_tab[which(abs(ms_tab[,z_col])>=z_thre),]
  tmp_ms_tab_up <- tmp_ms_tab[which(tmp_ms_tab[,z_col]>0),]
  tmp_ms_tab_down <- tmp_ms_tab[which(tmp_ms_tab[,z_col]<0),]
  tmp_ms_tab <- tmp_ms_tab[order(abs(tmp_ms_tab[,z_col]),decreasing=TRUE),]
  tmp_ms_tab_up <- tmp_ms_tab_up[order(tmp_ms_tab_up[,z_col],decreasing=TRUE),]
  tmp_ms_tab_down <- tmp_ms_tab_down[order(tmp_ms_tab_down[,z_col],decreasing=FALSE),]
  if(top_strategy=='Both') ms_tab <- tmp_ms_tab[1:min(top_num,nrow(tmp_ms_tab)),]
  if(top_strategy=='UP') ms_tab <- tmp_ms_tab_up[1:min(top_num,nrow(tmp_ms_tab_up)),]
  if(top_strategy=='DOWN') ms_tab <- tmp_ms_tab_down[1:min(top_num,nrow(tmp_ms_tab_down)),]
  return(ms_tab)
}
####
draw.driverTargetFuncEnrich <- function(use_driver=NULL,target_list=NULL,ms_tab=NULL,transfer_tab=NULL,
                                        use_gs2gene=NULL,exp_mat=NULL,use_ac_mat=NULL,use_p=NULL,use_z_mat=NULL,pdf_file=NULL,main='',cex_driver=0.9,cex_gs=0.5,
                                        scale_row=FALSE,color_bar='green2orange',top_gs_num=10,min_gs_size=50,max_gs_size=2000,Pv_thre=0.1,Pv_adj = 'none',label_z_col=NULL){
  use_driver2gs <- list()
  all_g <- c()
  all_mat1_draw_p <- list()
  all_gs_i <- list()
  use_driver <- rev(use_driver)
  uni_p <- unique(use_p)
  for(each_driver in use_driver){
    print(each_driver)
    ori_tar     <- target_list[[each_driver]]$target
    ori_tar <- intersect(ori_tar,rownames(exp_mat))
    use_tar <- ori_tar
    if(is.null(transfer_tab)==FALSE) use_tar <- get_name_transfertab(ori_tar,transfer_tab = transfer_tab,ignore_version = T)
    use_tar_wei <- target_list[[each_driver]]$MI*sign(target_list[[each_driver]]$spearman)
    res1 <- funcEnrich.Fisher(use_tar,gs2gene=use_gs2gene,min_gs_size = min_gs_size,max_gs_size = max_gs_size,Pv_thre=1,Pv_adj = Pv_adj)
    if(nrow(res1)==0 | is.na(res1[1,1])==TRUE) next
    if(nrow(res1)==1){
      all_g2s <- list(unlist(strsplit(res1$Intersected_items,';')))
      names(all_g2s) <- res1$`#Name`
    }
    if(nrow(res1)>1){
      w1 <- which(res1$Ori_P<=Pv_thre)
      if(length(w1)>3) res1 <- res1[w1,] else res1 <- res1[1:3,]
      all_g2s <- lapply(res1$Intersected_items,function(x1)unlist(strsplit(x1,';')))
      names(all_g2s) <- res1$`#Name`
      mat1 <- t(list2mat(all_g2s))
      h_gs <- hclust(dist(mat1,method='binary'))
      gs_cluster <- cutree(h_gs,h=0) ## remove gene sets with 100% percentage similarity
      if(length(unique(gs_cluster))>3){
        all_g2s <- all_g2s[sapply(unique(gs_cluster),function(x){
          n1 <- names(gs_cluster)[which(gs_cluster==x)]
          n1[which.min(res1[n1,]$Ori_P)]
        })]
      }
      all_g2s <- all_g2s[order(res1[names(all_g2s),'Ori_P'])[1:min(top_gs_num,length(all_g2s))]]
    }
    gs_m <- lapply(all_g2s,function(x){
      x0 <- which(use_tar %in% x)
      x1 <- exp_mat[ori_tar[x0],]*use_tar_wei[x0]
      if(length(x0)==1) x1 else colMeans(x1)
    })
    gs_m <- do.call(rbind,gs_m)
    tmp1 <- aggregate(t(gs_m),list(use_p),'mean')
    mat1_draw_p <- tmp1[,-1];
    if(is.null(dim(mat1_draw_p))==TRUE){
      mat1_draw_p <- as.data.frame(mat1_draw_p)
      names(mat1_draw_p) <- names(all_g2s)
    }
    rownames(mat1_draw_p) <- tmp1[,1]
    mat1_draw_p <- mat1_draw_p[uni_p,]
    mat1_draw_p <- t(mat1_draw_p)
    if(scale_row==TRUE) mat1_draw_p <- t(apply(mat1_draw_p,1,std))
    w1 <- which(is.na(mat1_draw_p[,1])==FALSE)
    mat1_draw_p <- mat1_draw_p[w1,]
    #
    if(is.null(dim(mat1_draw_p))==FALSE) mat1_draw_p <- mat1_draw_p[rev(1:nrow(mat1_draw_p)),]
    all_g <- c(all_g,rownames(mat1_draw_p))
    use_driver2gs[[each_driver]] <- rownames(mat1_draw_p)
    all_mat1_draw_p[[each_driver]] <- mat1_draw_p
    all_gs_i[[each_driver]] <- res1[rownames(mat1_draw_p),c('Num_item','Ori_P','Num_list_item','Intersected_items')]
  }
  mat1_draw_p <- do.call(rbind,all_mat1_draw_p)
  all_gi <- do.call(rbind,all_gs_i)
  # 1 2 3 4
  mm <- matrix(c(rep(1,20),rep(2,length.out=0.8*length(uni_p)),rep(3,length.out=length(uni_p)),rep(4,length.out=round(2.5*ncol(use_z_mat)))),nrow=1)
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=round(ncol(mm)*0.25),height=round(nrow(mat1_draw_p)/10))
  layout(mm)
  # get color bar
  if(color_bar=='spectral'){
    cc <- brewer.pal(11,'Spectral')
    cc1 <- colorRampPalette(c(cc[11:6],'white','white','white',cc[5:1]))(100)
  }else{
    if(color_bar=='blue2red'){
      cc <- brewer.pal(12,'Paired')
      cc1 <- colorRampPalette(c(cc[2],'white',cc[6]))(100)
    }else{ # green2orange
      cc <- brewer.pal(12,'Paired')
      cc1 <- colorRampPalette(c(cc[4],'white',cc[8]))(100)
    }
  }
  bb1 <- seq(0,max(abs(mat1_draw_p)),length.out=50)
  cc2 <- cc1[c(1,25,50,75,100)]
  bb2 <- seq(-max(abs(mat1_draw_p)),max(abs(mat1_draw_p)),length.out=5)
  bb2 <- format(bb2,digits=3)
  # 1
  par(mar=c(8,0,5,0))
  plot(1,xlim=c(0,1.1),ylim=c(0,1),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col='white')
  yy <- seq(0,1,length.out=nrow(mat1_draw_p))
  text(1.1,yy,paste0('(',all_gi$Num_list_item,')'),adj=1,xpd=TRUE,cex=cex_gs,xpd=TRUE) # number of intersection
  p1 <- 1.1-quantile(strwidth(paste0('(',all_gi$Num_list_item,')'),cex=cex_gs),probs=1)
  cc<-z2col(qnorm(1-all_gi$Ori_P),sig_thre=qnorm(0.9),blue_col = brewer.pal(9,'Blues')[5:8],red_col = brewer.pal(9,'Reds')[5:8])
  cc[which(cc=='white')] <- 'dark grey'
  text(p1,yy,format(all_gi$Ori_P,digits=2),adj=1,xpd=TRUE,cex=cex_gs,xpd=TRUE,col=cc) # pv
  p1 <- p1-1.15*quantile(strwidth(format(all_gi$Ori_P,digits=2),cex=cex_gs),probs=1)
  text(p1,yy,all_g,adj=1,xpd=TRUE,cex=cex_gs)
  segments(y0=0,y1=1,x0=p1,x1=p1,lwd=0.3,col='dark grey')
  gs_len <- unlist(lapply(use_driver2gs,length))
  gs_len_p <- cumsum(gs_len)
  if(length(gs_len_p)>1) gs_len_p1 <- gs_len/2+c(0,gs_len_p[1:(length(gs_len_p)-1)]) else gs_len_p1 <- gs_len/2
  p1 <- p1-1.005*quantile(strwidth(all_g,cex=cex_gs),probs=1) ## right position for the driver name
  segments(x0=p1,y0=c(0,yy[1+gs_len_p[c(1:c(length(gs_len_p)-1))]]),
           x1=p1,y1=yy[gs_len_p])
  segments(x0=p1,x1=1.1,y0=c(0-(yy[2]-yy[1])/2,yy[gs_len_p]+(yy[2]-yy[1])/2),y1=c(0-(yy[2]-yy[1])/2,yy[gs_len_p]+(yy[2]-yy[1])/2),lwd=0.3,col='dark grey')
  text(p1*0.99,yy[round(gs_len_p1)+1],labels=paste0(ms_tab[names(gs_len_p1),'gene_label'],'\n(Size:',ms_tab[names(gs_len_p1),]$Size,')'),adj=1,xpd=TRUE,cex=cex_driver)

  # 2
  par(mar=c(8,0,5,0))
  pp <- par()$usr
  graphics::image(t(mat1_draw_p),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n',bty='n',ylim=c(pp[3],pp[4]),main=main,xpd=TRUE)
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out = 1+ncol(mat1_draw_p));xxx <- xx[1:(length(xx)-1)]+(xx[2]-xx[1])/2
  text(xxx,0-c(pp[4]-pp[3])/nrow(mat1_draw_p),colnames(mat1_draw_p),srt=90,adj=1,xpd=TRUE)
  # legend
  xx_left <- seq(pp[1],pp[2],length.out=10)
  xx_mid  <- c(xx_left[3:7]+xx_left[4:8])/2
  rect(xleft=xx_left[3:7],xright=xx_left[4:8],ybottom=pp[4]+(pp[4]-pp[3])/nrow(mat1_draw_p),
       ytop=pp[4]+2*(pp[4]-pp[3])/nrow(mat1_draw_p),col=cc2,border='black',lwd=0.5,xpd=TRUE)
  text(xx_mid,pp[4]+(pp[4]-pp[3])/nrow(mat1_draw_p),adj=1,srt=60,labels=bb2,xpd=TRUE,cex=0.5)
  if(scale_row==TRUE) text(xx_mid[3],pp[4]+2.5*(pp[4]-pp[3])/nrow(mat1_draw_p),'Z-value\n(scale by row)',adj=0.5,xpd=TRUE,cex=0.5) else text(xx_mid[3],pp[4]+2.5*(pp[4]-pp[3])/nrow(mat1_draw_p),'Raw value',xpd=TRUE,cex=0.8)
  abline(h=c(0-(yy[2]-yy[1])/2,yy[gs_len_p]+(yy[2]-yy[1])/2),lwd=0.3,col='black')
  #
  # 3
  par(mar=c(8,1.5,5,0))
  plot(1,xlim=c(0,length(uni_p)),ylim=c(pp[3],pp[4]),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col='white')
  yy <- seq(pp[3],pp[4],length.out=nrow(mat1_draw_p))
  yyy <- (yy[2]-yy[1])/2
  yyyy <- yy[gs_len_p[1:(length(gs_len_p)-1)]]+yyy
  rect(xleft = 0,xright=length(uni_p),ybottom=pp[3]-yyy,ytop=pp[4]+yyy)
  segments(x0=0,x1=length(uni_p),y0=yyyy,y1=yyyy)
  segments(x0=0:length(uni_p),x1=0:length(uni_p),y0=pp[3]-yyy,y1=pp[4]+yyy,lwd=0.3,lty=2)
  # box
  get_b <- function(x){
    x <- as.numeric(x)
    x1 <- quantile(x,probs=c(0.25,0.5,0.75,0.0,0.99))
    x1 <- c(x1[1]-1.5*(x1[3]-x1[1]),x1[1],x1[2],x1[3],x1[3]+1.5*(x1[3]-x1[1]),x1[4],x1[5],mean(x1))
    if(x1[1]<x1[6]) x1[1]<-x1[6]
    if(x1[5]>x1[7]) x1[5]<-x1[7]
    x1
  }
  yyyy1 <- c(pp[3]-yyy,yyyy);   yyyy2 <- c(yyyy,pp[4]+yyy);
  for(i in 1:length(use_driver2gs)){
    x1 <- use_ac_mat[names(use_driver2gs)[i],]
    x20 <- get_b(x1)
    y1 <- yyyy1[i]; y2 <- yyyy2[i]; y21 <- y2-y1
    y1 <- y1+y21*0.05; y2 <- y2-y21*0.05;
    segments(x0=-0.125,x1=0,y0=c(y1,y2),y1=c(y1,y2),col='dark grey')
    text(-0.15,c(y1,y2),format(c(x20[6],x20[7]),digits=2),adj=1,xpd=TRUE,cex=0.4) ##
    for(j in 1:length(uni_p)){
      x1 <- use_ac_mat[names(use_driver2gs)[i],which(use_p==uni_p[j])]
      x2 <- get_b(x1)
      y3 <- y1+(x2[c(1:5,8)]-x20[6])/(x20[7]-x20[6])*(y2-y1)
      rect(xleft=j-1+0.2,xright=j-0.2,ybottom=y3[2],ytop=y3[4],border='light grey',col='light grey')
      segments(x0=j-0.5,x1=j-0.5,y0=y3[1],y1=y3[5],col='grey')
      segments(x0=j-1+0.2,x1=j-0.2,y0=y3[3],y1=y3[3],col='dark grey') ## median
      points(j-0.5,y3[6],pch=16) ## mean
    }
  }
  text(c(1:length(uni_p))-0.5,pp[3]-c(pp[4]-pp[3])/nrow(mat1_draw_p),colnames(mat1_draw_p),srt=90,adj=1,xpd=TRUE)
  pp <- par()$usr
  # 4
  par(mar=c(8,0,5,5))
  plot(1,xlim=c(0,ncol(use_z_mat)),ylim=c(pp[3],pp[4]),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col='white')
  yy <- seq(pp[3],pp[4],length.out=nrow(mat1_draw_p))
  yyy <- (yy[2]-yy[1])/2
  yyyy <- yy[gs_len_p[1:(length(gs_len_p)-1)]]+yyy
  rect(xleft = 0,xright=ncol(use_z_mat),ybottom=pp[3]-yyy,ytop=pp[4]+yyy,xpd=TRUE)
  segments(x0=0,x1=ncol(use_z_mat),y0=yyyy,y1=yyyy)
  segments(x0=0:ncol(use_z_mat),x1=0:ncol(use_z_mat),y0=pp[3]-yyy,y1=pp[4]+yyy,lwd=0.3,lty=1)
  #segments(x0=1:ncol(use_z_mat)-0.5,x1=1:ncol(use_z_mat)-0.5,y0=pp[3]-yyy,y1=pp[4]+yyy,lwd=0.3,lty=2)
  cn <- colnames(use_z_mat); cn <- gsub('-','Vs.',cn)
  text(1:ncol(use_z_mat)-0.5,pp[4]+yyy*2,cn,adj=0,xpd=TRUE,cex=1,srt=30)
  u1 <- use_z_mat[names(use_driver2gs),]
  y1 <- yy[gs_len_p1+1]
  y2 <- yy[gs_len_p1-1]
  if(length(use_driver2gs)==1){
    cc <- z2col(as.matrix(u1),sig_thre=qnorm(0.9),blue_col = brewer.pal(9,'Blues')[5:8],red_col = brewer.pal(9,'Reds')[5:8])
    cc[which(cc=='white')] <- 'dark grey'
    for(i in 1:ncol(use_z_mat)){
      text(x=i-0.5,y=y1,labels=get_z2p(u1[i],use_star=TRUE),cex=0.8,col=cc[i])
    }
  }else{
    cc <- matrix(z2col(as.matrix(u1),sig_thre=qnorm(0.9),blue_col = brewer.pal(9,'Blues')[5:8],red_col = brewer.pal(9,'Reds')[5:8]),nrow=nrow(u1))
    cc[which(cc=='white')] <- 'dark grey'
    for(i in 1:ncol(use_z_mat)){
      text(x=i-0.5,y=y1,labels=get_z2p(u1[,i],use_star=TRUE),cex=0.8,col=cc[,i])
    }
  }
  if(is.null(label_z_col)==FALSE) points(x=label_z_col-0.5,y=pp[3]-yyy*3,cex=2,pch=17,col='red',xpd=TRUE)
  if(is.null(pdf_file)==FALSE) dev.off()
}

#### merge master table
merge_ms_tab <- function(x){
  x1 <- rownames(x[[1]])
  w1 <- grep('.*_D[A|E]$',colnames(x[[1]]))
  w2 <- setdiff(1:ncol(x[[1]]),w1)
  w21 <- w2[which(w2<min(w1))]
  w22 <- w2[which(w2>max(w1))]
  x2 <- lapply(x,function(k){
    w1 <- grep('.*_D[A|E]$',colnames(k))
    x3 <- k[x1,w1]
    x4 <- x3[,grep('Z.',colnames(x3))]
    c1 <- unique(gsub('Z.(.*)_D[A|E]','\\1',colnames(x4)))
    c2 <- list(paste0('Z.',c1),paste0('AveExpr.',c1),paste0('logFC.',c1),paste0('P.Value.',c1))
    c3 <- lapply(c2,function(k){c(paste0(k,'_DA'),paste0(k,'_DE'))})
    c3 <- lapply(c3,function(k)intersect(k,colnames(x3)))
    lapply(c3,function(k)x3[,k])
  })
  x3 <- lapply(names(x2),function(k){
    x3 <- x2[[k]];
    x3 <- lapply(x3,function(k1){
      colnames(k1) <- paste0(gsub('ms_','',k),'|',colnames(k1))
      k1
    })
    x3
  })
  names(x3) <- names(x2)
  a1 <- list()
  for(i in 1:length(x3[[1]])){
    x4 <- do.call(cbind,lapply(x3,function(k){k[[i]]}))
    n4 <- unlist(lapply(x3,function(k)colnames(k[[i]])))
    colnames(x4) <- n4
    a1[[i]] <- x4
  }
  a1 <- do.call(cbind,a1)
  a2 <- cbind(x[[1]][,w21],a1,x[[1]][,w22])
  return(a2)
}
#
draw.topdriver2gs <- function(gene_list=NULL,gs2gene=NULL,clust_cex=0.8,low_Sum=5,low_K=3,outlier_cex=0.5,recluster_row=TRUE,plot=TRUE){
  gs2gene <- gs2gene[which(unlist(lapply(gs2gene,length))>0)]
  mat1 <- matrix(0,nrow=length(gene_list),ncol=length(gs2gene))
  for(i in 1:length(gene_list)){
    for(j in 1:length(gs2gene)){
      g1 <- gene_list[[i]]
#      print(g1)
      g1 <- gsub('(.*)_.*','\\1',g1)
      mat1[i,j] <- length(intersect(g1,gs2gene[[j]]))
    }
  }
  rownames(mat1) <- names(gene_list)
  colnames(mat1) <- names(gs2gene)
#  print(mat1)
  c1 <- hclust(dist(t(mat1)))
  if(recluster_row==TRUE) mat1 <- mat1[,rev(c1$order)]
  w1 <- colSums(mat1)
  t1 <- mat1[,which(w1>=low_Sum)]
  if(length(which(w1>=low_Sum))<2){message('lower the lowSum and re-try!');return(FALSE)}
  if(plot==FALSE) return(t1)
  par(mar=c(2,max(nchar(colnames(t1)))/2+5,8,6))
  image(t1,col=c('white',colorRampPalette(brewer.pal(8,'Reds'))(length(unique(as.numeric(t1)))-1)),bty='n',xaxt='n',yaxt='n')
  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
  xx <- seq(pp[1],pp[2],length.out = nrow(t1)+1)
  yy <- seq(pp[3],pp[4],length.out = ncol(t1)+1)
  xxx <- (xx[1:(length(xx)-1)]+xx[2:length(xx)])/2
  yyy <- (yy[1:(length(yy)-1)]+yy[2:length(yy)])/2
  abline(h=yy);abline(v=xx)
  text(pp[1],yyy,colnames(t1),adj=1,xpd=TRUE,col=1,cex=clust_cex)
  text(xxx,pp[4]+max(strheight(rownames(t1),units='inches',cex=1))/12,rownames(t1),adj=0,srt=30,xpd=TRUE,
       col=1,cex=clust_cex)
  for(i in 1:nrow(t1)){
    for(j in 1:ncol(t1)){
      v1 <- t1[i,j]
      if(v1==0) next
      if(v1>low_K){
        text(xxx[i],yyy[j],v1,cex=clust_cex)
      }else{
        g1 <- gene_list[[rownames(t1)[i]]]
        g2 <- gsub('(.*)_.*','\\1',g1)
        v2 <- g1[which(g2 %in% gs2gene[[colnames(t1)[j]]])]
        v2 <- paste(v2,collapse='\n')
        text(xxx[i],yyy[j],v2,cex=outlier_cex)
      }
    }
  }
  return(t1)
}
##
draw.heatmap.local <- function(mat,inner_line=FALSE,out_line=TRUE,n=20,col_srt=60,display_text_mat=NULL,row_cex=1,column_cex=1,text_cex=1){
  mat1 <- mat[nrow(mat):1,]
  if(is.null(display_text_mat)==FALSE) display_text_mat <- display_text_mat[nrow(display_text_mat):1,]
  cc1 <- colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(n*2)
  bb1 <- seq(0,max(abs(mat1)),length.out=n)
  if(out_line==TRUE) image(t(mat1),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n')
  if(out_line==FALSE) image(t(mat1),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out=1+ncol(mat1))
  yy <- seq(pp[3],pp[4],length.out=1+nrow(mat1))
  xxx <- xx[1:(length(xx)-1)]/2+xx[2:length(xx)]/2
  yyy <- yy[1:(length(yy)-1)]/2+yy[2:length(yy)]/2
  if(inner_line==TRUE){
    abline(v=xx)
    abline(h=yy)
  }
  text(pp[1]-max(strwidth(rownames(mat1))/10),yyy,rownames(mat1),adj=1,xpd=TRUE,cex=row_cex)
  text(xxx,pp[4]+max(strheight(colnames(mat1))/2)+max(strheight(colnames(mat1))/10),colnames(mat1),
       adj=0,xpd=TRUE,srt=col_srt,cex=column_cex)
  if(is.null(display_text_mat)==FALSE){
    for(i in 1:nrow(display_text_mat)){
      for(j in 1:ncol(display_text_mat)){
        text(xxx[j],yyy[i],display_text_mat[i,j],cex=text_cex)
      }
    }
  }
  return(TRUE)
}
##
draw.heatmap.TWO <- function(mat1,mat2,scale='row',column_ratio=c(ncol(mat1),ncol(mat2)),
                             cex_row=0.8,cex_column=1,cex_label=1,subtype2gene=NULL,
                             phenotype_info1=NULL,use_phe1=NULL,
                             phenotype_info2=NULL,use_phe2=NULL,col_names=TRUE){
  if(scale=='row'){
    mat1 <- t(apply(mat1,1,std))
    mat2 <- t(apply(mat2,1,std))
  }
  r1 <- column_ratio ## row ratio
  r2 <- c(1,15) ## column ratio
  l1 <- c(3,4)
  l2 <- c(1,2)
  b1 <- unlist(sapply(1:length(r1),function(i){rep(l1[i],r1[i])}))
  b2 <- unlist(sapply(1:length(r1),function(i){rep(l2[i],r1[i])}))
  layout(matrix(c(rep(b1,r2[1]),rep(b2,r2[2])),nrow=sum(r2),byrow=TRUE))
  # 1
  par(mar=c(10,6,0,0))
  nn <- 20
  cc1 <- colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(nn*2)
  bb1 <- seq(0,max(abs(mat1)),length.out=nn)
  image(t(mat1),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n')
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out=1+ncol(mat1))
  yy <- seq(pp[3],pp[4],length.out=1+nrow(mat1))
  xxx <- xx[1:(length(xx)-1)]/2+xx[2:length(xx)]/2
  yyy <- yy[1:(length(yy)-1)]/2+yy[2:length(yy)]/2
  text(pp[1],yyy,rownames(mat1),pos=2,xpd=TRUE,cex=cex_row)
  if(is.null(subtype2gene)==FALSE){
    subtype2gene <- lapply(subtype2gene,function(x)intersect(x,rownames(mat1)))
    s1 <- unlist(lapply(subtype2gene,length))
    s11 <- cumsum(s1)
    abline(h=yy[s11+1])
  }
  if(col_names==TRUE) text(xxx,yy[1]-yyy[2]+yyy[1],colnames(mat1),srt=60,xpd=TRUE,adj=1)
  # 2
  par(mar=c(10,4,0,4))
  nn <- 20
  cc1 <- colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(nn*2)
  bb1 <- seq(0,max(abs(mat2)),length.out=nn)
  image(t(mat2),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n')
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out=1+ncol(mat2))
  yy <- seq(pp[3],pp[4],length.out=1+nrow(mat2))
  xxx <- xx[1:(length(xx)-1)]/2+xx[2:length(xx)]/2
  yyy <- yy[1:(length(yy)-1)]/2+yy[2:length(yy)]/2
  text(pp[1],yyy,rownames(mat2),pos=2,xpd=TRUE,cex=cex_row)
  if(is.null(subtype2gene)==FALSE){
    subtype2gene <- lapply(subtype2gene,function(x)intersect(x,rownames(mat2)))
    s1 <- unlist(lapply(subtype2gene,length))
    s11 <- cumsum(s1)
    abline(h=yy[s11+1])
  }
  if(col_names==TRUE) text(xxx,yy[1]-yyy[2]+yyy[1],colnames(mat2),srt=60,xpd=TRUE,adj=1)
  # 3
  par(mar=c(0,6,2,0))
  plot(1,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='')
  pp <- par()$usr
  if(is.null(phenotype_info1)==FALSE){
    yy <- seq(0,1,length.out=1+length(use_phe1))
    xx <- seq(pp[1],pp[2],length.out=1+ncol(mat1))
    pp <- par()$usr
    #abline(h=yy);abline(v=c(pp[1:2]))
    for(i in 1:length(use_phe1)){
      use_p <- phenotype_info1[,use_phe1[i]]
      #use_c <- colorRampPalette(sample(brewer.pal(8,'Paired')))(length(use_p));
      use_c <- get.class.color(use_p)
      names(use_c)<-use_p
      rect(ybottom=yy[i],ytop=yy[i+1],xleft=xx[1:(length(xx)-1)],xright=xx[2:length(xx)],
           col=use_c[use_p],border=NA,xpd=TRUE)
      s1 <- table(use_p)[unique(use_p)]
      s11 <- cumsum(s1)
      xxx <- c(pp[1],xx[s11+1])
      text(xxx[1:(length(xxx)-1)]/2+xxx[2:(length(xxx))]/2,yy[i]/2+yy[i+1]/2,names(s1),adj=0.5,cex=cex_label,xpd=TRUE)
    }
  }
  ##
  # 4
  par(mar=c(0,4,2,4))
  plot(1,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='')
  pp <- par()$usr
  if(is.null(phenotype_info2)==FALSE){
    yy <- seq(0,1,length.out=1+length(use_phe2))
    xx <- seq(pp[1],pp[2],length.out=1+ncol(mat2))
    pp <- par()$usr
    #abline(h=yy);abline(v=c(pp[1:2]))
    for(i in 1:length(use_phe2)){
      use_p <- phenotype_info2[,use_phe2[i]]
     # use_c <- colorRampPalette(sample(brewer.pal(8,'Paired')))(length(use_p));
      use_c <- get.class.color(use_p)
      names(use_c)<-use_p
      rect(ybottom=yy[i],ytop=yy[i+1],xleft=xx[1:(length(xx)-1)],xright=xx[2:length(xx)],
           col=use_c[use_p],border=NA,xpd=TRUE)
      s1 <- table(use_p)[unique(use_p)]
      s11 <- cumsum(s1)
      xxx <- c(pp[1],xx[s11+1])
      text(xxx[1:(length(xxx)-1)]/2+xxx[2:(length(xxx))]/2,yy[i]/2+yy[i+1]/2,names(s1),adj=0.5,cex=cex_label)
    }
  }
  ##
  return(TRUE)
}
##
draw.heatmap.top <- function(mat,use_p=NULL,sample_label=NULL,subtype2gene=NULL,subtype2gs=NULL,int_gene_list=NULL,
                             column_ratio=c(6,10,2),sample_subtype_order = NULL,int_gene_list_label=NULL,std=TRUE,gene_cex=1){
  ##
  subtype2gene <- lapply(subtype2gene,function(x){intersect(x,rownames(mat))})
  subtype2gene <- subtype2gene[which(unlist(lapply(subtype2gene,length))>0)]
  subtype2gene <- subtype2gene[rev(names(subtype2gene))]
  gene_order <- unlist(subtype2gene)
  mat1 <- mat[gene_order,names(use_p)]
  if(std==TRUE) mat1 <- t(apply(mat1,1,std))
  if(is.null(int_gene_list_label)==TRUE) int_gene_list_label<-int_gene_list
  names(int_gene_list_label) <- int_gene_list
  int_gene_list <- intersect(int_gene_list,rownames(mat1))
  ## 0|2|0
  ## 3|1|4
  r1 <- column_ratio ## row ratio
  r2 <- c(1,10) ## column ratio
  l1 <- c(0,2,0)
  l2 <- c(3,1,4)
  b1 <- unlist(sapply(1:length(r1),function(i){rep(l1[i],r1[i])}))
  b2 <- unlist(sapply(1:length(r1),function(i){rep(l2[i],r1[i])}))
  layout(matrix(c(rep(b1,r2[1]),rep(b2,r2[2])),nrow=sum(r2),byrow=TRUE))
  # 1
  par(mar=c(6,0,0,0))
  nn <- 20
  cc1 <- colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(nn*2)
  bb1 <- seq(0,max(abs(mat1)),length.out=nn)
  image(t(mat1),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n')
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out=1+ncol(mat1))
  yy <- seq(pp[3],pp[4],length.out=1+nrow(mat1))
  s1 <- unlist(lapply(subtype2gene,length))
  s2 <- table(use_p)[unique(use_p)]
  s11 <- cumsum(s1)
  s22 <- cumsum(s2)
  abline(h=yy[s11+1])
  abline(v=xx[s22+1])
  text(xx[1:(length(xx)-1)]/2+xx[2:(length(xx))]/2,pp[3]-yy[2]+yy[1],sample_label,srt=90,xpd=TRUE,adj=1)
  pp0 <- pp
  if(length(sample_subtype_order)>0){
    xxx <- c(pp[1],xx[s22+1])
    xxx <- sort(xxx[1:(length(xxx)-1)]/2+xxx[2:(length(xxx))]/2)
    yyy <- c(pp[3],yy[s11+1])
    yyy <- yyy[1:(length(yyy)-1)]/2+yyy[2:(length(yyy))]/2
    yyy <- sort(yyy,decreasing = TRUE)
    text(x=xxx[sample_subtype_order],
         y=yyy,
         rev(names(subtype2gene)),xpd=TRUE,cex=1)
  }

  # 2
  par(mar=c(0,0,2,0))
  plot(1,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='')
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out=1+ncol(mat1))
  yy <- seq(pp[3],pp[4],length.out=1+nrow(mat1))
  yyy <- c(pp[1],xx[s22+1])
  rect(xleft=yyy[1:(length(yyy)-1)],xright=yyy[2:(length(yyy))],ybottom=pp[3],ytop=0.3,col=get.class.color(names(s22)),xpd=TRUE)
  text(yyy[1:(length(yyy)-1)]/2+yyy[2:(length(yyy))]/2,0.5,names(s22),xpd=TRUE,cex=1.2)
  # 3
  par(mar=c(6,0,0,0))
  plot(1,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='')
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out=1+ncol(mat1))
  yy <- seq(pp[3],pp[4],length.out=1+nrow(mat1))
  yyy <- c(pp[3],yy[s11+1])
  rect(ybottom=yyy[1:(length(yyy)-1)],ytop=yyy[2:(length(yyy))],xleft=1,xright=pp[2],col=get.class.color(names(s11)),xpd=TRUE)
  ##
  if(is.null(subtype2gs)==TRUE){
     text(x=0.99,y=yyy[1:(length(yyy)-1)]/2+yyy[2:(length(yyy))]/2,names(s11),xpd=TRUE,adj=1)
  }else{
    mm <- max(unlist(lapply(subtype2gs,length)))*2+1
    for(i in names(subtype2gs)){
      text(x=0.99,y=seq(yy[s11[i]-s1[i]/2-4],length.out = nrow(subtype2gs[[i]]),by=(pp[4]-pp[3])/(0.5*mm*length(subtype2gs))),
           paste0(subtype2gs[[i]]$`#Name`,' p=',format(subtype2gs[[i]]$`Ori_P`,digits=2,scientific=TRUE),sep=''),xpd=TRUE,adj=1)
    }
  }
  ## 4
  if(length(int_gene_list)>0){
    par(mar=c(6,0,0,0))
    plot(1,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='')
    pp <- par()$usr
    yy <- seq(pp[3],pp[4],length.out=1+nrow(mat1))
    w0 <- which(rownames(mat1) %in% int_gene_list)
    yyy <- c(pp[4]-pp[3])/(nrow(mat1)/3)
    s1 <- unlist(lapply(subtype2gene,length))
    s11 <- cumsum(s1)
    hh=c(0,yy[s11+1])
    for(i in 1:length(subtype2gene)){
      w1 <- which(subtype2gene[[names(s1)[i]]] %in% int_gene_list)
      w2 <- which(rownames(mat1) %in% subtype2gene[[names(s1)[i]]][w1])
      if(length(w2)>0){
        yy1 <- seq(hh[i]+yyy,hh[i+1]-yyy,length.out=length(w1))
        segments(x0=0,x1=0.1,y0=yy[w2],y1=yy1,xpd=TRUE)
        text(0.12,yy1,int_gene_list_label[rownames(mat1)[w2]],adj=0,cex=gene_cex,xpd=TRUE)
      }
    }
  }else{
    par(mar=c(6,0,0,0))
    plot(1,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='',xaxs='i',yaxs='i')
    pp <- par()$usr
    yy <- seq(0,1,length.out=1+nrow(mat1))
    text(0,yy[1:(length(yy)-1)]/2+yy[2:length(yy)]/2,rownames(mat1),pos=4,xpd=TRUE,cex=gene_cex)
  }
  ##
  return(TRUE)
}
##
draw.clustComp.2 <- function(pred_label, obs_label,strategy='ARI',
                           use_col=TRUE,low_K=5,
                           highlight_clust=NULL,
                           main=NULL,clust_cex=1,outlier_cex=0.3,pred_srt=30) {
  if(is.null(names(pred_label))==TRUE & is.null(names(obs_label))==FALSE){
    message('The names of pred_label will use the order of obs_label')
    names(pred_label) <- names(obs_label)
  }
  if(is.null(names(obs_label))==TRUE & is.null(names(pred_label))==FALSE){
    message('The names of obs_label will use the order of pred_label')
    names(obs_label) <- names(pred_label)
  }
  if(is.null(names(obs_label))==TRUE & is.null(names(pred_label))==TRUE){
    message('Assume pred_label and obs_label have the same order!')
    names(obs_label) <- as.character(1:length(obs_label))
    names(pred_label) <- names(obs_label)
  }
  nn <- names(pred_label)
  k1 <- get_clustComp(pred_label,obs_label,strategy=strategy)
  if(is.null(main)==TRUE){
    mm <- sprintf('%s:%s',strategy,format(k1,digits=4),
                  format(k1,digits=4))
  }else{
    mm <- main
  }
  t1 <- table(list(pred_label[nn],obs_label[nn]))
  layout(1)
  par(mar=c(2,12,10,6))
  if(use_col==TRUE){
    image(t1,col=c('white',colorRampPalette(brewer.pal(8,'Reds'))(length(unique(as.numeric(t1)))-1)),bty='n',xaxt='n',yaxt='n',
          main=mm)
  }else{
    image(t1,bty='n',xaxt='n',yaxt='n',
          main=mm,col='white')
  }

  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
  xx <- seq(pp[1],pp[2],length.out = nrow(t1)+1)
  yy <- seq(pp[3],pp[4],length.out = ncol(t1)+1)
  xxx <- (xx[1:(length(xx)-1)]+xx[2:length(xx)])/2
  yyy <- (yy[1:(length(yy)-1)]+yy[2:length(yy)])/2

  abline(h=yy);abline(v=xx)
  text(pp[1],yyy,colnames(t1),adj=1,xpd=TRUE,col=ifelse(colnames(t1) %in% highlight_clust,2,1),cex=clust_cex)
  text(xxx,pp[4]+max(strheight(rownames(t1),units='inches',cex=1))/12,rownames(t1),adj=0,srt=0,xpd=TRUE,
       col=ifelse(rownames(t1) %in% highlight_clust,2,1),cex=clust_cex,srt=pred_srt)

  for(i in 1:nrow(t1)){
    for(j in 1:ncol(t1)){
      v1 <- t1[i,j]
      if(v1==0) next
      if(v1>low_K){
        text(xxx[i],yyy[j],v1,cex=clust_cex)
      }else{
        v2 <- names(obs_label)[which(pred_label==rownames(t1)[i] & obs_label==colnames(t1)[j])]
        v2 <- paste(v2,collapse='\n')
        text(xxx[i],yyy[j],v2,cex=outlier_cex)
      }
    }
  }
  return(t1)
}

draw.multiBox.TWO <- function(label_info=NULL,box1=NULL,box2=NULL,main_box1='',main_box2='',
                              pval_box1=NULL,pval_box2=NULL,pv_col_box1=1,pv_col_box2=1){
  # box
  get_b <- function(x){
    if(length(x)==0) return(NULL)
    x <- as.numeric(x)
    x1 <- quantile(x,probs=c(0.25,0.5,0.75,0.0,0.99))
    x1 <- c(x1[1]-1.5*(x1[3]-x1[1]),x1[1],x1[2],x1[3],x1[3]+1.5*(x1[3]-x1[1]),x1[4],x1[5],mean(x1))
    if(x1[1]<x1[6]) x1[1]<-x1[6]
    if(x1[5]>x1[7]) x1[5]<-x1[7]
    x1
  }
  uni_p <- c(names(box1[[1]]),names(box2[[1]]))
  l1 <- length(names(box1[[1]]))
  l2 <- length(names(box2[[1]]))
  n1 <- length(label_info)
  plot(1,xlim=c(0,length(uni_p)),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col='white')
  pp <- par()$usr
  yy <- seq(pp[3],pp[4],length.out=n1)
  yyy <- (yy[2]-yy[1])/2
  yyyy <- yy-yyy
  rect(xleft = 0,xright=length(uni_p),ybottom=pp[3]-yyy,ytop=pp[4]+yyy,xpd=TRUE)
  segments(x0=0,x1=length(uni_p),y0=yyyy,y1=yyyy,xpd=TRUE)
  segments(x0=0:length(uni_p),x1=0:length(uni_p),y0=pp[3]-yyy,y1=pp[4]+yyy,lwd=0.3,lty=2,xpd=TRUE)
  text(pp[1]-strwidth(label_info)*0.1,yy,label_info,adj=1,xpd=TRUE)
  segments(x0=l1,x1=length(names(box1[[1]])),y0=pp[3]-yyy,y1=pp[4]+yyy,xpd=TRUE)
  text(l1/2,pp[4]+yyy+strheight(main_box1),main_box1,xpd=TRUE,adj=0.5)
  text(l1+l2/2,pp[4]+yyy+strheight(main_box2),main_box2,xpd=TRUE,adj=0.5)
  #
  yyyy <- c(yyyy,pp[4]+yyy)
  for(i in 1:length(box1)){
    x1 <- box1[[i]]
    x2 <- lapply(x1,get_b)
    x20 <- c(min(do.call(rbind,x2)[,6]),max(do.call(rbind,x2)[,7]))
    y1 <- yyyy[i]; y2 <- yyyy[i+1]
    y21 <- y2-y1
    y1 <- y1+y21*0.05; y2 <- y2-y21*0.15;
    segments(x0=-0.125,x1=0,y0=c(y1,y2),y1=c(y1,y2),col='dark grey',xpd=TRUE)
    text(-0.15,c(y1,y2),format(c(x20[1],x20[2]),digits=3),adj=1,xpd=TRUE,cex=0.5)
    for(j in 1:length(box1[[i]])){
      x2 <- get_b(x1[[j]])
      if(length(x2)>0){
        y3 <- y1+(x2[c(1:5,8)]-x20[1])/(x20[2]-x20[1])*(y2-y1)
        rect(xleft=j-1+0.2,xright=j-0.2,ybottom=y3[2],ytop=y3[4],border='light grey',col='light grey',xpd=TRUE)
        segments(x0=j-0.5,x1=j-0.5,y0=y3[1],y1=y3[5],col='grey',xpd=TRUE)
        segments(x0=j-1+0.2,x1=j-0.2,y0=y3[3],y1=y3[3],col='dark grey',xpd=TRUE) ## median
        points(j-0.5,y3[6],pch=16,xpd=TRUE) ## mean
      }
    }
    text(pv_col_box1,y2+y21*0.1,pval_box1[i],col=ifelse(as.numeric(gsub("\\*","",pval_box1[i]))<=0.1,2,1),cex=0.5,xpd=TRUE,adj=0.5)
    x1 <- box2[[i]]
    x2 <- lapply(x1,get_b)
    x20 <- c(min(do.call(rbind,x2)[,6]),max(do.call(rbind,x2)[,7]))
    y1 <- yyyy[i]; y2 <- yyyy[i+1]
    y21 <- y2-y1
    y1 <- y1+y21*0.05; y2 <- y2-y21*0.15;
    segments(x0=l1+l2+0.125,x1=l1+l2,y0=c(y1,y2),y1=c(y1,y2),col='dark grey',xpd=TRUE)
    text(l1+l2+0.15,c(y1,y2),format(c(x20[1],x20[2]),digits=3),adj=0,xpd=TRUE,cex=0.5)
    for(j in 1:length(box2[[i]])){
      x2 <- get_b(x1[[j]])
      if(length(x2)>0){
        y3 <- y1+(x2[c(1:5,8)]-x20[1])/(x20[2]-x20[1])*(y2-y1)
        rect(xleft=l1+j-1+0.2,xright=l1+j-0.2,ybottom=y3[2],ytop=y3[4],border='light grey',col='light grey',xpd=TRUE)
        segments(x0=l1+j-0.5,x1=l1+j-0.5,y0=y3[1],y1=y3[5],col='grey',xpd=TRUE)
        segments(x0=l1+j-1+0.2,x1=l1+j-0.2,y0=y3[3],y1=y3[3],col='dark grey',xpd=TRUE) ## median
        points(l1+j-0.5,y3[6],pch=16,xpd=TRUE) ## mean
      }
    }
    text(pv_col_box2+l1,y2+y21*0.1,pval_box2[i],col=ifelse(as.numeric(gsub("\\*","",pval_box2[i]))<=0.1,2,1),cex=0.5,xpd=TRUE,adj=0.5)
  }
  text((1:l1)-0.5,pp[3]-yyy-strheight(main_box1)/2,names(box1[[1]]),xpd=TRUE,srt=60,adj=1)
  text(l1+(1:l2)-0.5,pp[3]-yyy-strheight(main_box1)/2,names(box2[[1]]),xpd=TRUE,srt=60,adj=1)
  ###
  return(TRUE)
}

#############################################################################following are functions for inner use##########################################################################
##
RP.main_dir      <- '/Volumes/project_space' ## local path
RP.bash.main_dir <- '/research/projects/yu3grp' ## server path
RP.python3.path  <- '/hpcf/apps/python/install/3.6.1/bin/python3.6' ## required
RP.load          <- 'module load python/3.6.1' ## required
RP.MICA.main     <- '/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/scMINER-master/' ## required for MICA
RP.SJAR.main     <- '/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/SJARACNe-master/' ## required for SJAracne
RP.SJAR.pre_cmd  <- sprintf("export SJARACNE_PATH=%s \nexport PYTHON_PATH=%s \n%s",RP.SJAR.main,RP.python3.path,RP.load) ## required
##
#' Inner use: prepare for MICA
#' @param mat matrix
#' @param outdir character
#' @param prjname character
#' @param all_k a vector of integers
#' @param retransformation character
#' @param perplexity numeric
#' @export
SJ.MICA.prepare <- function(mat,outdir = '.',prjname = NULL,all_k=NULL,retransformation="False",perplexity=30) {
  if(is.null(prjname)){
    message('prjname should be input, please check and re-try !');return(FALSE)
  }
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  if(is.null(all_k)){
    all_k <- 2:12
  }
  input_exp <- file.path(outdir, 'input.exp')
  mat1 <- t(mat)
  mat1 <- cbind(rownames(mat1), mat1)
  colnames(mat1)[1] <- 'id'
  write.table(mat1,file = input_exp,row.names = FALSE,col.names = TRUE,sep = '\t',quote = FALSE)
  message(sprintf('Output exp file into %s',input_exp))
  run_file <- file.path(outdir, 'run_MICA_S1.sh')
  cmd <- sprintf('%s\nsh %s/run_MIE.sh %s %s/ %s',RP.load,RP.MICA.main,prjname,outdir,input_exp)
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  run_file <- file.path(outdir, 'run_MICA_S2.sh')
  cmd <-
    sprintf('%s\nsh %s/run_MICA.sh %s %s/ %s %s %s \"%s\"',RP.load,RP.MICA.main,prjname,
            outdir,input_exp,retransformation,perplexity,paste(all_k,collapse = ' ')
    )
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  return(TRUE)
}


########################## server-run functions (need local run version)
get_server_path <- function(x){
  gsub(RP.main_dir, RP.bash.main_dir, x)
}
#' Inner use: prepare for MICA with file
#' @param input_exp character
#' @param outdir character
#' @param prjname character
#' @param all_k a vector of integers
#' @param retransformation character
#' @param perplexity numeric
#' @export
SJ.MICA.prepare.withfile <- function(input_exp,outdir = '.',prjname = NULL,all_k=NULL,retransformation="False",perplexity=30) {
  if(is.null(prjname)){
    message('prjname should be input, please check and re-try !');return(FALSE)
  }
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  run_file <- file.path(outdir, 'run_MICA_S1.sh')
  cmd <- sprintf('%s\nsh %s/run_MIE.sh %s %s/ %s',RP.load,RP.MICA.main,prjname,outdir,input_exp)
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  run_file1 <- run_file
  run_file <- file.path(outdir, 'run_MICA_S2.sh')
  cmd <- sprintf('%s\nsh %s/run_MICA.sh %s %s/ %s %s %s \"%s\"',RP.load,RP.MICA.main,prjname,outdir,
                 input_exp,retransformation,perplexity,paste(all_k,collapse = ' '))
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  return(paste0('sh ',c(run_file1,run_file)))
}
##
#' Inner use: prepare for SJaracne run on SJ server
#' prepare SJAracne dataset for net-dataset
#' project_name expression_matrix hub_genes outdir
#' @param eset eSet object
#' @param use.samples a vector of characters
#' @param TF_list a vector of characters
#' @param SIG_list a vector of characters
#' @param SJAR.project_name character
#' @param SJAR.main_dir character
#' @param SJAR.bash.main_dir character
#' @param mem integer
#' @param IQR.thre numeric
#' @param IQR.loose_thre numeric
#' @export
SJ.SJAracne.prepare <-
  function(eset,use.samples = rownames(pData(eset)),TF_list=NULL,SIG_list=NULL,
           SJAR.project_name = "",
           SJAR.main_dir = NULL,
           SJAR.bash.main_dir = NULL,
           mem = 40960,IQR.thre=0.5,IQR.loose_thre=0.1) {
    RP.main_dir      <- '/Volumes/project_space'
    RP.bash.main_dir <- '/research/projects/yu3grp'
    RP.python3.path  <- '/home/xdong/bin/python3.6' ## required
    RP.load          <- 'module load python/3.6.1' ## required
    RP.MICA.main     <- '/home/xdong/softs/scMINER-master/' ## required for MICA
    RP.SJAR.main     <- '/home/xdong/softs/SJARACNe-master/' ## required for SJAracne
    RP.SJAR.pre_cmd  <- sprintf("export SJARACNE_PATH=%s \nexport PYTHON_PATH=%s \n%s",RP.SJAR.main,RP.python3.path,RP.load) ## required
    if(is.null(TF_list)==TRUE){
      message('Empty TF_list, please check and re-try !');return(FALSE)
    }
    if(is.null(SIG_list)==TRUE){
      message('Empty SIG_list, please check and re-try !');return(FALSE)
    }
    if(is.null(SJAR.project_name)==TRUE){
      message('Empty SJAR.project_name, please check and re-try !');return(FALSE)
    }
    SJAR.outdir <- file.path(SJAR.main_dir, SJAR.project_name)
    if (!file.exists(SJAR.outdir)) {
      dir.create(SJAR.outdir, recursive = TRUE)
    }
    SJAR.outdir.tf  <-
      file.path(SJAR.main_dir, SJAR.project_name, 'output_tf_')
    SJAR.outdir.sig <-
      file.path(SJAR.main_dir, SJAR.project_name, 'output_sig_')
    SJAR.expression_matrix <-
      file.path(SJAR.main_dir, SJAR.project_name, 'input.exp')
    SJAR.hub_genes.tf <-
      file.path(SJAR.main_dir, SJAR.project_name, 'tf.txt')
    SJAR.hub_genes.sig <-
      file.path(SJAR.main_dir, SJAR.project_name, 'sig.txt')
    SJAR.bash_file.tf <-
      file.path(SJAR.main_dir, SJAR.project_name, 'run_tf.sh')
    SJAR.bash_file.sig <-
      file.path(SJAR.main_dir, SJAR.project_name, 'run_sig.sh')
    ## process exp matrix
    d <- exprs(eset)[, use.samples]
    # filter genes with count=0
    d <- d[!apply(d == 0, 1, all), ]
    # filter genes with IQR
    choose1 <- IQR.filter(d, rownames(d),thre = IQR.thre,loose_thre=IQR.loose_thre,loose_gene=unique(c(TF_list,SIG_list)))
    d <- d[choose1, ]
    use.genes <- rownames(d)
    use.genes <- use.genes[which(is.na(use.genes)==FALSE)]
    w1 <- which(rownames(d) %in% use.genes)
    d <- d[w1,]
    use.genes <- rownames(d)
    # write exp data to exp format
    use.genes.symbol <- use.genes
    if('geneSymbol' %in% colnames(fData(eset))){use.genes.symbol <- fData(eset)[use.genes,'geneSymbol']}
    if('external_gene_name' %in% colnames(fData(eset))){use.genes.symbol <- fData(eset)[use.genes,'external_gene_name']}
    expdata <- data.frame(cbind(isoformId = use.genes, geneSymbol = use.genes.symbol, d))
    #
    write.table(
      expdata,
      file = SJAR.expression_matrix,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    ##
    cat(intersect(use.genes, TF_list),file = SJAR.hub_genes.tf,sep = '\n')
    cat(intersect(use.genes, SIG_list),file = SJAR.hub_genes.sig,sep = '\n')
    # write scripts to bash file for tf
    cmd_tf1 <-
      sprintf(
        '%s %s/generate_pipeline.py %s %s %s %s --run False --host CLUSTER',
        RP.python3.path,
        RP.SJAR.main,
        SJAR.project_name,
        SJAR.expression_matrix,
        SJAR.hub_genes.tf,
        SJAR.outdir.tf
      )
    SJAR.outdir.tf.script <-
      sprintf('%ssjaracne_%s_scripts_/',
              SJAR.outdir.tf,
              SJAR.project_name)
    cmd_tf2 <-
      sprintf(
        'cat %s/02_bootstrap_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.tf.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name
      )
    cmd_tf3 <-
      sprintf(
        'cat %s/03_getconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.tf.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name
      )
    cmd_tf4 <-
      sprintf(
        'cat %s/04_getenhancedconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.tf.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name
      )
    all_cmd <-
      paste(RP.SJAR.pre_cmd,
            cmd_tf1,
            cmd_tf2,
            cmd_tf3,
            cmd_tf4,
            sep = '\n')
    all_cmd <- get_server_path(all_cmd)
    cat(all_cmd, file = SJAR.bash_file.tf, sep = '\n')
    # write scripts to bash file for sig
    cmd_sig1 <-
      sprintf(
        '%s %s/generate_pipeline.py %s %s %s %s --run False --host CLUSTER',
        RP.python3.path,
        RP.SJAR.main,
        SJAR.project_name,
        SJAR.expression_matrix,
        SJAR.hub_genes.sig,
        SJAR.outdir.sig
      )
    SJAR.outdir.sig.script <-
      sprintf('%ssjaracne_%s_scripts_/',
              SJAR.outdir.sig,
              SJAR.project_name)
    cmd_sig2 <-
      sprintf(
        'cat %s/02_bootstrap_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.sig.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name
      )
    cmd_sig3 <-
      sprintf(
        'cat %s/03_getconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.sig.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name
      )
    cmd_sig4 <-
      sprintf(
        'cat %s/04_getenhancedconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.sig.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name
      )
    all_cmd <-
      paste(RP.SJAR.pre_cmd,
            cmd_sig1,
            cmd_sig2,
            cmd_sig3,
            cmd_sig4,
            sep = '\n')
    all_cmd <- get_server_path(all_cmd)
    cat(all_cmd, file = SJAR.bash_file.sig, sep = '\n')
    ##
    message(sprintf('Check %s',SJAR.bash_file.tf))
    message(sprintf('Check %s',SJAR.bash_file.sig))
    return(
      list(
        outdir.tf = SJAR.outdir.tf,
        outdir.sig = SJAR.outdir.sig,
        bash.tf = get_server_path(SJAR.bash_file.tf),
        bash.sig = get_server_path(SJAR.bash_file.sig)
      )
    )
  }

#' Inner use: auto generate bash files for Step1, Step2, Step3, Step4 for all bash files under one directory
#' @param out.dir.SJAR character
#' @export
SJ.SJAracne.step <- function(out.dir.SJAR) {
  tf_bash <- list.files(out.dir.SJAR, pattern = 'tf.sh', recursive = TRUE)
  sig_bash <-
    list.files(out.dir.SJAR, pattern = 'sig.sh', recursive = TRUE)
  s1 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 4 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  s2 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 5 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  s3 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 6 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  s4 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 7 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  cat(
    c(RP.SJAR.pre_cmd, unlist(s1)),
    file = sprintf('%s/step1.bash', out.dir.SJAR),
    sep = '\n'
  )
  cat(
    c(RP.SJAR.pre_cmd, unlist(s2)),
    file = sprintf('%s/step2.bash', out.dir.SJAR),
    sep = '\n'
  )
  cat(
    c(RP.SJAR.pre_cmd, unlist(s3)),
    file = sprintf('%s/step3.bash', out.dir.SJAR),
    sep = '\n'
  )
  cat(
    c(RP.SJAR.pre_cmd, unlist(s4)),
    file = sprintf('%s/step4.bash', out.dir.SJAR),
    sep = '\n'
  )
  message(get_server_path(sprintf('Check: %s', sprintf('%s/step1-4.bash', out.dir.SJAR))))
  return(TRUE)
}
##
#' GSEA (gene set enrichment analysis) plot for the Synergy Inference by data-driven Network-based Bayesian Analysis (SINBA) analysis results.
#'
#' \code{draw.GSEA.NetBID.SINBA} will generate a GSEA plot for Synergy Inference by data-driven Network-based Bayesian Analysis (SINBA)
#' analysis results. SINBA calculates the synergistic effect between a seed driver and a partner driver.
#' The plot includes the GSEA plot for the seed driver and the partner driver independently and
#' the GSEA plot for the combination for the seed driver to each partner driver.
#' The statistics on the plot include the differentiated expression (DE), differentiated activity (DA) for each driver,
#' and the different Z (deltaZ) showing the difference between the combination of the seed and the partner driver to the sum of the original Z statistics.
#'
#' This is a plot function to draw GSEA for synergistic effect prediction between the seed driver and a list of partner drivers.
#' User need to input the differentiated expression information, and choose to display the target genes in one row or two rows,
#' by selecting black color or red to blue color bar.
#'
#' @param DE data.frame,the differentiated expression results.
#' This data.frame could be generated by using \code{getDE.limma.2G} or \code{getDE.BID.2G}.
#' If user want to generate this data.frame by other strategies, the rownames must be the gene names or need one column to be the gene name
#' (set in \code{name_col}) and must contain the columns indicating the differentiated expression profile.
#' @param name_col character, the name of the column in \code{DE}, which contains the gene name. If NULL, will use the rownames of \code{DE}.
#' Default is NULL.
#' @param profile_col character, the name of the column in \code{DE}, which will be used as the differentiated expression profile.
#' If DE is created by \code{getDE.limma.2G} or \code{getDE.BID.2G}, this parameter could be 'logFC' or 't'.
#' @param profile_trend character, the choice of how to display the profile, from high/positive to low/negative ('pos2neg')
#' or low/negative to high/positive ('neg2pos').Default is 'pos2neg'.
#' @param seed_driver character, name for the seed driver.
#' @param partner_driver_list a vector of characters, name for the partner driver list.
#' @param seed_driver_label character, label for the seed driver displayed on the plot. Default is seed_driver.
#' @param partner_driver_label a vector of characters, label for the partner driver list displayed on the plot. Default is partner_driver_list
#' @param driver_DA_Z a vector of numeric values, the Z statistics of differentiated activity (DA) for the driver list.
#' Better to give name to the vector, otherwise will automatically use driver list (seed + partner) as the name.
#' @param driver_DE_Z a vector of numeric values, the Z statistics of differentiated expression (DE) for the driver list.
#' Better to give name to the vector, otherwise will automatically use driver list (seed + partner) as the name.
#' @param target_list a list for the target gene information for the drivers. The names for the list must contain the driver in driver list (seed + partner)
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param DA_Z_merge a vector of numeric values, the Z statistics of differentiated activity (DA) for the combination of the seed driver to partner drivers saperately.
#' Better to give name to the vector, otherwise will automatically use partner driver list as the name.
#' @param target_list_merge a list for the target gene information for thecombination of the seed driver to partner drivers saperately.
#' The names for the list must contain the driver in partner_driver_list
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list_merge} could be automatically generated by \code{merge_target_list}.
#' @param top_driver_number numeric, number for the top significant partner drivers to be displayed on the plot. Default is 10.
#' @param top_order character, choice of order pattern used to display the partner drivers. Two options,'merge' or 'diff'.
#' 'merge' means the partner drivers will be sorted by the combined Z statistics.
#' 'diff' means the partner drivers will be sorted by the delta Z statistics.
#' Default is 'merge'.
#' @param target_nrow numeric, number of rows for each driver display on the plot. Two options, 1 or 2.
#' If set to 1, the target genes' position on the profile will be displayed in one row.
#' If set to 2, the target genes' position on the profile will be displayed in two rows,
#' with positive regulated genes displayed on the first row and negative regulated genes displayed on the second row.
#' Default is 2.
#' @param target_col character, choice of color pattern used to display the targets. Two options,'black' or 'RdBu'.
#' If set to 'black', the lines will be colored in black.
#' If set to 'RdBu', the lines will be colored into Red to Blue color bar.
#' If \code{target_col_type} is set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If \code{target_col_type} is set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' with significant high set for red and low for blue. The significant threshold is set by \code{profile_sig_thre}.
#' Default is 'RdBu'.
#' @param target_col_type character, choice of the pattern used to display the color for target genes, only work when \code{target_col} is set as 'RdBu'.
#' Two options,'PN' or 'DE'.
#' If set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' Default is 'PN'.
#' @param left_annotation character, annotation displayed on the left of the figure representing left condition of the rank_profile. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure representing right condition of the rank_profile. Default is "".
#' @param main character, title for the plot. Default is "".
#' @param profile_sig_thre numeric, threshold for the absolute values in profile to be treated as significance.
#' Target genes without signifcant values in the profile will be colored in grey. Only work when \code{target_col_type} is set as "DE" and \code{target_col} is set as "RdBu".
#' Default is 0.
#' @param Z_sig_thre numeric, threshold for the Z statistics in \code{driver_DA_Z} and \code{driver_DE_Z} to be treated as signifcance.
#' Only signifcant values will have background color. Default is 1.64.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return logical value indicating whether the plot has been successfully generated

#' @examples
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' ## choose seed driver and partner driver list
#' seed_driver <- driver_list[1]
#' part_driver <- ms_tab$originalID_label
#' ## get merge target
#' merge_target <- lapply(part_driver,function(x){
#'   m1 <- merge_target_list(driver1=seed_driver,driver2=x,
#'                           target_list=analysis.par$merge.network$target_list)
#' })
#' names(merge_target) <- part_driver
#' ## get activity matrix for the merge target network
#' ac_combine_mat <- cal.Activity(all_target=merge_target,
#'                                cal_mat=exprs(analysis.par$cal.eset),
#'                                es.method='weightedmean')
#' ## get DA for the combined drivers
#' comp_name <- 'G4.Vs.others'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')]
#' # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')]
#' # get sample list for G1
#' DA_driver_combine <- getDE.limma.2G(eset=generate.eset(ac_combine_mat),
#'                                     G1=G1,G0=G0,
#'                                     G1_name='G4',G0_name='others')
#' ## or use: DA_driver_combine <- getDE.BID.2G(eset=generate.eset(ac_combine_mat),
#'                                     G1=G1,G0=G0,
#'                                     G1_name='G4',G0_name='others')
#' ## prepare for SINBA input
#' ori_part_Z <- analysis.par$DA[[comp_name]][part_driver,'Z-statistics']
#' ori_seed_Z <- analysis.par$DA[[comp_name]][seed_driver,'Z-statistics']
#' DE <- analysis.par$DE[[comp_name]]
#' driver_DA_Z <- analysis.par$DA[[comp_name]][,'Z-statistics']
#' names(driver_DA_Z) <- rownames(analysis.par$DA[[comp_name]])
#' driver_DE_Z <- analysis.par$DE[[comp_name]][,'Z-statistics']
#' names(driver_DE_Z) <- rownames(analysis.par$DE[[comp_name]])
#' DA_Z_merge <- DA_driver_combine[,'Z-statistics']
#' names(DA_Z_merge) <- rownames(DA_driver_combine)
#' target_list_merge <- merge_target
#' seed_driver_label <- ms_tab[seed_driver,'gene_label']
#' partner_driver_list <- part_driver
#' profile_col <- 't'
#' partner_driver_label <- ms_tab[partner_driver_list,'gene_label']
#' target_list <- analysis.par$merge.network$target_list
##
#' draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,
#'                        seed_driver=seed_driver,
#'                        partner_driver_list=partner_driver_list,
#'                        seed_driver_label=seed_driver_label,
#'                        partner_driver_label=partner_driver_label,
#'                        driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,
#'                        target_list=target_list,
#'                        DA_Z_merge=DA_Z_merge,
#'                        target_list_merge=target_list_merge,
#'                        top_driver_number=20,profile_trend='pos2neg',
#'                        top_order='merge',Z_sig_thre = 1.64,
#'                        target_nrow=1,target_col='RdBu',target_col_type='PN',
#'                        pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo1.pdf',
#'                        analysis.par$out.dir.PLOT))
#'}
#' @export
draw.GSEA.NetBID.SINBA <- function(DE=NULL,name_col=NULL,profile_col=NULL,profile_trend='pos2neg',
                                   seed_driver=NULL,partner_driver_list=NULL,
                                   seed_driver_label=seed_driver,partner_driver_label=partner_driver_list,
                                   driver_DA_Z=NULL,driver_DE_Z=NULL,target_list=NULL,
                                   DA_Z_merge=NULL,target_list_merge=NULL,
                                   top_driver_number=10,top_order='merge',target_nrow=2,target_col='RdBu',target_col_type='PN',
                                   left_annotation="",right_annotation="",main="",
                                   profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  if(!profile_col %in% colnames(DE)){
    message(sprintf('%s not in colnames of DE, please check and re-try!',profile_col))
    return(FALSE)
  }
  driver_list <- c(seed_driver,partner_driver_list)
  driver_list_gene <- gsub('(.*)_.*','\\1',driver_list)
  show_label <- c(seed_driver_label,partner_driver_label)
  # get names
  if(is.null(names(driver_DA_Z))) names(driver_DA_Z) <- driver_list
  if(is.null(names(driver_DE_Z))) names(driver_DE_Z) <- driver_list
  if(is.null(names(show_label))) names(show_label) <- driver_list
  if(is.null(names(DA_Z_merge))) names(DA_Z_merge) <- driver_list
  #
  ori_part_Z <- driver_DA_Z[partner_driver_list]
  ori_seed_Z <- driver_DA_Z[seed_driver]
  diff_Z <- 2*DA_Z_merge[partner_driver_list]-(ori_part_Z+ori_seed_Z)
  names(diff_Z) <- partner_driver_list
  driver_DA_Z <- driver_DA_Z[driver_list]
  driver_DE_Z <- driver_DE_Z[driver_list_gene]; names(driver_DE_Z) <- driver_list
  DA_Z_merge  <- DA_Z_merge[partner_driver_list]
  #
  if(top_order=='merge'){
    if(length(partner_driver_list)>top_driver_number){
      #partner_driver_list <- partner_driver_list[order(abs(DA_Z_merge[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)][1:top_driver_number] ## only consider positive part
    }
    if(profile_trend=='pos2neg')
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = FALSE)]
    else
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)]
  }else{
    if(length(partner_driver_list)>top_driver_number){
      #partner_driver_list <- partner_driver_list[order(abs(diff_Z[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = TRUE)][1:top_driver_number] ## only consider positive increase
    }
    if(profile_trend=='pos2neg')
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = FALSE)]
    else
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = TRUE)]
  }
  #
  driver_list <- c(seed_driver,partner_driver_list)
  show_label <- show_label[driver_list]
  driver_DA_Z <- driver_DA_Z[driver_list]
  driver_DE_Z <- driver_DE_Z[driver_list]
  diff_Z <- diff_Z[partner_driver_list]
  DA_Z_merge <- DA_Z_merge[partner_driver_list]
  if(is.null(name_col)==TRUE){
    DE <- cbind(DE[,setdiff(colnames(DE),'name')],name=rownames(DE),stringsAsFactors=FALSE)
    name_col <- 'name'
  }
  w1 <- which(is.na(DE[,profile_col])==FALSE)
  DE <- DE[w1,]
  if(profile_trend=='pos2neg') DE <- DE[order(DE[,profile_col],decreasing = TRUE),] else DE <- DE[order(DE[,profile_col],decreasing = FALSE),]
  DE_profile <- DE[,profile_col]
  DE_profile_name <- DE[,name_col]
  ##############################################
  ## calculate layout

  target_list <- lapply(target_list,function(x)x[which(x$target %in% DE_profile_name),])
  target_list_merge <- lapply(target_list_merge,function(x)x[which(x$target %in% DE_profile_name),])

  n_gene <- length(DE_profile)
  if(target_nrow==2){
    n_driver <- length(partner_driver_list)*4+2
    ratio1 <- ceiling(n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  } else {
    n_driver <- length(partner_driver_list)*2+1
    ratio1 <- ceiling(1.5*n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  }
  #
  if(is.null(pdf_file)==FALSE){
    pdf(pdf_file,width=(rr*2+ratio2)*2,height=(ratio1+rr)*2)
  }
  # get layout
  layout(matrix(c(rep(0,length.out=rr),rep(1,length.out=ratio2),rep(0,length.out=rr*1),
                  rep(c(rep(4,length.out=rr),rep(2,length.out=ratio2),rep(3,length.out=rr*1)),
                      length.out=ratio1*(ratio2+rr*2))),
                ncol=c(ratio2+rr*2),byrow=TRUE))
  ## plot 1
  par(mar=c(1.5,1.5,1.5,0))
  mm <- quantile(DE_profile,probs=c(0.0001,0.9999));
  mm <- max(abs(mm)); mm <- c(-mm,mm)
  y1 <- seq(mm[1],mm[2],length.out=5); y1 <- round(y1,1)
  unit <- n_gene/10; unit <- round(unit/100)*100
  x1 <- seq(0,n_gene,by=unit);x1 <- unique(x1); x1 <- c(x1,max(x1)+unit)
  par(usr=c(0,length(DE_profile),mm[1],mm[2]))
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
  pp <- par()$usr; rr <- (pp[2]-pp[1])/n_gene
  polygon(x=c(pp[1],c(1:n_gene)*rr+pp[1],pp[2]),y=c(0,DE_profile,0),col='grey',border='grey',xpd=TRUE,lwd=0.3)
  if(profile_trend=='pos2neg'){
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[2]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[1]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }else{
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[1]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[2]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }
  axis(side=2,at=y1,labels=y1)
  mtext(side=2,line = 2.5,profile_col,cex=1)
  segments(pp[1],mm[1],pp[2],mm[1],xpd=TRUE)
  segments(x1*rr,mm[1]-(mm[2]-mm[1])/50,x1*rr,mm[1],xpd=TRUE)
  text(x1*rr,mm[1]-(mm[2]-mm[1])/25,get_label_manual(x1),adj=0.5,xpd=TRUE)
  ## plot2
  par(mar=c(2,1.5,2,0))
  plot(1,col='white',xlab="",ylab="",xlim=c(0,n_gene),xaxt='n',yaxt='n')

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow

  rect(xleft = pp[1],xright=pp[2],ybottom = yy1[length(yy1)-target_nrow],ytop=yy1[length(yy1)],border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.3)) ## for seed rows
  rect(xleft = pp[1],xright=pp[2],ybottom = yy2,ytop=yy4,border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.2)) ## for combine rows

  segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
  segments(x0=pp[1],x1=pp[2],y0=yy3,y1=yy3,lwd=1,col='dark grey')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1,col='black')
  segments(x0=pp[1],x1=pp[2],y0=yy1[length(yy1)-target_nrow],y1=yy1[length(yy1)-target_nrow],lwd=1.5,col='black')

  # shorten yy1
  dyy <- yy1[2]-yy1[1]
  yy11 <- yy1-dyy*0.3
  # add columns
  use_target_list <- target_list[driver_list]
  use_merge_target_list <- target_list_merge[partner_driver_list]

  if(target_col_type=='DE'){
    cc <- z2col(DE_profile,sig_thre=profile_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[5:9],blue_col=brewer.pal(9,'Blues')[5:9],
                col_max_thre=max(abs(DE_profile)))
    #names(cc) <- names(DE_profile)
    cc[which(cc=='white')] <- 'light grey'
  }
  if(target_nrow==1){
    # for seed driver
    t1 <- use_target_list[[seed_driver]]
    w0 <- which(DE_profile_name %in% t1$target)
    w1 <- w0*rr+pp[1]
    if(target_col=='black'){
      segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
    }else{
      if(target_col_type=='DE'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=cc[w0])
      }else{
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      w0 <- which(DE_profile_name %in% t1$target)
      w1 <- w0*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
        }
      }
    }
    # for combine
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      w0 <- which(DE_profile_name %in% t1$target)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(DE_profile_name %in% t_over)
      w1_over <- w0_over*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
        }
      }
      points(w1_over,rep((yy11[2*i]+yy1[2*i])/2,length.out=length(w1_over)),pch='*',col='black')
    }
  }
  ###################
  if(target_nrow==2){
    # for seed driver
    t1 <- use_target_list[[seed_driver]]
    t11 <- t1[which(t1$spearman>=0),]$target
    t12 <- t1[which(t1$spearman<0),]$target
    w0 <- which(DE_profile_name %in% t11)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=pos_col,lwd=1.5)
        }
      }
    }
    w0 <- which(DE_profile_name %in% t12)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=neg_col,lwd=1.5)
        }
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(DE_profile_name %in% t11)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=pos_col,lwd=1.5)
          }
        }
      }
      w0 <- which(DE_profile_name %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=neg_col,lwd=1.5)
          }
        }
      }
    }
    # for each partner driver + seed combination
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(DE_profile_name %in% t11)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(DE_profile_name %in% t_over) ## setdiff(t_over,names(DE_profile)) !!!
      w1_over <- w0_over*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=pos_col,lwd=1.5)
          }
        }
      }
      points(intersect(w1_over,w1),rep((yy11[4*i-1]+yy1[4*i-1])/2,length.out=length(intersect(w1_over,w1))),pch='*',col='black')
      w0 <- which(DE_profile_name %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=neg_col,lwd=1.5)
          }
        }
      }
      points(intersect(w1_over,w1),rep((yy11[4*i-2]+yy1[4*i-2])/2,length.out=length(intersect(w1_over,w1))),pch='*',col='black')
    }
    ####
  }
  ###################
  ## plot 3
  par(mar=c(2,0.5,2,2))
  plot(1,col='white',xlab="",ylab="",xlim=c(0,3),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow
  xx1 <- seq(pp[1],pp[2],length.out=4)

  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white',xpd=TRUE)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey',xpd=TRUE)
  abline(v=xx1)

  ## add text
  mm_min <- min(min(abs(driver_DA_Z[partner_driver_list]),na.rm=TRUE)*0.9,min(abs(driver_DE_Z[partner_driver_list]),na.rm=TRUE)*0.9,
                min(abs(diff_Z[partner_driver_list]),na.rm=TRUE)*0.9,min(abs(DA_Z_merge[partner_driver_list]),na.rm=TRUE)*0.9)
  mm_min <- max(mm_min,Z_sig_thre)
  mm_max <- max(max(abs(driver_DA_Z[partner_driver_list]),na.rm=TRUE)*1.1,max(abs(driver_DE_Z[partner_driver_list]),na.rm=TRUE)*1.1,
                max(abs(diff_Z[partner_driver_list]),na.rm=TRUE)*1.1,max(abs(DA_Z_merge[partner_driver_list]),na.rm=TRUE)*1.1)
  c1 <- z2col(driver_DA_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c2 <- z2col(driver_DE_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c3 <- z2col(diff_Z[partner_driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c4 <- z2col(DA_Z_merge[partner_driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)

  # for seed driver
  yy1<-yy3
  z1 <- driver_DA_Z[seed_driver]
  z2 <- driver_DE_Z[seed_driver]
  p1 <- get_z2p(z1)
  p2 <- get_z2p(z2)
  rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[length(yy1)-1],ytop=yy1[length(yy1)],col=c1[1],border='dark grey',xpd=TRUE)
  rect(xright=xx1[3],xleft=xx1[2],ybottom=yy1[length(yy1)-1],ytop=yy1[length(yy1)],col=c2[1],border='dark grey',xpd=TRUE)
  text(x=sum(xx1[1:2])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,p1,adj=0.5)
  text(x=sum(xx1[3:4])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,p2,adj=0.5)
  text(x=sum(xx1[2:3])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,'-',adj=0.5)

  # for partner driver
  for(i in 1:length(partner_driver_list)){
    z1 <- driver_DA_Z[partner_driver_list[i]]
    z2 <- driver_DE_Z[partner_driver_list[i]]
    z3 <- diff_Z[partner_driver_list[i]]
    z4 <- DA_Z_merge[partner_driver_list[i]]
    p1 <- get_z2p(z1)
    p2 <- get_z2p(z2)
    p3 <- format(z3,digits=3)
    p4 <- get_z2p(z4)
    rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[2*i],ytop=yy1[2*i+1],col=c1[i+1],border='dark grey',xpd=TRUE) ## DA_Z
    rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[2*i-1],ytop=yy1[2*i],col=c4[i],border='dark grey',xpd=TRUE) ## merge_Z
    rect(xleft=xx1[3],xright=xx1[4],ybottom=yy2[i],ytop=yy2[i+1],col=c2[i+1],border='dark grey',xpd=TRUE) ## DE
    rect(xleft=xx1[2],xright=xx1[3],ybottom=yy2[i],ytop=yy2[i+1],col=c3[i],border='dark grey',xpd=TRUE) ## delta Z
    text(x=sum(xx1[1:2])/2,y=(yy1[2*i]+yy1[2*i+1])/2,p1,adj=0.5) ## DA_Z
    text(x=sum(xx1[1:2])/2,y=(yy1[2*i-1]+yy1[2*i])/2,p4,adj=0.5) ## merge_Z
    text(x=sum(xx1[3:4])/2,y=(yy2[i]+yy2[i+1])/2,p2,adj=0.5) ## DE
    text(x=sum(xx1[2:3])/2,y=(yy2[i]+yy2[i+1])/2,p3,adj=0.5) ## delta z
  }

  textheight <- strheight('DA',units='user',cex=1.5)
  text(sum(xx1[1:2])/2,pp[4]+textheight,'DA',xpd=TRUE,cex=1.5)
  textheight <- strheight('DE',units='user',cex=1.5)
  text(sum(xx1[3:4])/2,pp[4]+textheight,'DE',xpd=TRUE,cex=1.5)
  textheight <- strheight('deltaZ',units='user',cex=1.5)
  text(sum(xx1[2:3])/2,pp[4]+textheight,'deltaZ',xpd=TRUE,cex=1.5)

  ## plot 4
  par(mar=c(2,6,2,0.2))
  plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow
  xx1 <- seq(pp[1],pp[2],length.out=4)

  yy22 <- (yy2[1:(length(yy2)-1)]+yy2[2:length(yy2)])/2
  dyy22 <- yy22[2]-yy22[1]

  yy33 <- (yy3[1:(length(yy3)-1)]+yy3[2:length(yy3)])/2
  dyy33 <- yy33[2]-yy33[1]

  xleft <- pp[1]+(pp[2]-pp[1])*0.55
  tt <- pp[2]-xleft

  text(show_label[1],x=pp[1]+(pp[1]+pp[2])*0.53,y=(yy3[length(yy3)-1]+yy3[length(yy3)])/2,xpd=TRUE,adj=1,cex=1.2) ## label for seed
  text(show_label[2:length(show_label)],x=pp[1]+(pp[1]+pp[2])*0.52,y=yy22,xpd=TRUE,adj=1) ## label for partner

  # add target size
  use_target_list <- target_list[driver_list]
  use_merge_target_list <- target_list_merge[partner_driver_list]
  target_size <- do.call(rbind,lapply(use_target_list,function(x){
    x1 <- length(which(x$spearman>=0))
    x2 <- length(which(x$spearman<0))
    c(x1,x2)
  }))
  merge_target_size <- do.call(rbind,lapply(use_merge_target_list,function(x){
    x1 <- length(which(x$spearman>=0))
    x2 <- length(which(x$spearman<0))
    c(x1,x2)
  }))
  # for seed driver
  if(target_nrow==2){
    mm <- max(merge_target_size)
    i <- length(yy1)-1
    rect(xleft=xleft,xright=xleft+target_size[1,1]/mm*tt,
         ybottom=yy1[i],ytop=yy1[i]+dyy22/2*0.35,col=pos_col,border=NA)
    rect(xleft=xleft,xright=xleft+target_size[1,2]/mm*tt,
         ytop=yy1[i],ybottom=yy1[i]-dyy22/2*0.35,col=neg_col,border=NA)
  }else{
    target_size <- rowSums(target_size)
    merge_target_size <- rowSums(merge_target_size)
    mm <- max(merge_target_size)
    i <- length(yy1)
    rect(xleft=xleft,xright=xleft+target_size[1]/mm*tt,
         ybottom=(yy1[i]+yy1[i-1])/2-dyy22*0.2/2,ytop=(yy1[i]+yy1[i-1])/2+dyy22*0.2/2,col='dark grey',border=NA)
  }
  # for partner
  if(target_nrow==2){
    mm <- max(merge_target_size)
    for(i in 1:length(partner_driver_list)){
      rect(xleft=xleft,xright=xleft+target_size[i+1,1]/mm*tt,
           ybottom=yy33[2*i],ytop=yy33[2*i]+dyy33*0.35,col=pos_col,border=NA)
      rect(xleft=xleft,xright=xleft+target_size[i+1,2]/mm*tt,
           ytop=yy33[2*i],ybottom=yy33[2*i]-dyy33*0.35,col=neg_col,border=NA)
      # merge
      rect(xleft=xleft,xright=xleft+merge_target_size[i,1]/mm*tt,
           ybottom=yy33[2*i-1],ytop=yy33[2*i-1]+dyy33*0.35,col=pos_col,border=NA)
      rect(xleft=xleft,xright=xleft+merge_target_size[i,2]/mm*tt,
           ytop=yy33[2*i-1],ybottom=yy33[2*i-1]-dyy33*0.35,col=neg_col,border=NA)
    }
    segments(x0=xleft,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+xleft
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  #
  if(target_nrow==1){
    #target_size <- rowSums(target_size)
    #merge_target_size <- rowSums(merge_target_size)
    mm <- max(merge_target_size)
    for(i in 1:length(partner_driver_list)){
      rect(xleft=xleft,xright=xleft+target_size[i+1]/mm*tt,ybottom=yy33[i*2]-dyy33*0.2,ytop=yy33[i*2]+dyy33*0.2,col='dark grey',border=NA) #each
      rect(xleft=xleft,xright=xleft+merge_target_size[i]/mm*tt,ybottom=yy33[i*2-1]-dyy33*0.2,ytop=yy33[i*2-1]+dyy33*0.2,col='dark grey',border=NA) #merge
    }
    segments(x0=xleft,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+xleft
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  ## add lines
  segments(x0=xleft,x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='grey',xpd=TRUE)
  abline(v=xleft,col='grey')
  #abline(v=pp[2],col='grey')
  ## test for significant overlap
  total_possible_target <- unique(unlist(lapply(target_list,function(x)x$target)))
  for(i in 1:length(partner_driver_list)){
    res1 <- test.targetNet.overlap(seed_driver,partner_driver_list[i],
                                   target1=intersect(use_target_list[[seed_driver]]$target,DE_profile_name),
                                   target2=intersect(use_target_list[[partner_driver_list[i]]]$target,DE_profile_name),
                                   total_possible_target = total_possible_target)
    pv <- format(res1[1],digits=2,scientific=TRUE)
    ov <- round(res1[3])
    if(res1[1]<0.05){
      text(sprintf('Overlap:%d, P value:%s',ov,pv),x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.8,col='dark red')
    }else{
      if(ov==0){
        text('No overlap',x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.7)
      }else{
        text(sprintf('Overlap:%d, P value:%s',ov,pv),x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.7)
      }
    }
  }
  ##
  if(is.null(pdf_file)==FALSE) dev.off()
  layout(1);
  return(TRUE)
}

##
##
#' find.gsByGene
#' @param gene character
#' @param use_gs a vector of characters
#' @param return_type character, all, length, name
#' @export
find.gsByGene <- function(gene=NULL,use_gs=NULL,return_type='name',use_per=1){
  if(is.null(use_gs)==TRUE){
    all_gg <- unlist(all_gs2gene,recursive = FALSE)
  } else {
    if('all' %in% use_gs){
      all_gg <- unlist(all_gs2gene,recursive = FALSE)
    }else{
      all_gg <-  unlist(all_gs2gene[use_gs],recursive = FALSE)
    }
  }
  x1 <- unlist(lapply(all_gg,function(x){
    length(intersect(x,gene))
  }))
  x2 <- names(x1[which(x1>=use_per*length(gene))])
  x3 <- sort(unlist(lapply(all_gg[x2],function(x)length(x))))
  #print(sort(x3))
  if(return_type=='name') return(x2)
  if(return_type=='length') return(x3)
  if(return_type=='all') return(all_gg[names(x3)])
}
##
par.unit <- function(){
  user.range <- par("usr")[c(2,4)] - par("usr")[c(1,3)]
  region.pct <- par("plt")[c(2,4)] - par("plt")[c(1,3)]
  region.px <-
  dev.size(units="px") * region.pct
  px.per.xy <- region.px / user.range
  return(px.per.xy)
}
##
##
reOrder.sample <- function(obs_class_label,sample_class_order){
  w1 <- order(factor(obs_class_label,levels=sample_class_order))
}
reOrder.eset.sample <- function(eset,use_sample_col,sample_class_order){
  obs_class_label <- get_obs_label(pData(eset),use_sample_col)
  w1 <- order(factor(obs_class_label,levels=sample_class_order))
  generate.eset(exp_mat=exprs(eset)[,w1],pData(eset)[w1,],fData(eset))
}
##
library(easyPubMed)
search_pub <- function(q1){
  p1 <- list()
  for(x in q1){
    print(x);Sys.sleep(0.2);msg <- try(p1[[x]] <- get_pubmed_ids(x))
  }
  c1 <- unlist(lapply(p1,function(x)x$Count));c1 <- as.numeric(c1);
  names(c1) <- names(p1)
  return(list(p1,c1))
}
mat2cluster <- function(mat1,h_gene=0.9,h_gs=0.9,z_thre=-qnorm(0.01),gene_list=NULL,gs_list=NULL,strategy='max'){
  #mat1 <- t(mat1)
  h1_gs <- hclust(dist(mat1,method='binary'))
  h1_gene <- hclust(dist(t(mat1),method='binary'))
  gs_cluster <- cutree(h1_gs,h=h_gs)
  gene_cluster <- cutree(h1_gene,h=h_gene)
  if(is.null(gs_list)==TRUE){gs_list <- vec2list(gs_cluster);  names(gs_list) <- paste0('gs',names(gs_list))}
  if(is.null(gene_list)==TRUE){gene_list <- vec2list(gene_cluster);  names(gene_list) <- paste0('gene',names(gene_list))}
  list_net <- c()
  for(i in 1:length(gs_list)){
    for(j in 1:length(gene_list)){
      if(strategy=='max') v1 <- max(mat1[gs_list[[i]],intersect(gene_list[[j]],colnames(mat1))])
      if(strategy=='mean') v1 <- mean(mat1[gs_list[[i]],intersect(gene_list[[j]],colnames(mat1))])
      list_net <- rbind(list_net,c(names(gs_list)[i],names(gene_list)[j],v1))
    }
  }
  list_net <- list_net[which(list_net[,3]>=z_thre),]
  colnames(list_net) <- c('gs_list','gene_list','maxZ')
  return(list(gs_list=gs_list,gene_list=gene_list,list_net=list_net))
}

##### here arrow_direction is a vector
draw.targetNet.Direction <- function(source_label="",source_z=NULL,edge_score=NULL,
                           label_cex=0.7,source_cex=1,
                           pdf_file=NULL,arrow_direction=NULL,n_layer=1,alphabetical_order=FALSE){
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  edge_score<- sort(edge_score)
  tmp1 <- sapply(unique(names(edge_score)),function(x){
    x1 <- edge_score[which(names(edge_score)==x)]
    x1[which.max(abs(x1))]
  })
  names(tmp1) <- unique(names(edge_score))
  edge_score <- tmp1
  edge_score<- edge_score[order(edge_score,decreasing = TRUE)]
  g1 <- names(edge_score)
  print(edge_score)
  ec <- z2col(edge_score*100,sig_thre=0,n_len=length(edge_score),red_col=pos_col,blue_col=neg_col);names(ec) <- names(edge_score)
  ec <- get_transparent(ec,alpha=0.8); names(ec) <- names(edge_score)
  ew <- 2*label_cex*(abs(edge_score)-min(abs(edge_score)))/(max(abs(edge_score))-min(abs(edge_score)))+label_cex/2; names(ew) <- names(edge_score)
  t2xy <- function(tt,radius=1) {
    t2p <- pi*2 * tt
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  if(alphabetical_order==TRUE)  g1 <- sort(g1)
  plot_part <- function(ori=FALSE,before_off=FALSE){
    geneWidth <- max(strwidthMod(g1,'inches',cex=label_cex))
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=6+2*geneWidth*n_layer,height=6+2*geneWidth*n_layer)
    par(mai=c(1,1,1,1))
    plot(1,xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white',xlab="",ylab="",xaxt='n',yaxt='n')
    pp <- par()$usr
    ## add target label
    #tt <- seq(-0.5,0.5,length.out=length(g1)+1)[-1]; ## g1
    rad_v <- seq(from=0.5,to=1,length.out=n_layer)
    if(n_layer==1) rad_v=0.8
    #tmp1 <- ceiling(length(g1)/sum(rad_v)*rad_v)
    uu <- ceiling(length(g1)/length(rad_v))
    tmp1 <- rep(uu,length.out=length(rad_v))
    if(n_layer>1) tmp1[length(tmp1)] <- length(g1)-sum(tmp1[1:(length(tmp1)-1)]) else tmp1 <- length(g1)
    tmp1<-cumsum(tmp1)
    all_g1 <- list()
    for(i in 1:n_layer){
      if(i==1) all_g1[[i]] <- g1[1:tmp1[i]] else all_g1[[i]]  <- g1[(tmp1[i-1]+1):tmp1[i]]
    }
    if(alphabetical_order==FALSE) all_g1 <- lapply(all_g1,function(x)x[order(edge_score[x])])
    #all_tt <- lapply(1:n_layer,function(i)seq(-0.5,0.5,length.out=length(all_g1[[i]])+1)[-1])
    all_tt <- lapply(1:n_layer,function(i){
      if(i==1) return(seq(-0.5,0.5,length.out=uu+1)[-1][1:length(all_g1[[i]])])
      if(i>1) return(seq(-0.5-(i-1)/(n_layer*uu),0.5-(i-1)/(n_layer*uu),length.out=uu+1)[-1][1:length(all_g1[[i]])])
    })
    all_p <- lapply(1:n_layer,function(i)t2xy(all_tt[[i]],radius=rad_v[i]))
    # add line
    for(i in 1:n_layer){
      each_v <- rad_v[i]
      p1 <- all_p[[i]]
      g1_use <- all_g1[[i]]
      tt <- all_tt[[i]]
      geneWidth <- strwidthMod(source_label,'user',cex=source_cex)
      w1 <- which(arrow_direction[g1_use]==1) ## out
      p2<-t2xy(tt,radius=each_v-label_cex/36);
      p3<-t2xy(tt,radius=each_v-label_cex/48);
      arrows(x0=0,y0=0,x1=p2$x[w1],y1=p2$y[w1],col=ec[g1_use][w1],lwd=ew[g1_use][w1],angle=15,length=0.2*label_cex,xpd=TRUE);
      w1 <- which(arrow_direction[g1_use]== -1) ## in
      p2<-t2xy(tt,radius=each_v-label_cex/36);
      p3<-t2xy(tt,radius=each_v-label_cex/36);
      p4<-t2xy(tt,radius=geneWidth/2);
      arrows(x0=p2$x[w1],y0=p2$y[w1],x1=p4$x[w1],y1=p4$y[w1],col=ec[g1_use][w1],lwd=ew[g1_use][w1],angle=15,length=0.2*label_cex,xpd=TRUE);
      points(p3$x,p3$y,pch=16,col='dark grey',cex=label_cex)
    }
    # add label
    for(j in 1:n_layer){
      p1 <- all_p[[j]]
      g1_use <- all_g1[[j]]
      for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g1_use[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    }
    ## add source label
    if(is.null(source_z)==TRUE){
      draw.ellipse(0,0,a=1.05*geneWidth/2,b=1.05*geneWidth/2,col='light grey',border=NA)
    }else{
      draw.ellipse(0,0,a=1.05*geneWidth/2,b=1.05*geneWidth/2,col=z2col(source_z),border=NA)
    }
    text(0,0,source_label,adj=0.5,xpd=TRUE,cex=source_cex)
    legend(x=pp[1],y=pp[3],fill=c(pos_col,neg_col),c('Positively-regulated','Negatively-regulated'),bty='n',xpd=T,border=NA,cex=label_cex,horiz = TRUE)
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);dev.off();dev.off()} else {plot_part()}
  return(TRUE)
}
##
draw.targetNet.forDriver <- function(use_driver,use_driver_label=use_driver,network_dat=NULL,transfer_tab=NULL,use_label=TRUE,...){
  n1 <- which(gsub('(.*)_.*','\\1',network_dat$source) == use_driver)
  n2 <- which(network_dat$target == use_driver)
  use_net <- network_dat[c(n1,n2),]
  use_net$source.symbol <- get_name_transfertab(use_net$source.symbol,transfer_tab)
  use_net$target.symbol <- get_name_transfertab(use_net$target.symbol,transfer_tab)
  nn1 <- length(n1);nn2 <- length(n2)
  x1 <- network_dat[c(n1,n2),]$MI*sign(network_dat[c(n1,n2),]$spearman)
  x2 <- c(rep('source',length.out=length(n1)),rep('target',length.out=length(n2)))
  x3 <- c(network_dat[c(n1),'target'],network_dat[c(n2),'source'])
  x4 <- get_name_transfertab(gsub("(.*)_.*","\\1",x3),transfer_tab)
  if(use_label==TRUE) x4[(nn1+1):(nn1+nn2)] <-paste0(x4[(nn1+1):(nn1+nn2)],'_',gsub("(.*)_(.*)","\\2",x3[(nn1+1):(nn1+nn2)]))
  names(x1) <- x4;names(x2) <- x4
  draw.targetNet.Direction(use_driver_label,edge_score=x1,arrow_direction = ifelse(x2=='source',1,-1),...)
  return(use_net)
}
##
get.driverTarget.funcEnrich <- function(driver_list=NULL,target_list=NULL,
                                        gs2gene=NULL,use_gs=NULL,bg_list=NULL,
                                        transfer_tab=NULL,min_gs_size=50,max_gs_size=2000,Pv_thre=0.1,Pv_adj = 'none'){
  target_gene <- lapply(driver_list,function(x){
    x1 <- target_list[[x]]$target
    x1 <- x1[which(x1 %in% transfer_tab[,1])]
    x2 <- transfer_tab[which(transfer_tab[,1] %in% x1),]
    target <- unique(x2[,2])
    return(target)
  })
  names(target_gene) <- driver_list
  ##
  f_res <- lapply(target_gene,function(x){
    funcEnrich.Fisher(input_list=x,bg_list=bg_list,gs2gene=gs2gene,use_gs=use_gs,
                      min_gs_size=min_gs_size,max_gs_size=max_gs_size,
                      Pv_adj='none',Pv_thre=Pv_thre)
  })
  names(f_res) <- names(target_gene)
  return(f_res)
}
##

######################################
filter_target_top <- function(net1,source_col=1,target_col=2,value_col=3,top_num=3){
  all_s <- net1[,source_col]
  new_net <- c()
  for(i in unique(all_s)){
    w1 <- which(all_s==i)
    net2 <- net1[w1,]
    if(length(w1)>1){
      net2 <- net2[order(as.numeric(net2[,value_col]),decreasing = TRUE),]
      if(length(w1)>top_num) net2 <- net2[1:top_num,]
    }
    new_net <- rbind(new_net,net2)
  }
  return(new_net)
}
get_gs2mat <- function(use_gs,gene,min_size=30,max_size=300){
  w1 <- unlist(lapply(use_gs,length))
  use_gs <- use_gs[which(w1>=min_size & w1<=max_size)]
  use_gs1 <- lapply(use_gs,function(x)intersect(x,gene))
  w1 <- unlist(lapply(use_gs1,length))
  use_gs1 <- use_gs1[which(w1>0)]
  all_g <- unique(unlist(use_gs1))
  mat1 <- matrix(0,nrow=length(all_g),ncol=length(use_gs1))
  rownames(mat1) <- all_g; colnames(mat1) <- names(use_gs1)
  for(i in names(use_gs1)){
    mat1[use_gs1[[i]],i] <- 1
  }
  mat1
}
## this is not function!!!
draw.geneSet2Driver2geneSet <- function(r0,r1,main=main){
  cc <- get_transparent(get.class.color(i),0.05)
  # input r0 r1
  n1 <- r0$list_net; n2 <- r1$list_net; s2 <- unique(n2[,2]);
  n1 <- n1[which(n1[,2] %in% s2),]; s1 <- unique(n1[,1]); s3 <- unique(n2[,1]);
  s1 <- rev(s1); s2 <- rev(s2); s3 <- rev(s3);
  l1 <- unlist(lapply(r0$gs_list[s1],length)); l2 <- unlist(lapply(r0$gene_list[s2],length)); l3 <- unlist(lapply(r1$gs_list[s3],length))
  ls1 <- 1+cumsum(l1+1); ls2 <- 1+cumsum(l2+1); ls3 <- 1+cumsum(l3+1);
  plot(1,col='white',xlim=c(-2,6),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',bty='n',main=main)
  y1 <- seq(0,1,length.out=1+sum(l1)+length(l1)); y11 <- c(y1[1],y1[ls1]); dy1 <- y1[2]-y1[1]; y111<-dy1/2+y11[1:(length(y11)-1)]/2+y11[2:length(y11)]/2;
  y2 <- seq(0,1,length.out=1+sum(l2)+length(l2)); y22 <- c(y2[1],y2[ls2]); dy2 <- y2[2]-y2[1]; y222<-dy2/2+y22[1:(length(y22)-1)]/2+y22[2:length(y22)]/2;
  y3 <- seq(0,1,length.out=1+sum(l3)+length(l3)); y33 <- c(y3[1],y3[ls3]); dy3 <- y3[2]-y3[1]; y333<-dy3/2+y33[1:(length(y33)-1)]/2+y33[2:length(y33)]/2;
  names(y111) <-s1;names(y222) <-s2;names(y333) <-s3;
  cc2 <- sample(colorRampPalette(brewer.pal(9,'Paired'))(length(s2))); cc22 <- get_transparent(cc2,0.2)
  cc2 <- get_transparent(cc2,0.8); names(cc2) <- s2;
  #
  x1 <- c(-2,1.1,2.25); x2 <- c(0.75,1.4,6);
  segments(x0=x2[1],x1=x1[2],y0=y111[n1[,1]],y1=y222[n1[,2]],col=cc2[n1[,2]],lwd=as.numeric(n1[,3]))
  arrows(x0=x2[2],x1=x1[3],y1=y333[n2[,1]],y0=y222[n2[,2]],col=cc2[n2[,2]],lwd=as.numeric(n2[,3])-z_thre1+0.25,length=0.05)
  rect(xleft = x1[1],xright = x2[1],ybottom = y11[1:(length(y11)-1)]+dy1,ytop=y11[2:length(y11)],border = 'light grey',lwd=0.2,col=cc)
  rect(xleft = x1[2],xright = x2[2],ybottom = y22[1:(length(y22)-1)]+dy2,ytop=y22[2:length(y22)],border = 'light grey',lwd=0.2,col=cc22)
  rect(xleft = x1[3],xright = x2[3],ybottom = y33[1:(length(y33)-1)]+dy3,ytop=y33[2:length(y33)],border = 'light grey',lwd=0.2,col=cc)
  for(j in 1:length(s1)){
    text(x2[1],seq(from=y11[1:(length(y11)-1)][j]+1.5*dy1,to=-0.5*dy1+y11[2:length(y11)][j],
                   length.out=length(r0$gs_list[[s1[j]]])),r0$gs_list[[s1[j]]],cex=0.5,pos=2,xpd=TRUE)
  }
  for(j in 1:length(s3)){
    ggs1 <- r1$gs_list[[s3[j]]];
    rgg  <- gs_Z_rank[ggs1,grep(i,colnames(gs_Z_rank)),drop=F]; if(length(ggs1)==1) rgg <- sign(sum(rgg)) else rgg <- sign(rowSums(rgg))
    cc1 <- rep('grey',length.out=length(rgg)); cc1[which(rgg==1)] <- 'red'; cc1[which(rgg== -1)] <- 'blue';
    text(x1[3],seq(from=y33[1:(length(y33)-1)][j]+1.5*dy3,to=-0.5*dy3+y33[2:length(y33)][j],
                   length.out=length(r1$gs_list[[s3[j]]])),r1$gs_list[[s3[j]]],cex=0.5,pos=4,xpd=TRUE,col=cc1)
  }
  for(j in 1:length(s2)){
    gg <- r0$gene_list[[s2[j]]]; gg <- sort(gg,decreasing = TRUE);
    text(x1[2]/2+x2[2]/2,seq(from=y22[1:(length(y22)-1)][j]+1.5*dy2,to=-0.5*dy2+y22[2:length(y22)][j],
                             length.out=length(gg)),gg,
         cex=0.4,adj=0.5,xpd=TRUE)
    #if(length(gg)%%2==0){gg1 <- gg[1:(length(gg)/2)];gg2 <- gg[(length(gg)/2+1):length(gg)];}
    #if(length(gg)%%2==1){gg1 <- gg[1:((1+length(gg))/2)];gg2 <- gg[(((1+length(gg))/2)+1):length(gg)];}
    #text(x1[2]/2+x2[2]/2,seq(from=y22[1:(length(y22)-1)][j]+1.5*dy2,to=-0.5*dy2+y22[2:length(y22)][j],
    #                         length.out=length(gg1)),gg1,
    #     cex=0.4,adj=0,xpd=TRUE)
    #text(x1[2]/2+x2[2]/2,seq(from=y22[1:(length(y22)-1)][j]+1.5*dy2,to=-0.5*dy2+y22[2:length(y22)][j],
    #                         length.out=length(gg2)),gg2,
    #     cex=0.4,adj=1,xpd=TRUE)
  }
  text(x2[1],y1[1],'Enriched Gene Sets for Drivers',pos=2)
  text(x1[2]/2+x2[2]/2,y2[1],'Drivers')
  text(x1[3],y3[1],'Enriched Gene Sets for Drivers\' Targets',pos=4)
}


#### calculate activity
# mode: out, downstream targets
# mode: all, all related as targets
cal.Activity.mod <- function(target_list=NULL, igraph_obj = NULL, net_dat=NULL,cal_mat=NULL, es.method = 'weightedmean',
                             std=TRUE,memory_constrain=FALSE,mode='out') {
  #
  all_input_para <- c('cal_mat','es.method','std','memory_constrain')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('memory_constrain',c(TRUE,FALSE),envir=environment()),
                 check_option('std',c(TRUE,FALSE),envir=environment()),
                 check_option('es.method',c('mean','weightedmean','maxmean','absmean'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(target_list)==TRUE & is.null(igraph_obj)==TRUE & mode=='out'){
    message('Either target_list or igraph_obj is required, please check and re-try!');return(FALSE);
  }
  if(is.null(target_list)==FALSE & memory_constrain==TRUE){
    ac.mat <- cal.Activity.old(target_list=target_list, cal_mat=cal_mat, es.method = es.method,std=std)
    return(ac.mat)
  }
  if(is.null(target_list)==TRUE & memory_constrain==TRUE){
    message('Only accepts target_list input when memory_constrain=TRUE, please check and re-try!');return(FALSE);
  }
  if(nrow(cal_mat)==0){
    message('No genes in the cal_mat, please check and re-try!');return(FALSE);
  }
  if(std==TRUE) cal_mat <- apply(cal_mat, 2, do.std)
  if(mode=='all'){
    if(is.null(net_dat)==FALSE){
      igraph_obj   <- graph_from_data_frame(net_dat[,c('source','target')],directed=FALSE) ## add edge weight ???
      igraph_obj   <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
      igraph_obj   <- set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
    }else{
      message('need to input net_dat for mode set to all!')
    }
  }
  if(is.null(igraph_obj)==FALSE){
    gr <- igraph_obj
    mat1 <- get_igraph2matrix(gr,es.method=es.method)
    mat2 <- get_igraph2matrix(gr,es.method='mean')
    all_source <- get_gr2driver(gr,mode=mode)
    print(str(all_source))
  }else{
    if(is.null(target_list)==FALSE){
      mat1 <- get_target_list2matrix(target_list,es.method=es.method)
      mat2 <- get_target_list2matrix(target_list,es.method='mean')
      all_source <- names(target_list)
    }
  }
  ##
  mat1_source <- mat1[all_source,]
  w1 <- base::intersect(rownames(cal_mat),colnames(mat1_source))
  print(str(w1))
  if(base::length(w1)==0){
    message('No intersected genes cound for the cal_mat and target in the network, please check and re-try!');
    return(FALSE)
  }
  use_mat1_source <- mat1_source[,w1] ## network info
  use_mat2_source <- mat2[all_source,w1] ## network binary info

  ## weighted mean + mean
  if(es.method %in% c('weightedmean','mean')){
    use_cal_mat <- cal_mat[w1,] ## expression info
    out_mat <- use_mat1_source %*% use_cal_mat
    out_mat <- out_mat/Matrix::rowSums(use_mat2_source) ## get mean
  }
  #return(as.matrix(out_mat))
  ## absmean
  if(es.method == 'absmean'){
    use_cal_mat <- cal_mat[w1,] ## expression info
    out_mat <- use_mat1_source %*% abs(use_cal_mat)
    out_mat <- out_mat/Matrix::rowSums(use_mat2_source) ## get mean
  }
  ## maxmean
  if(es.method == 'maxmean'){
    use_cal_mat <- cal_mat[w1,] ## expression info
    use_cal_mat_pos <- use_cal_mat;use_cal_mat_pos[which(use_cal_mat_pos<0)] <- 0;
    use_cal_mat_neg <- use_cal_mat;use_cal_mat_neg[which(use_cal_mat_neg>0)] <- 0;
    out_mat_pos  <- use_mat1_source %*% use_cal_mat_pos
    out_mat_pos  <- out_mat_pos/Matrix::rowSums(use_mat2_source) ## get mean
    out_mat_neg  <- use_mat1_source %*% use_cal_mat_neg
    out_mat_neg  <- out_mat_neg/Matrix::rowSums(use_mat2_source) ## get mean
    out_mat_sign <- sign(abs(out_mat_pos)-abs(out_mat_neg))
    out_mat_sign_pos <- out_mat_sign; out_mat_sign_pos[out_mat_sign_pos!=1] <-0;
    out_mat_sign_neg <- out_mat_sign; out_mat_sign_neg[out_mat_sign_neg!= -1] <-0;
    out_mat <- out_mat_pos*out_mat_sign_pos-out_mat_neg*out_mat_sign_neg
  }
  ## median, min, max , not supported
  # output
  ac.mat <- as.matrix(out_mat)
  print(str(ac.mat))
  w1 <- which(is.na(ac.mat[,1])==FALSE)
  print(str(w1))
  if(base::length(w1)==0){
    message('Fail in calculating activity, please check the ID type in cal_mat and target_list and try again !')
  }
  ac.mat <- ac.mat[w1,]
  return(ac.mat)
}
##
get_net2target_list.all <- function(net_dat){
  all_input_para <- c('net_dat')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  all_source <- base::unique(c(net_dat$source,net_dat$target))
  all_target <- lapply(all_source, function(x) {
    print(x)
    t1 <- net_dat[which(net_dat$source == x), base::intersect(c('target', 'MI', 'spearman'),colnames(net_dat))]
    t2 <- net_dat[which(net_dat$target == x), base::intersect(c('source', 'MI', 'spearman'),colnames(net_dat))]
    colnames(t2) <- colnames(t1)
    tmp1 <- rbind(t1,t2)
    tmp2 <- lapply(unique(tmp1[,1]),function(x1){
      x2 <- tmp1[which(tmp1[,1]==x1),]
      x2[which.max(x2[,2]),]
    })
    tmp2 <- do.call(rbind,tmp2)
    rownames(tmp2) <- tmp2[,1]
    return(tmp2)
  })
  names(all_target) <- all_source
  return(all_target)
}
###### modify function enrichment results
## prepare gs_name2info
gs_name2info <- lapply(names(all_gs2gene),function(x){
  print(x)
  w1 <- which(all_gs2gene_info$Category == x | all_gs2gene_info$`Sub-Category`==x)
  if(length(w1)==1) return(cbind(all_gs2gene_info[w1,],names(all_gs2gene[[x]]),stringsAsFactors=F))
})
gs_name2info <- do.call(rbind,gs_name2info)
colnames(gs_name2info)[7] <- 'ID'
##
modify_funcEnrichRes <- function(res,gs_name2info){
  res1 <- base::merge(x=res,y=gs_name2info,by.x='#Name',by.y='ID')
  return(res1)
}












