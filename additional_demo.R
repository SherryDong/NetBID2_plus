## This script aims to do following analysis ##
############# preparations ####################

library(NetBID2)
source('/Volumes/project_space/Network_JY/yu3grp/NetBID2/NetBID2-R/R//pipeline_functions.R') ## need to source the inner functions
source('../additional_functions/additional_functions.R')
#source('/Volumes/project_space/Network_JY/yu3grp/NetBID2//NetBID2-R_internal/additional_functions.R')


# set the directory to get the pre-saved RData
analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")

############## I: driver Vs. gene set Vs. boxplot (draw.gsBygroup(),draw.topdriver2gs(),draw.driverTargetFuncEnrich())
# RData for ms table
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
analysis.par$out.dir.PLOT <- '../additional_functions/PLOT' ## directory to save the plot figures
##
mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
                  Group3=c('IMPG2','GABRA5','EGFL11','NRL','MAB21L2','NPR3','MYC'),
                  Group4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
mark_col <- get.class.color(names(mark_gene))
##
db.preload(use_level='gene',use_spe='human',update=FALSE)
gs.preload()
exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames must be the originalID

## get expression matrix for the transfered gene name, if original is gene-based expression matrix, just use the exp_mat
exp_mat_gene <- exp_mat
## calculate activity for all genesets
use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
ac_gs <- cal.Activity.GS(use_gs2gene = c(use_gs2gene,mark_gene),cal_mat = exp_mat_gene)
## get DA
phe_info <- pData(analysis.par$cal.eset)
DA_gs<-list()
for(each in c('WNT','SHH','G4')){
  G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each)] # get sample list for G0
  G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each)] # get sample list for G1
  DA_gs[[each]] <- getDE.limma.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,G1_name=each,G0_name='others')
}
use_z_mat <- lapply(DA_gs,function(x)x[DA_gs$G4$ID,'Z-statistics'])
use_z_mat <- do.call(cbind,use_z_mat)
rownames(use_z_mat) <- DA_gs$G4$ID

## function: draw.gsBygroup()
use_p=get_obs_label(phe_info,'subgroup')
use_p <- use_p[order(factor(use_p,levels=c('WNT','SHH','G4')))]
draw.gsBygroup(mark_gene,exp_mat_gene[,names(use_p)],
               use_p=use_p,
               use_z_mat = use_z_mat,scale_row=TRUE)

## function: draw.topdriver2gs()
all_top <- lapply(analysis.par$DA,function(x)x$ID[1:200])
draw.topdriver2gs(all_top,mark_gene,low_Sum = 1)
all_top <- lapply(analysis.par$DA,function(x)x$ID[1:100])
draw.topdriver2gs(all_top,all_gs2gene[['H']],low_Sum = 5)
##
use_z_mat <- lapply(analysis.par$DA,function(x)x[analysis.par$DA$`G4.Vs.WNT`$ID,'Z-statistics'])
use_z_mat <- do.call(cbind,use_z_mat)
rownames(use_z_mat) <- analysis.par$DA$`G4.Vs.WNT`$ID

## function: draw.driverTargetFuncEnrich()
draw.driverTargetFuncEnrich(use_driver = all_top$G4.Vs.others[1:5],
                            target_list = analysis.par$merge.network$target_list,
                            ms_tab=analysis.par$final_ms_tab,
                            use_gs2gene = use_gs2gene,
                            exp_mat=exp_mat_gene,
                            use_ac_mat=exprs(analysis.par$merge.ac.eset),
                            use_p=use_p,use_z_mat=use_z_mat,top_gs_num=8)


############## II: SINBA plot for synergistic effect (draw.GSEA.NetBID.SINBA())
ms_tab <- analysis.par$final_ms_tab ## get the master table data frame
ms_tab <- ms_tab[which(ms_tab$Size>=30 & ms_tab$Size <=1000),]
comp_name <- 'G4.Vs.others' ## get the comparison name
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                               Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.4,Pv_thre=1e-8,
                               main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,label_cex = 1)

seed_driver <- sig_driver[which(sig_driver$logFC.G4.Vs.others_DA>0)[1],1]
part_driver <- ms_tab$originalID_label

######
comp_name <- 'G4.Vs.others' ## each comparison must give a name !!!
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # get sample list for G0
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # get sample list for G1
##
merge_target <- lapply(part_driver,function(x){
  m1 <- merge_target_list(driver1=seed_driver,driver2=x,target_list=analysis.par$merge.network$target_list)
})
names(merge_target) <- part_driver
ac_combine_mat <- cal.Activity(target_list=merge_target,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
DA_driver_combine <- getDE.BID.2G(eset=generate.eset(ac_combine_mat),G1=G1,G0=G0,G1_name='G4',G0_name='others')
ori_part_Z <- analysis.par$DA[[comp_name]][part_driver,'Z-statistics']
ori_seed_Z <- analysis.par$DA[[comp_name]][seed_driver,'Z-statistics']
#diff_Z <- 2*DA_driver_combine[part_driver,'Z-statistics']-(ori_part_Z+ori_seed_Z)
#names(diff_Z) <- part_driver
comp <- 'G4.Vs.others'
DE <- analysis.par$DE[[comp]]
driver_DA_Z <- analysis.par$DA[[comp_name]][,'Z-statistics']
names(driver_DA_Z) <- rownames(analysis.par$DA[[comp_name]])
driver_DE_Z <- analysis.par$DE[[comp_name]][,'Z-statistics']
names(driver_DE_Z) <- rownames(analysis.par$DE[[comp_name]])
DA_Z_merge <- DA_driver_combine[,'Z-statistics']
names(DA_Z_merge) <- rownames(DA_driver_combine)
target_list_merge <- merge_target
seed_driver_label <- ms_tab[seed_driver,'gene_label']
partner_driver_list <- part_driver
profile_col <- 't'
partner_driver_label <- ms_tab[partner_driver_list,'gene_label']
target_list <- analysis.par$merge.network$target_list

## function: draw.GSEA.NetBID.SINBA()
draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=20,profile_trend='pos2neg',top_order='merge',Z_sig_thre = 1.64,
                       target_nrow=1,target_col='RdBu',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo1.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=15,profile_trend='pos2neg',top_order='merge',Z_sig_thre = 1.64,
                       target_nrow=2,target_col='RdBu',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo2.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=15,profile_trend='pos2neg',top_order='diff',Z_sig_thre = 1.64,
                       target_nrow=2,target_col='RdBu',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo3.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=5,profile_trend='pos2neg',top_order='merge',Z_sig_thre = 1.64,
                       target_nrow=1,target_col='black',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo4.pdf',analysis.par$out.dir.PLOT))


## Part III. cal.Activity.mod(), use mode="all" will use all edges, and calculate for all genes
# the network file must not be marked by _TF/_SIG, otherwise only pre-defined drivers will be calculated
if(exists('analysis.par')==TRUE) rm(analysis.par)
network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
network.project.name <- 'project_2019-02-14' # demo project name
project_main_dir <- 'test/'
project_name <- 'test_driver'
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
                                            project_name=project_name,
                                            network_dir=network.dir,
                                            network_project_name=network.project.name)
tf_net <- get.SJAracne.network(analysis.par$tf.network.file)
sig_net <- get.SJAracne.network(analysis.par$sig.network.file)
net_dat <- rbind(tf_net$network_dat,sig_net$network_dat)
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
ac_mat <- cal.Activity.mod(net_dat = net_dat,
                       cal_mat=Biobase::exprs(analysis.par$cal.eset),
                       es.method='maxmean',mode='all')


## Part IV. target network with defined direction
# draw.targetNet.forDriver <- function(use_driver,use_driver_label=use_driver,network_dat=NULL,transfer_tab=NULL,use_label=TRUE){
network_dat <- analysis.par$merge.network$network_dat
use_driver <- 'TP53'
transfer_tab <- data.frame('external_gene_name'=rownames(exprs(analysis.par$cal.eset)),'external_gene_name.1'=rownames(exprs(analysis.par$cal.eset)),stringsAsFactors = F)
draw.targetNet.forDriver(use_driver,network_dat=network_dat,transfer_tab=transfer_tab,n_layer=2)

## Part V. well, complicated function, draw.geneSet2Driver2geneSet()
db.preload(use_level = 'gene')
gs.preload()
#
use_ms_tab <- analysis.par$final_ms_tab
g1 <- use_ms_tab$originalID_label[1:30] ## draw for top 30
g0 <- use_ms_tab$geneSymbol[1:30] ## gene ysmbol
transfer_tab <- get_IDtransfer2symbol2type(rownames(exprs(analysis.par$cal.eset)),from_type = 'external_gene_name')
z_thre1 <- 2.36
# sig_gs
use_gs <- merge_gs(all_gs2gene,use_gs=c('BP','CP:BIOCARTA','H','CP:REACTOME','CP:KEGG'))
use_gs_Z <- cal.Activity.GS(use_gs,cal_mat=exprs(analysis.par$cal.eset))
phe1 <- pData(analysis.par$cal.eset)
gs_DA <- getDE.limma.2G(eset=generate.eset(use_gs_Z,phe1),
                        G1=rownames(phe1)[which(phe1$subgroup=='G4')],G0=rownames(phe1)[which(phe1$subgroup!='G4')])
gs_Z_rank <- as.matrix(rank(gs_DA$`Z-statistics`)-nrow(gs_DA)/2);rownames(gs_Z_rank) <- gs_DA$ID;colnames(gs_Z_rank) <- i
#
i <- 'G4.Vs.others'
driver2geneset <- draw.bubblePlot(driver_list = g1,target_list=analysis.par$merge.network$target_list,
                                  Z_val=use_ms_tab[g1,sprintf('Z.%s_DA',i)],use_gs=c('BP','CP:BIOCARTA','H','CP:REACTOME','CP:KEGG'),
                                  min_gs_size=30,max_gs_size=300,transfer2symbol2type=transfer_tab,top_driver_number = length(g1),
                                  top_geneset_number = 100,
                                  only_return_mat=TRUE) ## time consuming
mat1 <- get_gs2mat(use_gs,g0,30,300); mat1 <- mat1[,which(colSums(mat1)>=round(length(g0)/100))] # get gene set annotation for top drivers
d2s <- g0;names(d2s) <- g1; # get driver to gene set transfer vector
p1 <- driver2geneset; # row driver, column gene sets
rownames(p1) <- as.character(d2s[rownames(p1)]) # transfer driver name to gene symbol
p1 <- p1[,which(apply(p1,2,max)>z_thre1)];  w1 <- intersect(rownames(mat1),rownames(p1))
p1 <- p1[w1,]; mat1 <- mat1[w1,]
## transfer matrix to cluster
r1 <- mat2cluster(mat1=t(p1),z_thre = z_thre1,gene_list=NULL,h_gene=0.5,h_gs=0.5,strategy='max')
r0 <- mat2cluster(mat1=t(mat1),z_thre=0.5,gene_list=r1$gene_list[unique(r1$list_net[,'gene_list'])],
                  strategy='max',h_gene=0.5,h_gs=0.5)
# get top
r0$list_net <- filter_target_top(r0$list_net,source_col=2,target_col=1,value_col=3,top_num=1)
r0$list_net <- filter_target_top(r0$list_net,source_col=1,target_col=2,value_col=3,top_num=2)
r1$list_net <- filter_target_top(r1$list_net,source_col=2,target_col=1,value_col=3,top_num=3) ## driver's top
r1$list_net <- filter_target_top(r1$list_net,source_col=1,target_col=2,value_col=3,top_num=2) ## target gene set's top
r0$list_net <- r0$list_net[order(factor(r0$list_net[,2], levels= unique(r1$list_net[,2]))),];
draw.geneSet2Driver2geneSet(r0,r1,main='test')
#


##### modify function enrichment results
analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
ms_tab <- analysis.par$final_ms_tab
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
                               logFC_col='logFC.G4.Vs.others_DA',
                               Pv_col='P.Value.G4.Vs.others_DA',
                               logFC_thre=0.4,
                               Pv_thre=1e-7,
                               main='Volcano Plot for G4.Vs.others_DA',
                               show_label=FALSE,
                               label_type = 'origin',
                               label_cex = 0.5)
gs.preload(use_spe='Homo sapiens',update=FALSE)
res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
                          bg_list=ms_tab[,'geneSymbol'],
                          use_gs=c('H','C5'),
                          Pv_thre=0.1,Pv_adj = 'none')
new_res1 <- modify_funcEnrichRes(res=res1,gs_name2info=gs_name2info)




