library(uwot)
library(dplyr)
library(Matrix)
library(cowplot)
library(data.table)
library(magrittr)
library(tidyverse)
library(RCurl)
library(purrr)
library(pheatmap)
library(MAST)
library(RColorBrewer)
library(grid)
library(gtable)
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(harmony)
library(clustree)
library(cowplot)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(ggrepel)
##################################
levels(MergeF$celltype)
APC <- subset(MergeF, idents= c("moDC", "Mono", "migDC", "LC", "Mac_1", "Mac_2"))
Idents(APC) <- APC$celltype
APC$celltype <- droplevels(APC$celltype)

DimPlot(APC, label = T, label.size = 5, group.by = "celltype")+NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(6))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=15),
        axis.title.y=element_text(colour='black', size=15),
        axis.text=element_text(colour='black',size=15),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=15),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


saveRDS(APC, "H:/BP data analysis/Integrated with SI data/APC/APC.rds")
###################################Percentage change
#分组箱线图1
my_comparisons = list(c("HC", "AD"), c("HC", "Pso"), c("HC", "BP"), c("AD", "BP"), c("AD", "Pso"), c("Pso", "BP")) 

Singlecellratio_plotstat(APC, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'group',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =6)

########################DEGs
Ident <- c("moDC", "Mono", "migDC", "LC", "Mac_1", "Mac_2")
set1 = list()
for (i in Ident){
  tryCatch({
    print(i)
    gc() 
    set1[[paste0(i)]] <- FindMarkers(APC, ident.1 = "BP", ident.2 ="NML", group.by = "dis", verbose =TRUE,
                                         assay = "RNA", slot = "data", subset.ident = i, test.use = "MAST")
  }, error=function(e){cat("ERROR :",GroupitionMessage(e), "\n")})
}
library(openxlsx)
write.xlsx(set1, file = "H:/BP data analysis/Integrated with SI data/APC/BPvsNML_DEG.xlsx", rowNames = TRUE, colNames=TRUE)


Ident <- c("moDC", "Mono", "migDC", "LC", "Mac_1", "Mac_2")
set1 = list()
for (i in Ident){
  tryCatch({
    print(i)
    gc() 
    set1[[paste0(i)]] <- FindMarkers(APC, ident.1 = "AD", ident.2 ="NML", group.by = "group", verbose =TRUE,
                                     assay = "RNA", slot = "data", subset.ident = i, test.use = "MAST")
  }, error=function(e){cat("ERROR :",GroupitionMessage(e), "\n")})
}
library(openxlsx)
write.xlsx(set1, file = "H:/BP data analysis/Integrated with SI data/APC/ADvsHC_DEG.xlsx", rowNames = TRUE, colNames=TRUE)


Ident <- c("moDC", "Mono", "migDC", "LC", "Mac_1", "Mac_2")
set1 = list()
for (i in Ident){
  tryCatch({
    print(i)
    gc() 
    set1[[paste0(i)]] <- FindMarkers(APC, ident.1 = "Pso", ident.2 ="NML", group.by = "group", verbose =TRUE,
                                     assay = "RNA", slot = "data", subset.ident = i, test.use = "MAST")
  }, error=function(e){cat("ERROR :",GroupitionMessage(e), "\n")})
}
library(openxlsx)
write.xlsx(set1, file = "H:/BP data analysis/Integrated with SI data/APC/PsovsHC_DEG.xlsx", rowNames = TRUE, colNames=TRUE)

#######################################Volcanoplot#########################
library(openxlsx)
library(ggVolcano)
nrDEG <- read.xlsx('H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx', colNames=T, rowNames=T, sheet = 1)
colnames(nrDEG)
add_regulate <- function (data, log2FC_name = "avg_log2FC", fdr_name = "p_val_adj", 
                          log2FC = 0.1, fdr = 0.05) 
{
  colnames(data)[colnames(data) == log2FC_name] <- "avg_log2FC"
  colnames(data)[colnames(data) == fdr_name] <- "p_val_adj"
  data$regulate <- "NS"
  loc_up <- intersect(which(data$avg_log2FC > log2FC), 
                      which(data$p_val_adj < fdr))
  loc_down <- intersect(which(data$avg_log2FC < (-log2FC)), 
                        which(data$p_val_adj < fdr))
  data$regulate[loc_up] <- "Up"
  data$regulate[loc_down] <- "Down"
  return(data)
}
data <- add_regulate(nrDEG, log2FC_name="avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.3, fdr = 0.05)
data$row <- rownames(data)

for_label <- data %>%   filter(abs(avg_log2FC) > 4 & (-log10(p_val_adj))> -log10(0.05))

ggplot(data = data,aes(x = avg_log2FC,y = (-log10(p_val_adj)))) +
  geom_point(alpha=0.4, size=3.5,aes(color=regulate)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw()+  
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(aes(label = row),data = for_label,color="black", max.overlaps = 20)+
  ggtitle("DEGs in moDC")

###################
setwd("H:/BP data analysis/Integrated with SI data/APC/")
#差异基因====================================================================
library(Seurat)
library(tibble)
library(ggplot2)
library(ggrepel)


KS_scRNA_multiVlnvo_plot <- function(Seurat_object,#单细胞seurat对象
                                     DEGs_outdir,#差异基因结果保存文件路径
                                     DEGs_list=F,#是否已经已经进行了DEGs分析
                                     dges,#通过KS_scRNA_multiVlnvo_plot函数分析得到得cluster的差异基因list
                                     min.pct,#FindMarkers差异基因分析参数
                                     logfc.threshold,#FindMarkers差异基因分析参数
                                     test.use,#FindMarkers差异基因分析参数
                                     group,#FindMarkers差异基因分析参数
                                     ident.1,#FindMarkers差异基因分析参数
                                     ident.2,#FindMarkers差异基因分析参数
                                     logFC_cut,#需要展示的基因logFC阈值
                                     top_gene=T,#T，表示展示差异基因top10
                                     label_gene,#top_gene=F,这里展示需要展示的基因
                                     xlabel=NULL, #x轴cluster标签，F为默认。可自行设置，顺序与idents_level一致
                                     idents_level=NULL,#x轴cluster顺序，F为默认。可自行设置，顺序与xlabel一致
                                     point_color=NULL,#散点颜色（区分显著与不显著）
                                     box_color=NULL,#x轴cluster的底色
                                     cellText_color=NULL,#x轴cluster标签的颜色
                                     text_size,#标记基因的大小
                                     height#x轴cluster框的高度
){
  require(Seurat)
  require(tibble)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  require(dittoSeq)
  
  cluster <- as.character(unique(Idents(Seurat_object)))  
  
  if(DEGs_list==F){
    
    DEGs <- list()
    
    for (i in 1:length(cluster)) {
      
      sce = Seurat_object[, Idents(Seurat_object) %in% cluster[i]]
      
      dges <- FindMarkers(sce, min.pct = min.pct, 
                          logfc.threshold = logfc.threshold,
                          group.by = group,
                          ident.1 = ident.1,
                          ident.2= ident.2,
                          test.use= test.use)
      dges = as_tibble(cbind(gene = rownames(dges), dges))
      dges$idents <- cluster[i]
      
      DEGs[[i]] <- dges
      names(DEGs)[i] <- paste0(cluster[i],"_DEGs")
      
    }
    
    
    saveRDS(DEGs, file = "sce_DEGs.rds")
    
    
    diff <- do.call(rbind, DEGs)
    
  }else{
    
    diff <- do.call(rbind, dges)
    
  }
  
  
  
  
  
  #显著性
  diff <- diff[which(abs(diff$avg_log2FC)>=logFC_cut),]
  diff$label <- ifelse(diff$p_val_adj<=0.05,"sig","unsig")
  
  
  if(top_gene==T){
    
    ##label gene
    top_gene = diff %>% group_by(idents) %>% top_n(n = 10, wt = abs(avg_log2FC))
    
  }else{
    
    top_gene = diff[diff$gene %in% label_gene,]
    
  }
  
  
  ###bakground
  max_fc <- c()
  min_fc <- c()
  for (i in 1:length(cluster)) {
    
    df <- subset(diff, idents == cluster[i])
    
    df_max <- max(df$avg_log2FC)
    df_min <- min(df$avg_log2FC)
    max_fc <- append(max_fc, df_max)
    min_fc <- append(min_fc, df_min)
  }
  
  bakground_FC <- data.frame(x=cluster,
                             max_fc = max_fc+0.2,
                             min_fc = min_fc-0.2)
  
  
  
  ###cell label
  if(is.null(xlabel)){
    
    cluster_label <-data.frame(x=cluster,
                               y=0,
                               label=cluster)
  }else{
    
    cluster_label <-data.frame(x=cluster,
                               y=0,
                               label=xlabel)
    
  }
  
  
  if(is.null(idents_level)){
    
    diff = diff
    
  }else{
    
    diff$idents <- factor(diff$idents, levels = idents_level)
    bakground_FC$x <- factor(bakground_FC$x, levels = idents_level)
    
  }
  
  #作图
  if(is.null(point_color)){
    
    point_color = c("#E64B357F","#00A0877F")
  }
  
  if(is.null(box_color)){
    
    box_color = dittoColors()[1:length(cluster)]
  }
  
  
  if(is.null(cellText_color)){
    
    cellText_color="black"
  }
  
  
  p <- ggplot()+
    geom_col(data = bakground_FC,
             mapping = aes(x = x,y = max_fc),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = bakground_FC,
             mapping = aes(x = x,y = min_fc),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = diff,
                aes(x = idents, y = avg_log2FC, color = label),
                size = 1,
                width =0.4)+
    scale_color_manual(values = point_color)+
    geom_tile(data = cluster_label,
              aes(x=x,y=y),
              height=height,
              color = "black",
              show.legend = F,
              fill=box_color,
              linewidth=0.5)+
    geom_text(data=cluster_label,
              aes(x=x,y=y,label=label),
              size =4,
              fontface="bold",
              color = cellText_color)+
    geom_text_repel(data=top_gene,
                    aes(x=idents,y=avg_log2FC,label=gene),
                    color="black", size=text_size, fontface="bold.italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")),
                    box.padding = 0.2,
                    point.padding = 0.3, 
                    segment.color = 'black', 
                    segment.size = 0.3, force = 1, 
                    max.iter = 3e3)+
    theme_minimal()+
    labs(y="average logFC")+
    theme(axis.title.y = element_text(size = 12,
                                      color = "black",
                                      face = "bold"),
          axis.title.x = element_blank(),
          axis.line.y = element_line(color = "black",
                                     size = 0.8),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position.inside =  c(10, 9),
          legend.direction = "horizontal",
          legend.justification = c(1,0),
          legend.text = element_text(size = 12))
  return(p)
  
}


###########################laod DEGs with multiple sheets########################
multiplesheets <- function(fname) { 
  
  # getting info about all excel sheets 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  
  # assigning names to data frames 
  names(data_frame) <- sheets 
  
  # print data frame 
  print(data_frame) 
} 

# specifying the path name 
path <- "H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx"
sce_DEGs <- multiplesheets(path)

#############################multiple-Vlnplot
KS_scRNA_multiVlnvo_plot(Seurat_object = APC,
                         DEGs_list=T,
                         dges = sce_DEGs,
                         logFC_cut = 1,
                         text_size = 3,
                         xlabel = NULL,
                         idents_level = NULL,
                         height = 1.5,
                         top_gene = F,
                         label_gene = NULL #c("CCL26", "CCL14", "CCL24", "CCL27", "CCL23", "IL13", "IL18", "CXCL5","CXCL8","CXCL12")
)+ggtitle("DEGs of BP vs HC")


#########################AD vs HC
path <- "H:/BP data analysis/Integrated with SI data/APC/ADvsHC_DEG.xlsx"
sce_DEGs <- multiplesheets(path)

#############################multiple-Vlnplot
KS_scRNA_multiVlnvo_plot(Seurat_object = APC,
                         DEGs_list=T,
                         dges = sce_DEGs,
                         logFC_cut = 1,
                         text_size = 3,
                         xlabel = NULL,
                         idents_level = NULL,
                         height = 1.5,
                         top_gene = F,
                         label_gene = NULL #c("CCL26", "CCL14", "CCL24", "CCL27", "IL13", "CXCL5","CXCL8","CXCL12")
)+ggtitle("DEGs of AD vs HC")

#########################Pso vs HC
path <- "H:/BP data analysis/Integrated with SI data/APC/PsovsHC_DEG.xlsx"
sce_DEGs <- multiplesheets(path)

#############################multiple-Vlnplot
KS_scRNA_multiVlnvo_plot(Seurat_object = APC,
                         DEGs_list=T,
                         dges = sce_DEGs,
                         logFC_cut = 1,
                         text_size = 3,
                         xlabel = NULL,
                         idents_level = NULL,
                         height = 1.5,
                         top_gene = F,
                         label_gene = NULL
)+ggtitle("DEGs of Pso vs HC")

#################################geneX expression across the cell types########################
#########Mac1#############
#######load Ad and trim dataframe
AD_Mac1 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/ADvsHC_DEG.xlsx", sheet = 5, colNames = T, rowNames = T)
AD_Mac1 <- AD_Mac1[,c(2,5)]
colnames(AD_Mac1) <- c("AD", "p_val_adj_1")

#######load BP and trim dataframe
BP_Mac1 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 5, colNames = T, rowNames = T)
BP_Mac1 <- BP_Mac1[,c(2,5)]
colnames(BP_Mac1) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(AD_Mac1, BP_Mac1, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE AD",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$AD) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE AD" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(AD, BP, color=threshold))+
  geom_text_repel(aes(AD, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("Mac_1 between BP and AD")


############Pso vs BP########################
#######load Ad and trim dataframe
Pso_Mac1 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/PsovsHC_DEG.xlsx", sheet = 5, colNames = T, rowNames = T)
Pso_Mac1 <- Pso_Mac1[,c(2,5)]
colnames(Pso_Mac1) <- c("Pso", "p_val_adj_1")

#######load BP and trim dataframe
BP_Mac1 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 5, colNames = T, rowNames = T)
BP_Mac1 <- BP_Mac1[,c(2,5)]
colnames(BP_Mac1) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(Pso_Mac1, BP_Mac1, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE Pso",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$Pso) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE Pso" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(Pso, BP, color=threshold))+
  geom_text_repel(aes(Pso, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("Mac_1 between BP and Pso")

################################# Mac_2###############geneX expression across the cell types########################
#######load Ad and trim dataframe
AD_Mac2 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/ADvsHC_DEG.xlsx", sheet = 6, colNames = T, rowNames = T)
AD_Mac2 <- AD_Mac2[,c(2,5)]
colnames(AD_Mac2) <- c("AD", "p_val_adj_1")

#######load BP and trim dataframe
BP_Mac2 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 6, colNames = T, rowNames = T)
BP_Mac2 <- BP_Mac2[,c(2,5)]
colnames(BP_Mac2) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(AD_Mac2, BP_Mac2, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE AD",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$AD) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE AD" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(AD, BP, color=threshold))+
  geom_text_repel(aes(AD, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("Mac_2 between BP and AD")


############Pso vs BP########################
#######load Ad and trim dataframe
Pso_Mac2 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/PsovsHC_DEG.xlsx", sheet = 6, colNames = T, rowNames = T)
Pso_Mac2 <- Pso_Mac2[,c(2,5)]
colnames(Pso_Mac2) <- c("Pso", "p_val_adj_1")

#######load BP and trim dataframe
BP_Mac2 <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 6, colNames = T, rowNames = T)
BP_Mac2 <- BP_Mac2[,c(2,5)]
colnames(BP_Mac2) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(Pso_Mac2, BP_Mac2, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE Pso",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$Pso) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE Pso" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(Pso, BP, color=threshold))+
  geom_text_repel(aes(Pso, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("Mac_2 between BP and Pso")

################################# moDC###############geneX expression across the cell types########################
#######load Ad and trim dataframe
AD_moDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/ADvsHC_DEG.xlsx", sheet = 1, colNames = T, rowNames = T)
AD_moDC <- AD_moDC[,c(2,5)]
colnames(AD_moDC) <- c("AD", "p_val_adj_1")

#######load BP and trim dataframe
BP_moDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 1, colNames = T, rowNames = T)
BP_moDC <- BP_moDC[,c(2,5)]
colnames(BP_moDC) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(AD_moDC, BP_moDC, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE AD",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$AD) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE AD" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(AD, BP, color=threshold))+
  geom_text_repel(aes(AD, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("moDC between BP and AD")


############Pso vs BP########################
#######load Ad and trim dataframe
Pso_moDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/PsovsHC_DEG.xlsx", sheet = 1, colNames = T, rowNames = T)
Pso_moDC <-Pso_moDC[,c(2,5)]
colnames(Pso_moDC) <- c("Pso", "p_val_adj_1")

#######load BP and trim dataframe
BP_moDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 1, colNames = T, rowNames = T)
BP_moDC <- BP_moDC[,c(2,5)]
colnames(BP_moDC) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(Pso_moDC, BP_moDC, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE Pso",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$Pso) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE Pso" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(Pso, BP, color=threshold))+
  geom_text_repel(aes(Pso, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("moDC between BP and Pso")

################################# migDC###############geneX expression across the cell types########################
#######load Ad and trim dataframe
AD_migDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/ADvsHC_DEG.xlsx", sheet = 3, colNames = T, rowNames = T)
AD_migDC <- AD_migDC[,c(2,5)]
colnames(AD_migDC) <- c("AD", "p_val_adj_1")

#######load BP and trim dataframe
BP_migDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 3, colNames = T, rowNames = T)
BP_migDC <- BP_migDC[,c(2,5)]
colnames(BP_migDC) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(AD_migDC, BP_migDC, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE AD",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$AD) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE AD" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(AD, BP, color=threshold))+
  geom_text_repel(aes(AD, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("migDC between BP and AD")


############Pso vs BP########################
#######load Ad and trim dataframe
Pso_migDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/PsovsHC_DEG.xlsx", sheet = 3, colNames = T, rowNames = T)
Pso_migDC <- Pso_migDC[,c(2,5)]
colnames(Pso_migDC) <- c("Pso", "p_val_adj_1")

#######load BP and trim dataframe
BP_migDC <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/BPvsHC_DEG.xlsx", sheet = 3, colNames = T, rowNames = T)
BP_migDC <- BP_migDC[,c(2,5)]
colnames(BP_migDC) <- c("BP", "p_val_adj_2")


Cluster0 <- merge(Pso_migDC, BP_migDC, by = "row.names")
#use values from first column as row names
#rownames(Cluster0) <- Cluster0[,1]

#remove first column from data frame
#Cluster0 <- Cluster0[,-1]

########Make scatterplot
Cluster0$threshold <- ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE both", 
                             ifelse(Cluster0$p_val_adj_1 < 0.05 & Cluster0$p_val_adj_2 > 0.05, "DE Pso",
                                    ifelse(Cluster0$p_val_adj_1 > 0.05 & Cluster0$p_val_adj_2 < 0.05, "DE BP", "NS")))
Cluster0$genelabels <- ""

Cluster0$genelabels <- ifelse(abs(Cluster0$Pso) > 5
                              |abs(Cluster0$BP) > 5, TRUE, FALSE)

cols <- c("DE both" = "purple", "DE Pso" ="blue", "DE BP" = "red",  "NS" = "grey")


ggplot(Cluster0, aes(Pso, BP, color=threshold))+
  geom_text_repel(aes(Pso, BP, label= ifelse(Cluster0$genelabels, as.character(Cluster0$Row.names), "")), color = "black", size = 5, vjust = 2, box.padding = unit(0.2, "lines"))+
  geom_point(size=3)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme() +  
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  ggtitle("migDC between BP and Pso")

#############################sorted cytokines#################################
gene <- c("CCL14",	"CCL15",	"CCL2",	"CCL22",	"CCL23",	"CCL27",	"CCL3",	"CXCL1",	"CXCL12",	"CXCL14",	"CXCL16",
          "CXCL3",	"CXCL5",	"CXCL8",	"IFNG",	"IL10",	"IL1B",	"IL20",	"IL23A",	"IL33",	"TNF",	"IL32",
          "CCL18",	"CCL3L1",	"IL6ST",	"IL2RG")

marker <- FindMarkers(APC, ident.1 = "BP", ident.2 = "NML", subset.ident = "Mac_1", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/BPvsNML_sort_cytokines.xlsx", colNames=T, rowNames=T)


marker <- FindMarkers(APC, ident.1 = "Pso", ident.2 = "NML", subset.ident = "Mac_1", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/PsovsNML_sort_cytokines.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "AD", ident.2 = "NML", subset.ident = "Mac_1", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/ADvsNML_sort_cytokines.xlsx", colNames=T, rowNames=T)

library(ComplexHeatmap)
library(circlize)
data <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 1, colNames=T, rowNames=T)
data <- as.matrix(data)
display <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 2, colNames=T, rowNames=T)
mycol<-colorRampPalette(c("blue","white","tomato"))(800)

data_mark=display
# 新建mark矩阵

for(i in 1:21){
  for(j in 1:3){
    if(display[i,j] <= 0.001)
    {
      data_mark[i,j]="***"
    }
    else if(display[i,j] <= 0.01 && display[i,j] > 0.001)
    {
      data_mark[i,j]="**"
    }
    else if(display[i,j] <= 0.05 && display[i,j] > 0.01)
    {
      data_mark[i,j]="*"
    }
    else
    {
      data_mark[i,j]=""
    }
  }
}


data_mark <- as.matrix(data_mark)######非常重要，否则热图上显示不出来

pheatmap(data,  
         color=mycol, 
         border_color = "black",  
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T, 
         annotation_row = annotation_row, 
         legend = TRUE, 
         legend_breaks = c(-1, 0, 1), 
         legend_labels = c("low","","heigh"), 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize = 8,
         display_numbers = data_mark,  #将预先定义的矩阵传递给display_numbers参数
         number_color = "black",
         fontsize_number = 10,#*大小
         fontsize_col=10, fontsize_row=10,
         fontface_row ="italic",
         main = "Mac_1")

#########################
#############################sorted cytokines#################################
gene <- c("IL1B",	"CCL3",	"CCL4",	"CCL5",	"CCR5",	"CXCL8",	"CXCR3",	"IL18",	"IFNG",	"IL12A",	"CXCL10",	"CXCL11",	"CXCL16",
          "CCR2", "IFNGR1",	"IFNGR2",	"CXCL9",	"CCR7",	"XCR1",	"XCL1",	"IL2",	"IL6",	"IL15",	"CSF2",	"TNF",	"STAT1", "CCL18",
          "CCL17",	"CCL13",	"CCL14",	"CTSC",	"CCL26",	"IL4",	"IL5",	"IL19",	"CCL11",	"CCR3", "CCR4",	"CCR8",	"IL13",	"TSLP",
          "IL25",	"IL33",	"GATA3",	"IL9",	"IL31", "IL17A",	"IL17F",	"IL22",	"CCL20",	"CCR6",	"CXCL2",	"CXCR1",	"CXCR2",	"RORC")

#############BP
marker <- FindMarkers(APC, ident.1 = "BP", ident.2 = "NML", subset.ident = "Mac_1", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/BPvsNML_Mac_1_type123.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "BP", ident.2 = "NML", subset.ident = "Mac_2", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/BPvsNML_Mac_2_type123.xlsx", colNames=T, rowNames=T)


marker <- FindMarkers(APC, ident.1 = "BP", ident.2 = "NML", subset.ident = "moDC", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/BPvsNML_moDC_type123.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "BP", ident.2 = "NML", subset.ident = "migDC", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/BPvsNML_migDC_type123.xlsx", colNames=T, rowNames=T)
#############Pso
marker <- FindMarkers(APC, ident.1 = "Pso", ident.2 = "NML", subset.ident = "Mac_1", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/PsovsNML_Mac_1_type123.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "Pso", ident.2 = "NML", subset.ident = "Mac_2", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/PsovsNML_Mac_2_type123.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "Pso", ident.2 = "NML", subset.ident = "moDC", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/PsovsNML_moDC_type123.xlsx", colNames=T, rowNames=T)


marker <- FindMarkers(APC, ident.1 = "Pso", ident.2 = "NML", subset.ident = "migDC", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/PsovsNML_migDC_type123.xlsx", colNames=T, rowNames=T)
#############AD
marker <- FindMarkers(APC, ident.1 = "AD", ident.2 = "NML", subset.ident = "Mac_1", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/ADvsNML_Mac_1_type123.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "AD", ident.2 = "NML", subset.ident = "Mac_2", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/ADvsNML_Mac_2_type123.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "AD", ident.2 = "NML", subset.ident = "moDC", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/ADvsNML_moDC_type123.xlsx", colNames=T, rowNames=T)

marker <- FindMarkers(APC, ident.1 = "AD", ident.2 = "NML", subset.ident = "migDC", group.by = "dis", verbose =TRUE,
                      assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                      logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = gene)
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/APC/ADvsNML_migDC_type123.xlsx", colNames=T, rowNames=T)


library(ComplexHeatmap)
library(circlize)

#################MaC_1
data <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 3, colNames=T, rowNames=T)
data <- as.matrix(data)
display <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 4, colNames=T, rowNames=T)
mycol<-colorRampPalette(c("blue","white","tomato"))(800)
#######
annotation_row = data.frame(
 Genes = c("IL1B",	"CCL3",	"CCL4",	"CCL5",	"CCR5",	"CXCL8",	"CXCR3",	"IL18",	"IFNG",	"IL12A",	"CXCL10",	"CXCL11",	"CXCL16",	"CCR2",
           "IFNGR1",	"IFNGR2",	"CXCL9",	"CCR7",	"XCR1",	"XCL1",	"IL2",	"IL6",	"IL15",	"CSF2",	"TNF",	"STAT1",
           "CCL18",	"CCL17",	"CCL13",	"CCL14",	"CTSC",	"CCL26",	"IL4",	"IL5",	"IL19",	"CCL11",	"CCR3",
           "CCR4",	"CCR8",	"IL13",	"TSLP",	"IL25",	"IL33",	"GATA3",	"IL9",	"IL31",
           "IL17A",	"IL17F",	"IL22",	"CCL20",	"CCR6",	"CXCL2",	"CXCR1",	"CXCR2",	"RORC"),
 Type = c(rep("type 1",26), rep("type 2",20),rep("type 3",9))
)
rownames(annotation_row) <- annotation_row$Genes
#define columns to exclude
cols <- names(annotation_row) %in% c('Genes')

#exclude points column
anno_row <- annotation_row[!cols]
#################

data_mark=display
# 新建mark矩阵

for(i in 1:55){
  for(j in 1:3){
    if(display[i,j] <= 0.001)
    {
      data_mark[i,j]="***"
    }
    else if(display[i,j] <= 0.01 && display[i,j] > 0.001)
    {
      data_mark[i,j]="**"
    }
    else if(display[i,j] <= 0.05 && display[i,j] > 0.01)
    {
      data_mark[i,j]="*"
    }
    else
    {
      data_mark[i,j]=""
    }
  }
}


data_mark <- as.matrix(data_mark)######非常重要，否则热图上显示不出来

p1 <- pheatmap(data,  
         color=mycol, 
         border_color = "black",  
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T, 
         annotation_row = anno_row,
         legend = TRUE, 
         legend_breaks = c(-1, 0, 1), 
         legend_labels = c("low","","heigh"), 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_number = 10,#*大小
         fontsize_col=10, 
         fontsize_row=10,
         fontface_row ="italic",
         display_numbers = data_mark,  #将预先定义的矩阵传递给display_numbers参数
         number_color = "black",
         main = "Mac_1")

#################MaC_2
data <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 5, colNames=T, rowNames=T)
data <- as.matrix(data)
display <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 6, colNames=T, rowNames=T)
mycol<-colorRampPalette(c("blue","white","tomato"))(800)

data_mark=display
# 新建mark矩阵

for(i in 1:55){
  for(j in 1:3){
    if(display[i,j] <= 0.001)
    {
      data_mark[i,j]="***"
    }
    else if(display[i,j] <= 0.01 && display[i,j] > 0.001)
    {
      data_mark[i,j]="**"
    }
    else if(display[i,j] <= 0.05 && display[i,j] > 0.01)
    {
      data_mark[i,j]="*"
    }
    else
    {
      data_mark[i,j]=""
    }
  }
}


data_mark <- as.matrix(data_mark)######非常重要，否则热图上显示不出来

p2 <- pheatmap(data,  
               color=mycol, 
               border_color = "black",  
               scale = "row", 
               cluster_rows = T, 
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               annotation_row = anno_row,
               fontsize = 8,
               display_numbers = data_mark,  #将预先定义的矩阵传递给display_numbers参数
               number_color = "black",
               fontsize_number = 10,#*大小
               fontsize_col=10, fontsize_row=10,
               fontface_row ="italic",
               main = "Mac_2")

#################moDC
data <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 7, colNames=T, rowNames=T)
data <- as.matrix(data)
display <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 8, colNames=T, rowNames=T)
mycol<-colorRampPalette(c("blue","white","tomato"))(800)

data_mark=display
# 新建mark矩阵

for(i in 1:55){
  for(j in 1:3){
    if(display[i,j] <= 0.001)
    {
      data_mark[i,j]="***"
    }
    else if(display[i,j] <= 0.01 && display[i,j] > 0.001)
    {
      data_mark[i,j]="**"
    }
    else if(display[i,j] <= 0.05 && display[i,j] > 0.01)
    {
      data_mark[i,j]="*"
    }
    else
    {
      data_mark[i,j]=""
    }
  }
}


data_mark <- as.matrix(data_mark)######非常重要，否则热图上显示不出来

p3 <- pheatmap(data,  
               color=mycol, 
               border_color = "black",  
               scale = "row", 
               cluster_rows = T, 
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE,
               annotation_row = anno_row,
               fontsize = 8,
               display_numbers = data_mark,  #将预先定义的矩阵传递给display_numbers参数
               number_color = "black",
               fontsize_number = 10,#*大小
               fontsize_col=10, fontsize_row=10,
               fontface_row ="italic",
               main = "moDC")

#################MigDC
data <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 9, colNames=T, rowNames=T)
data <- as.matrix(data)
display <- read.xlsx("H:/BP data analysis/Integrated with SI data/APC/Integrated.xlsx", sheet = 10, colNames=T, rowNames=T)
mycol<-colorRampPalette(c("blue","white","tomato"))(800)

data_mark=display
# 新建mark矩阵

for(i in 1:55){
  for(j in 1:3){
    if(display[i,j] <= 0.001)
    {
      data_mark[i,j]="***"
    }
    else if(display[i,j] <= 0.01 && display[i,j] > 0.001)
    {
      data_mark[i,j]="**"
    }
    else if(display[i,j] <= 0.05 && display[i,j] > 0.01)
    {
      data_mark[i,j]="*"
    }
    else
    {
      data_mark[i,j]=""
    }
  }
}


data_mark <- as.matrix(data_mark)######非常重要，否则热图上显示不出来

p4 <- pheatmap(data,  
               color=mycol, 
               border_color = "black",  
               scale = "row", 
               cluster_rows = T, 
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               annotation_row = anno_row,
               fontsize = 8,
               display_numbers = data_mark,  #将预先定义的矩阵传递给display_numbers参数
               number_color = "black",
               fontsize_number = 10, #*大小
               fontsize_col=10, fontsize_row=10,
               fontface_row = "italic",
               main = "migDC")



ComplexHeatmap::draw(p1 + p2 + p3 + p4, ht_gap = unit(0.5, "cm"))

plot_grid(p1, p2, p3, p4)

###########################Vlnplot
VlnPlot(APC, features = "LAMP3",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CCL17",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CCL19",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CCL22",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CCR7",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CD200",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CD80",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CD274",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CD1B",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "IL15",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "IL32",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

#####################
###########################Vlnplot
VlnPlot(APC, features = "LAMP3",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CCL17",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CCL19",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


VlnPlot(APC, features = "CCL22",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CCR7",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CD200",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CD80",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(APC, features = "CD274",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

########################
FeaturePlot(APC, features = "IFNG", min.cutoff = 0.1,  order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(),
        plot.title = element_text(color='black',family="Arial",size=13, face = "italic"),axis.ticks.length = unit(0.4,"lines"),  
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=13),
        axis.title.y=element_text(colour='black', size=13),
        axis.text=element_text(colour='black',size=13),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=13),
        legend.key=element_blank())

FeaturePlot(APC, features = "IL13", min.cutoff = 0.1,  order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(),
        plot.title = element_text(color='black',family="Arial",size=13, face = "italic"),axis.ticks.length = unit(0.4,"lines"),  
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=13),
        axis.title.y=element_text(colour='black', size=13),
        axis.text=element_text(colour='black',size=13),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=13),
        legend.key=element_blank())

FeaturePlot(APC, features = "IL17A", min.cutoff = 0.1,  order = F, pt.size = 0.5)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(),
        plot.title = element_text(color='black',family="Arial",size=13, face = "italic"),axis.ticks.length = unit(0.4,"lines"),  
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=13),
        axis.title.y=element_text(colour='black', size=13),
        axis.text=element_text(colour='black',size=13),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=13),
        legend.key=element_blank())
