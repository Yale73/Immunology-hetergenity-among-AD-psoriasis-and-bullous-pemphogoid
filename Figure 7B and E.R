##################################4B############################################
genes <- c('CD207','CD1A','FCGBP','HPGDS','CLEC9A','IRF8', 'XCR1','WDFY4','LAMP3','CCR7','CD1B','BIRC3',
           'CLEC10A','IL1B','CXCL8','CXCL2','FCGR3A','HES1','CXCL9','CXCL10','CD14','CD163','C1QC','RNASE1','SELENOP')
DotPlot(APC, features = genes, cols=c("#DDA0DD", "#6A0DAD", "#502380", "#290916"), assay = "RNA", col.min = 0.1, col.max = 1, dot.min=0.05, dot.scale = 1,
        cluster.idents=T, group.by = "celltype")+
  scale_size(range = c(0, 5))+ 
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, family="TT Times New Roman", face = "italic"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 9, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) +
  labs(x="", y="")

##################################4D---BPvsHC############################################
# First, you need to analyse the DEG for each individual cluster. Below is the example for one cluster, but you need to apply it to all your clusters...
######c0
C0 <- subset(APC, idents = "moDC")

Idents(C0) <- C0@meta.data$dis
DefaultAssay(C0) <- "RNA"

DEG.Data.C0 <- FindMarkers(C0, ident.1 = "BP", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C0$Gene <- rownames(DEG.Data.C0)

######C1
C1 <- subset(APC, idents = "Mac_2")

Idents(C1) <- C1@meta.data$dis
DefaultAssay(C1) <- "RNA"

DEG.Data.C1 <- FindMarkers(C1, ident.1 = "BP", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C1$Gene <- rownames(DEG.Data.C1)

######C2
C2 <- subset(APC, idents = "migDC")

Idents(C2) <- C2@meta.data$dis
DefaultAssay(C2) <- "RNA"

DEG.Data.C2 <- FindMarkers(C2, ident.1 = "BP", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C2$Gene <- rownames(DEG.Data.C2)

######C3
C3 <- subset(APC, idents = "LC")

Idents(C3) <- C3@meta.data$dis
DefaultAssay(C3) <- "RNA"

DEG.Data.C3 <- FindMarkers(C3, ident.1 = "BP", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C3$Gene <- rownames(DEG.Data.C3)

######C4
C4 <- subset(APC, idents = "Mac_1")

Idents(C4) <- C4@meta.data$dis
DefaultAssay(C4) <- "RNA"

DEG.Data.C4 <- FindMarkers(C4, ident.1 = "BP", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C4$Gene <- rownames(DEG.Data.C4)

######C5
C5 <- subset(APC, idents = "Mono")

Idents(C5) <- C5@meta.data$dis
DefaultAssay(C5) <- "RNA"

DEG.Data.C5 <- FindMarkers(C5, ident.1 = "BP", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C5$Gene <- rownames(DEG.Data.C5)



#stores each cluster and each cell in an object
DEG.Data <- FetchData(APC, "celltype")


# Once its done, you run the following (just adapt it to the number of clusters you have). In this case we decided to go with DEG selection criteria of LogFC>0.25 and adj_Pvalue <0.05.
# But you can change these as well. We had 25 clusters in this example, and clusters 22, 23, 24 and 25 had 0 DEG that we entered manually as you can see at the end of this code. 
#Again, just adapt to the number of clusters you have and pay attention at the number of parenthesis to close at the end of the line. 


DEG.Data$nDEG <- ifelse(DEG.Data$celltype=="moDC",  sum(abs(DEG.Data.C0$avg_log2FC)>0.25 & DEG.Data.C0$p_val_adj <0.05),
                        ifelse(DEG.Data$celltype=="Mac_2", sum(abs(DEG.Data.C1$avg_log2FC)>0.25 & DEG.Data.C1$p_val_adj <0.05),
                               ifelse(DEG.Data$celltype=="migDC", sum(abs(DEG.Data.C2$avg_log2FC)>0.25 & DEG.Data.C2$p_val_adj <0.05),
                                      ifelse(DEG.Data$celltype=="LC", sum(abs(DEG.Data.C3$avg_log2FC)>0.25 & DEG.Data.C3$p_val_adj <0.05),
                                             ifelse(DEG.Data$celltype=="Mac_1", sum(abs(DEG.Data.C4$avg_log2FC)>0.25 & DEG.Data.C4$p_val_adj <0.05),
                                                    ifelse(DEG.Data$celltype=="Mono", sum(abs(DEG.Data.C5$avg_log2FC)>0.25 & DEG.Data.C5$p_val_adj <0.05),
                                                           0))))))
#Statch the number of DEG for each cluster into the meta.data slot of your object (here under the name of nDEG)

APC@meta.data$nDEG <- DEG.Data$nDEG


#Project the result on the UMAP plot

minimal<- min(DEG.Data$nDEG)#109
maximal<- max(DEG.Data$nDEG)#6241

P1 <- FeaturePlot(APC, features = "nDEG", pt.size = 0.5)+ggtitle("nDEG of BPvsHC")+
  scale_color_gradient2(low = 'grey', high = 'red',limits=c(109, 6241))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=12),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

##################################4D---ADvsHC############################################
# First, you need to analyse the DEG for each individual cluster. Below is the example for one cluster, but you need to apply it to all your clusters...
######c0
C0 <- subset(APC, idents = "moDC")

Idents(C0) <- C0@meta.data$dis
DefaultAssay(C0) <- "RNA"

DEG.Data.C0 <- FindMarkers(C0, ident.1 = "AD", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C0$Gene <- rownames(DEG.Data.C0)

######C1
C1 <- subset(APC, idents = "Mac_2")

Idents(C1) <- C1@meta.data$dis
DefaultAssay(C1) <- "RNA"

DEG.Data.C1 <- FindMarkers(C1, ident.1 = "AD", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C1$Gene <- rownames(DEG.Data.C1)

######C2
C2 <- subset(APC, idents = "migDC")

Idents(C2) <- C2@meta.data$dis
DefaultAssay(C2) <- "RNA"

DEG.Data.C2 <- FindMarkers(C2, ident.1 = "AD", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C2$Gene <- rownames(DEG.Data.C2)

######C3
C3 <- subset(APC, idents = "LC")

Idents(C3) <- C3@meta.data$dis
DefaultAssay(C3) <- "RNA"

DEG.Data.C3 <- FindMarkers(C3, ident.1 = "AD", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C3$Gene <- rownames(DEG.Data.C3)

######C4
C4 <- subset(APC, idents = "Mac_1")

Idents(C4) <- C4@meta.data$dis
DefaultAssay(C4) <- "RNA"

DEG.Data.C4 <- FindMarkers(C4, ident.1 = "AD", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C4$Gene <- rownames(DEG.Data.C4)

######C5
C5 <- subset(APC, idents = "Mono")

Idents(C5) <- C5@meta.data$dis
DefaultAssay(C5) <- "RNA"

DEG.Data.C5 <- FindMarkers(C5, ident.1 = "AD", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C5$Gene <- rownames(DEG.Data.C5)



#stores each cluster and each cell in an object
DEG.Data <- FetchData(APC, "celltype")


# Once its done, you run the following (just adapt it to the number of clusters you have). In this case we decided to go with DEG selection criteria of LogFC>0.25 and adj_Pvalue <0.05.
# But you can change these as well. We had 25 clusters in this example, and clusters 22, 23, 24 and 25 had 0 DEG that we entered manually as you can see at the end of this code. 
#Again, just adapt to the number of clusters you have and pay attention at the number of parenthesis to close at the end of the line. 


DEG.Data$nDEG <- ifelse(DEG.Data$celltype=="moDC",  sum(abs(DEG.Data.C0$avg_log2FC)>0.25 & DEG.Data.C0$p_val_adj <0.05),
                        ifelse(DEG.Data$celltype=="Mac_2", sum(abs(DEG.Data.C1$avg_log2FC)>0.25 & DEG.Data.C1$p_val_adj <0.05),
                               ifelse(DEG.Data$celltype=="migDC", sum(abs(DEG.Data.C2$avg_log2FC)>0.25 & DEG.Data.C2$p_val_adj <0.05),
                                      ifelse(DEG.Data$celltype=="LC", sum(abs(DEG.Data.C3$avg_log2FC)>0.25 & DEG.Data.C3$p_val_adj <0.05),
                                             ifelse(DEG.Data$celltype=="Mac_1", sum(abs(DEG.Data.C4$avg_log2FC)>0.25 & DEG.Data.C4$p_val_adj <0.05),
                                                    ifelse(DEG.Data$celltype=="Mono", sum(abs(DEG.Data.C5$avg_log2FC)>0.25 & DEG.Data.C5$p_val_adj <0.05),
                                                           0))))))
#Statch the number of DEG for each cluster into the meta.data slot of your object (here under the name of nDEG)

APC@meta.data$nDEG <- DEG.Data$nDEG


#Project the result on the UMAP plot

P2 <- FeaturePlot(APC, features = "nDEG", pt.size = 0.5)+ggtitle("nDEG of ADvsHC")+
  scale_color_gradient2(low = 'grey', high = 'red',limits=c(109, 6241))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=12),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

##################################4D---PsovsHC############################################
# First, you need to analyse the DEG for each individual cluster. Below is the example for one cluster, but you need to apply it to all your clusters...
######c0
C0 <- subset(APC, idents = "moDC")

Idents(C0) <- C0@meta.data$dis
DefaultAssay(C0) <- "RNA"

DEG.Data.C0 <- FindMarkers(C0, ident.1 = "Pso", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C0$Gene <- rownames(DEG.Data.C0)

######C1
C1 <- subset(APC, idents = "Mac_2")

Idents(C1) <- C1@meta.data$dis
DefaultAssay(C1) <- "RNA"

DEG.Data.C1 <- FindMarkers(C1, ident.1 = "Pso", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C1$Gene <- rownames(DEG.Data.C1)

######C2
C2 <- subset(APC, idents = "migDC")

Idents(C2) <- C2@meta.data$dis
DefaultAssay(C2) <- "RNA"

DEG.Data.C2 <- FindMarkers(C2, ident.1 = "Pso", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C2$Gene <- rownames(DEG.Data.C2)

######C3
C3 <- subset(APC, idents = "LC")

Idents(C3) <- C3@meta.data$dis
DefaultAssay(C3) <- "RNA"

DEG.Data.C3 <- FindMarkers(C3, ident.1 = "Pso", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C3$Gene <- rownames(DEG.Data.C3)

######C4
C4 <- subset(APC, idents = "Mac_1")

Idents(C4) <- C4@meta.data$dis
DefaultAssay(C4) <- "RNA"

DEG.Data.C4 <- FindMarkers(C4, ident.1 = "Pso", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C4$Gene <- rownames(DEG.Data.C4)

######C5
C5 <- subset(APC, idents = "Mono")

Idents(C5) <- C5@meta.data$dis
DefaultAssay(C5) <- "RNA"

DEG.Data.C5 <- FindMarkers(C5, ident.1 = "Pso", ident.2 = "NML",logfc.threshold = 0.25)
DEG.Data.C5$Gene <- rownames(DEG.Data.C5)



#stores each cluster and each cell in an object
DEG.Data <- FetchData(APC, "celltype")


# Once its done, you run the following (just adapt it to the number of clusters you have). In this case we decided to go with DEG selection criteria of LogFC>0.25 and adj_Pvalue <0.05.
# But you can change these as well. We had 25 clusters in this example, and clusters 22, 23, 24 and 25 had 0 DEG that we entered manually as you can see at the end of this code. 
#Again, just adapt to the number of clusters you have and pay attention at the number of parenthesis to close at the end of the line. 


DEG.Data$nDEG <- ifelse(DEG.Data$celltype=="moDC",  sum(abs(DEG.Data.C0$avg_log2FC)>0.25 & DEG.Data.C0$p_val_adj <0.05),
                        ifelse(DEG.Data$celltype=="Mac_2", sum(abs(DEG.Data.C1$avg_log2FC)>0.25 & DEG.Data.C1$p_val_adj <0.05),
                               ifelse(DEG.Data$celltype=="migDC", sum(abs(DEG.Data.C2$avg_log2FC)>0.25 & DEG.Data.C2$p_val_adj <0.05),
                                      ifelse(DEG.Data$celltype=="LC", sum(abs(DEG.Data.C3$avg_log2FC)>0.25 & DEG.Data.C3$p_val_adj <0.05),
                                             ifelse(DEG.Data$celltype=="Mac_1", sum(abs(DEG.Data.C4$avg_log2FC)>0.25 & DEG.Data.C4$p_val_adj <0.05),
                                                    ifelse(DEG.Data$celltype=="Mono", sum(abs(DEG.Data.C5$avg_log2FC)>0.25 & DEG.Data.C5$p_val_adj <0.05),
                                                           0))))))
#Statch the number of DEG for each cluster into the meta.data slot of your object (here under the name of nDEG)

APC@meta.data$nDEG <- DEG.Data$nDEG


#Project the result on the UMAP plot

P3 <- FeaturePlot(APC, features = "nDEG", pt.size = 0.5)+ggtitle("nDEG of PsovsHC")+
  scale_color_gradient2(low = 'grey', high = 'red',limits=c(109, 6241))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=12),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

wrap_plots(P1, P2, P3, ncol = 1)
