library(Seurat)
library(ggplot2)
library(RColorBrewer)
MergeF <- readRDS("H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")

Trm <- subset(MergeF, idents="Trm")
DimPlot(Trm, reduction = "umap",  pt.size = .6, group.by="Ident", label = T, repel = F) + NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(21))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=15,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


Idents(Trm) <- Trm$orig.ident
table(Idents(Trm))
############Figure3A________CD3 and cytokines#############
CD3 <- subset(Trm, CD3D>0|CD3G>0)
table(Idents(CD3))

IL4 <- subset(CD3, IL4>0)
table(Idents(IL4))

IL5 <- subset(CD3, IL5>0)
table(Idents(IL5))


IL9 <- subset(CD3, IL9>0)
table(Idents(IL9))


IL13<- subset(CD3, IL13>0)
table(Idents(IL13))


IL17A<- subset(CD3, IL17A>0)
table(Idents(IL17A))


IL17F<- subset(CD3, IL17F>0)
table(Idents(IL17F))


IL21 <- subset(CD3, IL21>0)
table(Idents(IL21))



IL22 <- subset(CD3, IL22>0)
table(Idents(IL22))


IL23A <- subset(CD3, IL23A>0)
table(Idents(IL23A))


IL26 <- subset(CD3, IL26>0)
table(Idents(IL26))


CXCL13<- subset(CD3, CXCL13>0)
table(Idents(CXCL13))

IFNG <- subset(CD3, IFNG>0)
table(Idents(IFNG))

IL12A <- subset(CD3, IL12A>0)
table(Idents(IL12A))

TNF <- subset(CD3, TNF>0)
table(Idents(TNF))
############Figure3B________CD3CD4 and cytokines#############
CD3CD4 <- subset(Trm, (CD3D>0|CD3G>0)&CD4>0)
table(Idents(CD3CD4))


IL4 <- subset(CD3CD4, IL4>0)
table(Idents(IL4))

IL5 <- subset(CD3CD4, IL5>0)
table(Idents(IL5))


IL9 <- subset(CD3CD4, IL9>0)
table(Idents(IL9))


IL13<- subset(CD3CD4, IL13>0)
table(Idents(IL13))


IL17A<- subset(CD3CD4, IL17A>0)
table(Idents(IL17A))


IL17F<- subset(CD3CD4, IL17F>0)
table(Idents(IL17F))


IL21 <- subset(CD3CD4, IL21>0)
table(Idents(IL21))



IL22 <- subset(CD3CD4, IL22>0)
table(Idents(IL22))


IL23A <- subset(CD3CD4, IL23A>0)
table(Idents(IL23A))


IL26 <- subset(CD3CD4, IL26>0)
table(Idents(IL26))


CXCL13<- subset(CD3CD4, CXCL13>0)
table(Idents(CXCL13))


IFNG <- subset(CD3CD4, IFNG>0)
table(Idents(IFNG))


IL12A <- subset(CD3CD4, IL12A>0)
table(Idents(IL12A))

TNF <- subset(CD3CD4, TNF>0)
table(Idents(TNF))

############Figure3C________CD3CD4CD8 and cytokines#############
CD3CD8 <- subset(Trm, (CD3D>0|CD3G>0)& (CD8A>0|CD8B>0))
table(Idents(CD3CD8))

IL4 <- subset(CD3CD8, IL4>0)
table(Idents(IL4))

IL5 <- subset(CD3CD8, IL5>0)
table(Idents(IL5))


IL9 <- subset(CD3CD8, IL9>0)
table(Idents(IL9))


IL13<- subset(CD3CD8, IL13>0)
table(Idents(IL13))


IL17A<- subset(CD3CD8, IL17A>0)
table(Idents(IL17A))


IL17F<- subset(CD3CD8, IL17F>0)
table(Idents(IL17F))


IL21 <- subset(CD3CD8, IL21>0)
table(Idents(IL21))



IL22 <- subset(CD3CD8, IL22>0)
table(Idents(IL22))


IL23A <- subset(CD3CD8, IL23A>0)
table(Idents(IL23A))


IL26 <- subset(CD3CD8, IL26>0)
table(Idents(IL26))


CXCL13<- subset(CD3CD8, CXCL13>0)
table(Idents(CXCL13))


IFNG <- subset(CD3CD8, IFNG>0)
table(Idents(IFNG))


IL12A <- subset(CD3CD8, IL12A>0)
table(Idents(IL12A))

TNF <- subset(CD3CD8, TNF>0)
table(Idents(TNF))
############################Figure 4D
genes <- c("IL4", "IL5", "IL13", "IL17A", "IL17F", "CXCL13", "IFNG", "IL9", "IL22")
genes <- c('IFNG', 'TBX21', 'TNF', 'STAT4', 'CXCR3','IL4', 'IL5', 'IL13', 'GATA3', 'STAT6', 'CCR4','IL17A', 'IL17F', 'IL21', 'IL22', "IL23A", 'RORC', 'STAT3', 'CCR6','CXCL13',
           'BCL6', 'ICOS', 'PDCD1', 'CXCR5','IL9', 'SPI1', 'IRF4', 'CCR10', 'AHR')




Idents(Trm) <- Trm$Ident

set1 <- FindMarkers(Trm, ident.1 = "BP", ident.2 ="NML", group.by = "dis", verbose =TRUE,
                                       assay = "RNA", slot = "data", test.use = "MAST", min.cells.feature = 0, min.cells.group = 0, 
                                       subset.ident = "Trm", logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, features = genes)

write.xlsx(set1, file = paste0("H:/BP data analysis/Integrated with SI data/80 percentage DEG/", "BP", "-", "Cytokine-TF", ".xlsx"), colNames=T, rowNames=T)


######### Figure 3D
F3D <- read.xlsx("H:/BP data analysis/Integrated with SI data/80 percentage DEG/Integrated cytokines.xlsx", sheet = 1, colNames = T, rowNames = T)
library(ComplexHeatmap)
library(circlize)
Heatmap(as.matrix(F3D),
        name = "log2FC", #### annotate legend
        col = colorRamp2(c(-2,0, 2), c("blue", "white", "red")), #### set the color scales
        row_names_gp = gpar(fontsize = 12, fontface="italic"),
        column_names_gp = gpar(fontsize = 12),
        clustering_distance_columns = "canberra",
        column_title = "Cytokines and transcription factors",
        show_column_names = T,
        cluster_columns = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12),
        row_km = 3,
        column_km = 2)
