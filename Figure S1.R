Merge <- readRDS("H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")

######################################Figure S1A################################
DimPlot(Merge, reduction = "umap",  pt.size = .6, group.by="celltype", label = T, repel = F, split.by = "source") + NoLegend()+
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


######################################Figure S1B################################
Idents(Merge) <- Merge$source
levels(Merge$source)
Liu <- subset(Merge, idents="Liu et al")
Idents(Liu) <- Liu$Ident1

DimPlot(Liu, reduction = "umap",  pt.size = .6, group.by="Ident1", label = T, repel = T) + NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(23))+
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
  guides(colour = guide_legend(override.aes = list(size=5)))+
  ggtitle("Idents from Liu et al")

################################Figure S1C
Idents(Merge) <- Merge$celltype
levels(Merge$celltype)


Idents(Merge) <- Merge$source
levels(Merge$source)
Liu <- subset(Merge, idents="Liu et al")


data <- Liu@meta.data[, c("celltype", "Ident1")]
colnames(data) <- c("Manuscript", "Liu et al")

#load psych package
library(psych)

#create pairs plot
pairs.panels(data)

################################Figure S1E
DimPlot(Merge, reduction = "umap",  pt.size = 1, group.by="celltype", split.by = "orig.ident", label = F, repel = F, ncol = 7) +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(21))+
  theme(text=element_text(family="Arial",size=12)) +
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
        legend.key=element_blank(),
        legend.position="bottom")+
  theme(plot.title = element_text(size=15,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5), nrow = 2))
