gene  <- list(Type_1=c("IL1B",	"CCL3",	"CCL4",	"CCL5",	"CCR5",	"CXCL8",	"CXCR3",	"IL18",	"IFNG",	"IL12A",	"CXCL10",	"CXCL11",	"CXCL16",	"CCR2",
                                                   "IFNGR1",	"IFNGR2",	"CXCL9",	"CCR7",	"XCR1",	"XCL1",	"IL2",	"IL6",	"IL15",	"CSF2",	"TNF",	"STAT1"),
                     Type_2=c("CCL18",	"CCL17",	"CCL13",	"CCL14",	"CTSC",	"CCL26",	"IL4",	"IL5",	"IL19",	"CCL11",	"CCR3",
                                                   "CCR4",	"CCR8",	"IL13",	"TSLP",	"IL25",	"IL33",	"GATA3",	"IL9",	"IL31"),
                    Type_3=c("IL17A",	"IL17F",	"IL22",	"CCL20",	"CCR6",	"CXCL2",	"CXCR1",	"CXCR2",	"RORC")
)

Trm <- subset(Merge, idents="Trm")
Trm$dis <- factor(Trm$dis, levels = c("HC", "AD", "Pso", "BP"))
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


DotPlot(Trm, features = gene, split.by = "dis", cols = c(RColorBrewer::brewer.pal(9, "Oranges")[7], 
                                                                RColorBrewer::brewer.pal(9, "Blues")[7], 
                                                                RColorBrewer::brewer.pal(9, "PuRd")[8],
                                                                RColorBrewer::brewer.pal(9, "Greens")[8]), 
        assay = "RNA", col.min = 0.1, col.max = 1, dot.min = 0.05, cluster.idents=F)+
  scale_size(range = c(0, 10))+
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, family="TT Times New Roman", face = "italic"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 9, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) +
  labs(x="", y="")
#+scale_color_gradientn(colours = viridis::magma(20), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)')

DotPlot(Trm, features = gene, split.by = "dis", group.by = "celltype", cols = c(RColorBrewer::brewer.pal(9, "Oranges")[7], 
                                                         RColorBrewer::brewer.pal(9, "Blues")[7], 
                                                         RColorBrewer::brewer.pal(9, "PuRd")[8],
                                                         RColorBrewer::brewer.pal(9, "Greens")[8]), 
        assay = "RNA", col.min = 0, col.max = 1, dot.min =0.01, cluster.idents=F)+
  scale_size(range = c(0, 10))+
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, family="TT Times New Roman", face = "italic"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 9, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) +
  labs(x="", y="") #+coord_flip()


#######No legend
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 10),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0, face = "italic"))
}
Idents(Trm) <- Trm$dis
my_comparisons <- list(c("BP", "HC"), c("AD", "HC"), c("Pso", "HC"), c("AD", "BP"), c("AD", "Pso"),c("BP", "Pso")) 
#####################Figure 4E
library(ggpubr)
VlnPlot(Trm, features = 'IL17A', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'IL17F', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'IL4', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'IL13', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较


VlnPlot(Trm, features = 'IFNG', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较


VlnPlot(Trm, features = 'TNF', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较
#####################Figure S4
VlnPlot(Trm, features = 'IL4R', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'IL13RA1', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'IL13RA2', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'IL2RG', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'IL17RA', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

VlnPlot(Trm, features = 'TNFRSF1A', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较
