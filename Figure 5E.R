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

Idents(Trm) <- Trm$dis
## create fake 'Reponse' variable
library(ggpubr)
Trm@meta.data$Response <- Trm@meta.data$dis

vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(Trm, features = signature,
            pt.size = 0.1, 
            group.by = "Response", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 3) )#+3 to get all the hidden p-val
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 16, height = 14)
}

gene_sig <- c('NPDC1', 'UCP2', 'LIME1', 'PLAAT4', 'ACTB', 'ACTG1', 'HSPA1B', 'ARHGAP4', 'RAC2', 'S100A4', 'RPS4Y1', 'ITGB2', 'LDHB',
              'HCST', 'C9orf16', 'LIMD2', 'ATP5MC2', 'CORO1A', 'HSPA1A', 'IFITM3', 'SEPTIN1', 'ZAP70', 'ARHGDIB', 'TMSB10', 'HSPA8', 'PFN1', 'CD74', 'CLIC1', 'ARPC1B', 'ATP5MG', 'CAPZB', 'RBM3', 'CRIP1', 'LAPTM5', 'ATP5F1E', 'COMMD6', 'TAF10', 'HIGD2A','MYL12A', 'RPL39', 'CD3E')

comparisons <- list(c("NML", "PV"), c("NML", "AD"), c("NML", "BP"), c("AD", "BP"), c("AD", "Pso"), c("BP", "Pso"))
vp_case1(gene_signature = gene_sig, file_name = "H:/BP data analysis/Integrated with SI data/80 percentage DEG/BP_specific_genes", test_sign = comparisons)


############################Figure 3E
#########with legend
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

#######No legend
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 10),
          axis.title.x = element_text(size = 0),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 10, vjust = 1, hjust = 0, face = "italic"))
}
#########
Trm$dis <- factor(Trm$dis, levels = c("HC", "AD", "Pso", "BP"))

VlnPlot(Trm, features = 'CORO1A', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+NoLegend()

VlnPlot(Trm, features = 'HCST', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+NoLegend()

VlnPlot(Trm, features = 'LAPTM5', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+NoLegend()

VlnPlot(Trm, features = 'RAC2', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+NoLegend()

VlnPlot(Trm, features = 'UCP2', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+NoLegend()

VlnPlot(Trm, features = 'LIME1', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+NoLegend()



gene_sig <- c('NPDC1', 'UCP2', 'LIME1', 'PLAAT4', 'ACTB', 'ACTG1', 'HSPA1B', 'ARHGAP4',
              'RAC2', 'S100A4', 'RPS4Y1', 'ITGB2', 'LDHB', 'HCST', 'C9orf16', 'LIMD2',
              'ATP5MC2', 'CORO1A', 'HSPA1A', 'IFITM3', 'SEPTIN1', 'ZAP70', 'ARHGDIB',
              'TMSB10', 'HSPA8', 'PFN1', 'CD74', 'CLIC1', 'ARPC1B', 'ATP5MG', 'CAPZB',
              'RBM3', 'CRIP1', 'LAPTM5', 'ATP5F1E', 'COMMD6', 'TAF10', 'HIGD2A','MYL12A',
              'RPL39', 'CD3E')


VlnPlot(Trm, features = 'UCP2', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+
  NoLegend()


VlnPlot(Trm, features = 'LIME1', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+
  NoLegend()


VlnPlot(Trm, features = 'HSPA1B', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+
  NoLegend()

VlnPlot(Trm, features = 'IFITM3', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+
  NoLegend()

VlnPlot(Trm, features = 'SEPTIN1', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+
  NoLegend()

VlnPlot(Trm, features = 'CRIP1', split.by = "dis", ncol=1) +
  stat_summary(fun.y = median, geom='crossbar', width=0.2, colour = "black")+
  geom_boxplot(aes(colour = "dis"), width = 0.2)+  # 添加箱线图图层
  theme_niwot()+
  NoLegend()

