Merge <- RenameIdents(Merge, `0` = "Trm", `1` = "Treg", `2` = "Tcm", `3` = "CTLac", `4` = "CTLex_1", `5` = "Tmm", `6` = "Trm",
                      `7` = "Tcm", `8` = "Tet", `9` = "moDC", `10` = "NK", `11` = "Mast", `12` = "Mono", 
                      `13` = "migDC", `14` = "ILC", `15` = "CTLac", `16` = "CTLex_2",  `17` = "Baso", 
                      `18` = "LC", `19` = "moDC", `20` = "Mac_1", `21` = "Cycling_1", `22` = "B",
                      `23` = "Mac_2", `24` = "Mono", `25` = "Cycling_2")

#Tsubset <- subset(MergeF, idents=c( "Trm", "Treg", "Tcm", "CTLac", "CTLex_1", "Tmm", "Trm", "Tcm", "Tet", "NK", "ILC", "CTLac", "CTLex_2"))

#######################Figure 2C
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


gene_sig <- c("IL17A", "IL17F", "IL21", "IL22", "IL23A", "IL26", "IL36A", "CXCL13", "IL4", "IL5", "IL9", "IL13")

comparisons <- list(c("NML", "BP"), c("NML", "AD"), c("NML", "Pso"), c("BP", "Pso"), c("BP", "AD"), c("AD", "Pso"))
vp_case1(gene_signature = gene_sig, file_name = "H:/BP data analysis/Integrated with SI data/Type2_17", test_sign = comparisons)

#############################DGE for BP vs NML for Trm clusters
marker <- FindMarkers(MergeF, ident.1 = "BP", ident.2 = "NML", group.by = "dis", subset.ident = "Trm")
write.xlsx(marker, "H:/BP data analysis/Integrated with SI data/BPvsNML_Trm_DEG.xlsx", colNames=T, rowNames=T)
