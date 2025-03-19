Idents(Merge) <- Merge$orig.ident
Merge <- RenameIdents(Merge,`Skin150` = "NML", `Skin154` = "NML", `Skin155` = "NML",`Skin169` = "NML",`Skin195` = "NML", `Skin204` = "NML", 
                      `Skin207` = "NML", `Skin170` = "AD",`Skin198` = "AD", `Skin230` = "AD", `Skin231` = "AD", `Skin232` = "AD", 
                      `Skin233` = "AD", `Skin236` = "AD", `Skin165` = "Pso", `Skin173` = "Pso", `Skin194` = "Pso", `Skin199` = "Pso", 
                      `Skin211` = "Pso", `Skin222` = "Pso", `Skin234` = "Pso", `Skin235` = "Pso", `Skin167` = "RashX", `Skin174` = "RashX", 
                      `Skin192` = "RashX", `Skin200` = "RashX", `Skin202` = "RashX",  `Skin219` = "RashX", `Skin141` = "RashX", `Skin163` = "RashX", 
                      `Skin175` = "BP", `Skin203` = "RashX", `HC1` = "NML",`HC2` = "NML",`HC3` = "NML",`HC4` = "NML",`HC5` = "NML", 
                      `BP1` = "BP",`BP2` = "BP",`BP3` = "BP",`BP4` = "BP",`BP5` = "BP",`BP6` = "BP",`BP7` = "BP")
Merge[["dis"]] <- Idents(object = Merge)

Idents(Merge) <- Merge$orig.ident
levels(Merge$orig.ident)
Merge <- RenameIdents(Merge,`NML1` = "HC1", `NML2` = "HC2", `NML3` = "HC3",`NML4` = "HC4",`NML5` = "HC5", `NML6` = "HC6", 
                      `NML7` = "HC7", `NML8` = "HC8", `NML9` = "HC9",`NML10` = "HC10",`NML11` = "HC11", `NML12` = "HC12", 
                      `AD1` = "AD1", `AD2` = "AD2", `AD3` = "AD3", `AD4` = "AD4", `AD5` = "AD5", `AD6` = "AD6", `AD7` = "AD7", 
                      `Pso1` = "Pso1", `Pso2` = "Pso2", `Pso3` = "Pso3", `Pso4` = "Pso4", `Pso5` = "Pso5", `Pso6` = "Pso6", `Pso7` = "Pso7", `Pso8` = "Pso8", 
                      `BP1` = "BP1",`BP2` = "BP2",`BP3` = "BP3",`BP4` = "BP4",`BP5` = "BP5",`BP6` = "BP6",`BP7` = "BP7",`BP8` = "BP8")
Merge[["orig.ident"]] <- Idents(object = Merge)


Idents(Merge) <- Merge$dis
Merge <- subset(Merge, idents=c("NML", "AD", "Pso", "BP"))
Idents(Merge) <- Merge$Ident
saveRDS(Merge, "H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")


############
Merge <- readRDS("H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")

Idents(Merge) <- Merge$group

Merge <- RenameIdents(Merge, `NML` = "HC", `AD` = "AD", `Pso` = "Pso", `BP` = "BP")

Merge[["group"]] <- Idents(object = Merge)
saveRDS(Merge, "H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")
################change source name
########Figure 1B
DimPlot(Merge, reduction = "umap",  pt.size = .6, group.by="celltype", label = T, repel = F) + NoLegend()+
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

########Figure 1C
Merge$dis <- factor(Merge$dis, levels = c("HC", "AD", "Pso", "BP"))
DimPlot(Merge, reduction = "umap",  pt.size = .6, group.by= "dis", label = F, repel = T, split.by = "dis", cols = c("#2c91e0", "#ff9a9b", "#c99bff", "#3abf99"), ncol = 2) + NoLegend()+
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
########Figure 1D
rashx_markers <- c("CD3D",
                   "CCR7", "SELL", "KLF2","CD69","TCF7", "S1PR1", #Tn/Tcm/Tmm
                   "ITGAE","CXCR6", #Trm
                   "CD4","FOXP3", "TIGIT", "CTLA4", "IL2RA", #Treg
                   "CD8A","CD8B","GZMB", "NKG7", "CCL5", #CTL
                   "KLRB1", #NK/ILC
                   "GZMA", "GZMK", "KLRC1", "EOMES", "PRF1","KLRD1","GNLY", "CCR8", #NK
                   "IKZF3", "TBX21", #ILC1
                   "GATA3","MAF", "PTGDR2", "HPGDS", #ILC2
                   "RORC", "IL23R", "IL1R1", "KIT", "AHR", #ILC3
                   "TNFRSF4", "CD96","TNFRSF18","PRDM1","PDCD1","LAG3", #ac or ex
                   "BATF", "SNHG12", "ZFAS1", #Tet
                   "TRAT1","RORA","IL7R",  #not specific
                   "HLA-DRA","HLA-DRB1","CD83","IDO1", #APC
                   "CD207","EPCAM", #LC
                   "CD68","CEBPB","FCER1G","C1QB","C1QC","CD163", "KLF4","CD14","S100A9", #Mac
                   "MS4A7","LYZ","SERPINA1", "IL23A","IL1B", "CXCL3", #Mono
                   "CD1C", "CLEC10A","CLEC9A", "XCR1", "FSCN1", "LAMP3", #DC
                   "CSF3R", "FCGR3B",  "NCF2",  #neutrophils "CXCR2","CD177", "CEACAM8", "ELANE",
                   "CD79A","MS4A1","CD19","IGKC","JCHAIN",#B
                   "MS4A2", "TPSB2","TPSAB1", #Mast
                   "MKI67","TOP2A", "CENPF", "UBE2C", #cycling
                   "THBD","SIRPA","F13A1", #not specific
                   "IL4","IL13", "IL17A","IL17F", "IFNG", #not specific
                   "CLC", "FCER1A", "HDC", 'CCR3' #Eosinophils
)

rashx_markers <- c("CD3D",
                   "CCR7", "SELL", "KLF2","CD69", #Tn/Tcm/Tmm
                   "ITGAE","CXCR6", #Trm
                   "CD4","FOXP3", "TIGIT", "CTLA4", "IL2RA", #Treg
                   "CD8A","CD8B","GZMB", "NKG7", "CCL5", #CTL
                   "KLRB1", #NK/ILC
                   "GZMA", "GZMK", "KLRC1", "PRF1","KLRD1","GNLY", "CCR8", #NK
                   "IKZF3", "TBX21", "EOMES", #ILC1
                   "GATA3","MAF", "PTGDR2", "HPGDS", #ILC2
                   "RORC", "IL23R", "IL1R1", "KIT", "AHR", #ILC3
                   "TNFRSF4", "CD96","TNFRSF18","PRDM1","PDCD1","LAG3", #ac or ex
                   "BATF", "SNHG12", "ZFAS1", #Tet
                   "TRAT1","RORA","IL7R",  #not specific
                   "HLA-DRA","HLA-DRB1","CD83","IDO1", #APC
                   "CD207","EPCAM", #LC
                   "CD68","CEBPB","FCER1G","C1QB","C1QC","CD163", "KLF4","CD14","S100A9", #Mac
                   "MS4A7","LYZ","SERPINA1", "IL23A","IL1B", "CXCL3", #Mono
                   "CD1C", "CLEC10A","CLEC9A", "XCR1", "FSCN1", "LAMP3", #DC
                   "CSF3R", "FCGR3B",  "NCF2",  #neutrophils "CXCR2","CD177", "CEACAM8", "ELANE",
                   "CD79A","MS4A1","CD19","IGKC","JCHAIN",#B
                   "MS4A2", "TPSB2","TPSAB1", #Mast
                   "MKI67","TOP2A", "CENPF", "UBE2C", #cycling
                   "CLC", "FCER1A", "HDC", 'CCR3' #Eosinophils
)
DotPlot(Merge, features = rashx_markers, cols=c("#DDA0DD", "#6A0DAD", "#502380", "#290916"), assay = "RNA", col.min = 0.1, col.max = 1, dot.min=0.1, dot.scale = 1,
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
#+scale_color_gradientn(colours = viridis::magma(20), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)')
