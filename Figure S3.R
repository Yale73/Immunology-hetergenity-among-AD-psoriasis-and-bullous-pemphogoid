Trm <- subset(MergeF, idents="Trm")
Idents(Trm) <- Trm$dis

#######No legend
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 12),
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

VlnPlot(Trm, features = 'IL4R', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()

VlnPlot(Trm, features = 'IL2RG', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()

VlnPlot(Trm, features = 'IL13RA1', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()

VlnPlot(Trm, features = 'IL13RA2', split.by = "dis", ncol=1) +
  theme_niwot()+
  NoLegend()
