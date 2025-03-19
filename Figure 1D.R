####################### Figure 2C
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(Seurat)

Merge <- readRDS("H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")

plot_df <- as.data.frame.matrix(Merge@meta.data)

plot_df$sample_name <- plot_df$orig.ident

test_df <- plot_df %>% group_by(orig.ident, celltype, dis) %>% 
  summarise(n = n())
test_df <- test_df %>% group_by(orig.ident, dis) %>% mutate(Freq = n / sum(n))

test_df$sample_name_condition <- paste0(test_df$orig.ident)



test_df$dis <- factor(test_df$dis, levels = c("HC", "AD", "Pso", "BP"))

#test_df$celltype <- factor(test_df$celltype, levels = c("Tcm","Tex","Tmm","Treg","Trm","CTL_1","CTL_2","CTL_3","ILC","NK", "LC","Mac_1","Mac_2","Mono","migDC","moDC","Mast","Baso","B", "Cycling_1","Cycling_2"))
test_df$celltype <- factor(test_df$celltype, levels = c("Trm","Treg","Tcm","CTL_1","CTL_2","Tmm","CTL_3","moDC","NK","Mast","Mac_2","migDC","ILC","Tex","Baso","LC","Mac_1","Cycling_1","B","Mono","Cycling_2"))
ggplot(data=test_df, aes(x=sample_name_condition, y = Freq, fill=celltype)) +
  geom_bar(stat="identity") + 
  facet_grid(cols = vars(dis), 
             labeller = label_wrap_gen(width = 16,multi_line = TRUE), 
             scales="free", space = "free") +
  xlab("Individual Sample") + ylab("Proportion of CD45+ Cells for each sample") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(21)) +
  theme_bw()+
  theme(legend.position="bottom", axis.text.x  = element_blank(), axis.ticks.x = element_blank())+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))
