library(Seurat)
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
library(Seurat)
library(cowplot)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
##########################
MergeF <- readRDS("H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")
DimPlot(MergeF, group.by = "celltype")
Idents(MergeF) <- MergeF$celltype
##########################################T cell integration ########################################
Mast <- subset(Merge, idents=c("Mast", "Baso", "Cycling_2"))
DimPlot(Mast, group.by = "celltype")
Mast$celltype <- droplevels(Mast$celltype)

Mast$dis <- droplevels(Mast$dis)
Mast$celltype <- droplevels(Mast$celltype)
saveRDS(Mast, "H:/BP data analysis/Integrated with SI data/Mast/Mast.rds")

#################################Figure 6A##########################################
DimPlot(Mast, label = T, label.size = 4.5, group.by = "celltype")+NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(3))+
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


################Figure S4###########
DimPlot(Mast, label = T, label.size = 4.5, group.by = "celltype", split.by = "dis", ncol = 2, repel = T)+NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(3))+
  theme(text=element_text(family="Arial",size=12)) +
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
        legend.key=element_blank())+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))

#############################Figure 6B#########################################
rashx_markers <- list(Mast=c("TPSB2","TPSAB1",'CPA3', 'KIT'), #Mast
                      Baso=c("FCER1A", "CCR3", "IL4", 'IL13'), #eosinophil
                      Cycling=c("MKI67","UBE2C") #cycling
)

DotPlot(Mast, features = rashx_markers, cols=c("#DDA0DD", "#6A0DAD", "#502380", "#290916"), assay = "RNA", col.min = 0.1, col.max = 1, dot.min=0.02, dot.scale = 1,
        cluster.idents=T, group.by = "celltype")+
  scale_size(range = c(0, 5))+ 
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, face = "italic"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 11 ),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9),
        # 分面标题
        strip.background = element_rect(color="#DDA0DD", fill="#DDA0DD"),
        strip.text = element_text(margin=margin(b=3, unit="mm")),
        strip.placement = 'outlet') +
  labs(x="", y="")
#+scale_color_gradientn(colours = viridis::magma(20), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)')

################################Figure 6C##########################################
my_comparisons = list(c("HC", "AD"), c("HC", "Pso"), c("HC", "BP"), c("AD", "BP"), c("AD", "Pso"), c("Pso", "BP")) 


library(dittoSeq)
#分组箱线图1
Singlecellratio_plotstat(Mast, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'group',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =3)

###########################Figure 6D######################################
a1 <- data.frame(
  counts = c(-354, -613, -403, -1, -159, 0, 0, 0, 0), #将a1显示在x轴下方,下调基因数目
  type1 =   rep(c('Mast', 'Baso', 'Cycling'), each=3),
  type2 =rep(c('AD', 'BP', 'Pso'), each=1)
)
a1$counts2 <- abs(a1$counts) #为了让label 显示正值

a1$type1 <- factor(a1$type1, levels = c("Mast", "Baso", "Cycling"))


a2 <- data.frame(
  counts = c(223, 1745, 193, 1, 257,0, 0, 0, 0),#将a2显示在x轴上方,上调基因数目
  type1 =   rep(c('Mast', 'Baso', 'Cycling'), each=3),
  type2 =rep(c('AD', 'BP', 'Pso'), each=1)
)

a2$type1 <- factor(a2$type1, levels = c("Mast", "Baso", "Cycling"))


ggplot() + 
  geom_col(data = a2, aes(type1, counts, fill = type2),
           position = "dodge") +
  geom_col(data = a1, 
           aes(type1, counts, fill = type2),
           position = "dodge") +
  geom_text(data = a2, 
            aes(type1, counts,fill = type2, label = counts),
            position = position_dodge2(0.9), vjust = -0.1) + 
  geom_text(data = a1, 
            aes(type1, counts,fill = type2, label = counts2),
            position = position_dodge2(0.9), vjust = 1) + 
  ggthemes::theme_economist() + 
  # 在这个网站找颜色https://colorbrewer2.org/
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Blues")[7], 
                               RColorBrewer::brewer.pal(9, "PuRd")[8],
                               RColorBrewer::brewer.pal(9, "Greens")[8])) + 
  labs(x = NULL, y="DEG numbers") + 
  theme( axis.ticks.x = element_blank())+ 
  theme_bw()

#######################################Volcanoplot#########################
library(openxlsx)
library(ggVolcano)
nrDEG <- read.xlsx('H:/BP data analysis/Integrated with SI data/Mast/BPvsHC_DEG.xlsx', colNames=T, rowNames=T, sheet = 2)
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

data_keep_rows <- c("CLEC12A", "CLC", "CD9", "CCR3", "IL17RA", "CXCL8", "CCL4", "CXCR4", "CCL22") 
for_label <- data[rownames(data) %in% data_keep_rows, ] 

ggplot(data = data,aes(x = avg_log2FC,y = (-log10(p_val_adj)))) +
  geom_point(alpha=0.4, size=3.5,aes(color=regulate)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw()+  
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(aes(label = row), data = for_label,color="black", max.overlaps = 20)+
  ggtitle("DEGs in Basophils")


##########################Figure 6F## KEGG plot
##########eosinophils
library(stringi)       #处理表格数据的包
library(GOplot)
#数据载入#
#气泡图——多组数据#
go_enrich = read.xlsx("H:/BP data analysis/Integrated with SI data/Mast/Eosinophil_logFC_more than 1.xlsx",sheet= 3,sep=',') 
head(go_enrich)

#数据处理#
go_enrich$term <- paste(go_enrich$Category, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term,ordered = T) #转成因子，防止重新排列
head(go_enrich)
#绘图#
#纵向柱状图-根据ONTOLOGY类型绘制#
p1 <- ggplot(go_enrich,
             aes(x=Description,y=`Log(q-value)`, fill=Category)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#ff9a9b", "#9998ff", "#f3b169", "#37ab78") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Log(q-value)") +  #y轴标签
  labs(title = "GO Terms Enrich in Basophils")+  #设置标题
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(),
        axis.text=element_text(size=12, colour = "black"),
        axis.title=element_text(size=12, colour = "black"))
p1
#根据ONTOLOGY分类信息添加分组框#
p1+facet_grid(Category~., scale = 'free_y', space = 'free_y')


############绘制横向
#横向柱状图-根据ONTOLOGY类型绘制#
p2 <- ggplot(go_enrich, 
             aes(x=Description,y=`Log(q-value)`, fill=Category)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#ff9a9b", "#9998ff", "#f3b169", "#37ab78") ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(),
        axis.text=element_text(size=12, colour = "black", angle = 90),
        axis.title=element_text(size=12, colour = "black"),
        axis.text.x=element_text(color="black",angle = 90,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）

p2 
#根据ONTOLOGY分类信息添加分组框#
p2 + facet_grid(.~Category, scale = 'free_x', space = 'free_x')

#########################No legend
Basophil <- subset(Mast, idents="Baso")
gene <- list(Cytokine_Signaling_in_Immune_system=c("ADAR", "BIRC3",	"CSF2RB",	"DUSP6",	"EGR1",	"EIF4A1",	"FLT3", "IFI6",	"HLA-A",	"HLA-C",	
             "IFIT2",	"IFIT3",	"IL3RA",		"IRF1",	"LCP1",		"PSMB9", "CCL4",	"SP100", "ELOB",
             "TRIM25",	"TNFSF11", "ISG15",	"BCL2L11",	"IRF9",	"IFITM2",		"PDCD4",	"XAF1",	"RSAD2"),
             IL_17_signal=c("CASP3",	"FOSB",	"HSP90AA1",	"CXCL8",	"JUND",	"TNFAIP3",	"IL17RA","RUNX1",	"NFIL3"),
             IL4_signal=c("CEBPB", "MAPK14", "FOS", "NFKBIA", "PTPN6", "STAT1",	"IRS2",	"DOK2"))


gene  <- list(Type_1=c("IL1B",	"CCL3",	"CCL4",	"CCL5",	"CCR5",	"CXCL8",	"CXCR3",	"IL18",	"IFNG",	"IL12A",	"CXCL10",	"CXCL11",	"CXCL16",	"CCR2",
                       "IFNGR1",	"IFNGR2",	"CXCL9",	"CCR7",	"XCR1",	"XCL1",	"IL2",	"IL6",	"IL15",	"CSF2",	"TNF",	"STAT1"),
              Type_2=c("CCL18",	"CCL17",	"CCL13",	"CCL14",	"CTSC",	"CCL26",	"IL4",	"IL5",	"IL19",	"CCL11",	"CCR3",
                       "CCR4",	"CCR8",	"IL13",	"TSLP",	"IL25",	"IL33",	"GATA3",	"IL9",	"IL31"),
              Type_3=c("IL17A",	"IL17F",	"IL22",	"CCL20",	"CCR6",	"CXCL2",	"CXCR1",	"CXCR2",	"RORC")
)


Basophil$dis <- factor(Basophil$dis, levels = c("HC", "AD", "Pso", "BP"))

DotPlot(Basophil, features = gene, split.by = "dis", group.by = "celltype", cols = c(RColorBrewer::brewer.pal(9, "Oranges")[7], 
                                                                                       RColorBrewer::brewer.pal(9, "Blues")[7], 
                                                                                       RColorBrewer::brewer.pal(9, "PuRd")[8],
                                                                                       RColorBrewer::brewer.pal(9, "Greens")[8]), 
        assay = "RNA", col.min = 0.1, col.max = 1, dot.min = 0.01, cluster.idents=F)+
  scale_size(range = c(0, 10))+
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, family="TT Times New Roman", face = "italic"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 9, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9),
        legend.position = "bottom") +
  labs(x="", y="")

#################
Baso <- subset(Mast, idents="Baso")
Idents(Baso) <- Baso$orig.ident
IL4 <- subset(Baso, IL4>0)
table(Idents(IL4))

IL13<- subset(Baso, IL13>0)
table(Idents(IL13))