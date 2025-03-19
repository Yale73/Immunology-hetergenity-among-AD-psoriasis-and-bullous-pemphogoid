Merge <- readRDS("H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")

##################################Figure 2A#####################################
Idents(Merge) <- Merge$celltype
       
levels(Merge$celltype)
#"Trm"       "Treg"      "Tcm"       "CTL_1"     "CTL_2"     "Tmm"       "CTL_3"     "moDC"      "NK" 
#"Mast"      "Mac_2"     "migDC"     "ILC"       "Tex"       "Baso"      "LC"        "Mac_1"     "Cycling_1"
#"B"         "Mono"      "Cycling_2"
Tsubset <- subset(Merge, idents= c("Trm", "Treg", "Tcm",  "CTL_1", "CTL_2", "Tmm", "CTL_3", "NK", "ILC", "Tex", "Cycling_1"))
Idents(Tsubset) <- Tsubset$celltype
Tsubset$celltype <- droplevels(Tsubset$celltype)

DimPlot(Tsubset, label = T, label.size = 4, group.by = "celltype")+NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
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
        legend.text=element_text(family="Arial", size=13),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


DimPlot(Tsubset, label = T, label.size = 4, group.by = "celltype", split.by = "dis", ncol = 2)+NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
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
        legend.text=element_text(family="Arial", size=13),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))


#######################Figure 2B
FeaturePlot(Tsubset, features = "FOXP3", order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = c("grey", "red", "darkred"))+
  theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title= element_text(color="black",hjust = 0.5, face="italic", family = "Arial", size=12),
        panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=11),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=11),
        axis.title.y=element_text(colour='black', size=11),
        axis.text=element_text(colour='black',size=11),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=11),
        legend.key=element_blank())

FeaturePlot(Tsubset, features = "CD3D", order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = c("grey", "red", "darkred"))+
  theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title= element_text(color="black",hjust = 0.5, face="italic", family = "Arial", size=12),
        panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=11),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=11),
        axis.title.y=element_text(colour='black', size=11),
        axis.text=element_text(colour='black',size=11),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=11),
        legend.key=element_blank())

FeaturePlot(Tsubset, features = "CD8A", order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = c("grey", "red", "darkred"))+
  theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title= element_text(color="black",hjust = 0.5, face="italic", family = "Arial", size=12),
        panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=11),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=11),
        axis.title.y=element_text(colour='black', size=11),
        axis.text=element_text(colour='black',size=11),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=11),
        legend.key=element_blank())


FeaturePlot(Tsubset, features = "CD4", order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = c("grey", "red", "darkred"))+
  theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title= element_text(color="black",hjust = 0.5, face="italic", family = "Arial", size=12),
        panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=11),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=11),
        axis.title.y=element_text(colour='black', size=11),
        axis.text=element_text(colour='black',size=11),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=11),
        legend.key=element_blank())

FeaturePlot(Tsubset, features = "CD69", order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = c("grey", "red", "darkred"))+
  theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title= element_text(color="black",hjust = 0.5, face="italic", family = "Arial", size=12),
        panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=11),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=11),
        axis.title.y=element_text(colour='black', size=11),
        axis.text=element_text(colour='black',size=11),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=11),
        legend.key=element_blank())

FeaturePlot(Tsubset, features = "ITGAE", order = T, pt.size = 0.5)+
  scale_colour_gradientn(colours = c("grey", "red", "darkred"))+
  theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title= element_text(color="black",hjust = 0.5, face="italic", family = "Arial", size=12),
        panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=11),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=11),
        axis.title.y=element_text(colour='black', size=11),
        axis.text=element_text(colour='black',size=11),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=11),
        legend.key=element_blank())

##################################Figure 2C Percentage change###################
#分组箱线图1
my_comparisons = list(c("HC", "AD"), c("HC", "Pso"), c("HC", "BP"), c("AD", "BP"), c("AD", "Pso"), c("Pso", "BP")) 

Singlecellratio_plotstat(Tsubset, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'group',
                         group_by.point = "orig.ident",label.x = 1, pt.size =0.2,
                         label = 'p.format', ncol = 2)

##################################Figure 2D DEG number###################
a1 <- data.frame(
  counts = c(-1779, -1951, -1671, -1287, -1000, -1042, -2291, -550, -2083, -1237, -816, -1715, -1051, -1055, -1378, -788, -580, -1172, 
             -524, -534, -726, -741, -397, -963, -288, -2390, -416, -39, -688, -111, 0, 0, 0), #将a1显示在x轴下方,下调基因数目
  type1 =   rep(c("Trm", "Treg", "Tcm", "CTL_3", "CTL_2", "Tmm", "CTL_1", "NK", "ILC", "Tex", "Cycling_1"), each=3),
  type2 =rep(c('AD', 'BP', 'Pso'), 11)
)
a1$counts2 <- abs(a1$counts) #为了让label 显示正值

a1$type1 <- factor(a1$type1, levels = c("Trm", "Treg", "Tcm", "CTL_3", "CTL_2", "Tmm", "CTL_1", "NK", "ILC", "Tex", "Cycling_1"))


a2 <- data.frame(
  counts = c(2373, 2188, 2928, 1367, 1807, 1820, 998, 1855, 12032, 2135, 1614, 2302, 2221, 1099, 2944, 1640, 1993, 13254, 
             663, 1071, 2140, 2296, 722, 1608, 550, 347, 518, 18, 802, 55, 0, 0, 0),#将a2显示在x轴上方,上调基因数目
  type1 =  rep(c("Trm", "Treg", "Tcm", "CTL_3", "CTL_2", "Tmm", "CTL_1", "NK", "ILC", "Tex", "Cycling_1"), each=3),
  type2 =rep(c('AD', 'BP', 'Pso'), 11)
)

a2$type1 <- factor(a2$type1, levels = c("Trm", "Treg", "Tcm", "CTL_3", "CTL_2", "Tmm", "CTL_1", "NK", "ILC", "Tex", "Cycling_1"))


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
  theme_bw() +   #与下面的theme顺序很重要，否则字体无法设置
  theme(panel.background = element_rect(color ='black',fill = 'transparent'),
        panel.border = element_blank(),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(size = 1, color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black", size = 0.5), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=13),
        legend.key=element_blank())+
  geom_hline(yintercept = 0, size=0.5)


################################Figure 2E type1,2,3 signature###################
library(Seurat)
library(UCell)
library(ggrastr)
library(stringr)
library(ggplot2)
library(viridis)
library(ggpubr)

DimPlot(Tsubset,label = T)+NoLegend()

#UCell评分函数AddModuleScore_UCell可以提供方seurat对象
#评分基因set是以list的形式提供
markers <- list()
markers$Type_1_signature <- c("IL1B",	"CCL3",	"CCL4",	"CCL5",	"CCR5",	"CXCL8",	"CXCR3",	"IL18",	"IFNG",	"IL12A",	"CXCL10",	"CXCL11",	"CXCL16",	"CCR2",
                              "IFNGR1",	"IFNGR2",	"CXCL9",	"CCR7",	"XCR1",	"XCL1",	"IL2",	"IL6",	"IL15",	"CSF2",	"TNF",	"STAT1")
markers$Type_2_signature <- c("CCL18",	"CCL17",	"CCL13",	"CCL14",	"CTSC",	"CCL26",	"IL4",	"IL5",	"IL19",	"CCL11",	"CCR3",
                              "CCR4",	"CCR8",	"IL13",	"TSLP",	"IL25",	"IL33",	"GATA3",	"IL9",	"IL31")
markers$Type_3_signature <- c("IL17A",	"IL17F",	"IL22",	"CCL20",	"CCR6",	"CXCL2",	"CXCR1",	"CXCR2",	"RORC")

###########################prepare the objects#################################
Idents(Tsubset) <- Tsubset$dis
Tsubset$dis <- droplevels(Tsubset$dis)
HC <- subset(Tsubset, idents="HC")
HC$dis <- droplevels(HC$dis)
Idents(HC) <- HC$celltype

AD <- subset(Tsubset, idents="AD")
AD$dis <- droplevels(AD$dis)
Idents(AD) <- AD$celltype

BP <- subset(Tsubset, idents="BP")
BP$dis <- droplevels(BP$dis)
Idents(BP) <- BP$celltype

Pso <- subset(Tsubset, idents="Pso")
Pso$dis <- droplevels(Pso$dis)
Idents(Pso) <- Pso$celltype

############################# type 1 for HC
library(UCell)
library(ggrastr)
genes <- list(c( "IFNG",	"IL12A", "TNF"))
names(genes) <- 'Type_1_signature'
gene_score <- AddModuleScore_UCell(HC,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_1_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype, y=Type_1_signature_score, fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 1 signature score in HC")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较
############################# type 1 for AD
genes <- list(c("IFNG",	"IL12A", "TNF"))
names(genes) <- 'Type_1_signature'
gene_score <- AddModuleScore_UCell(AD,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_1_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_1_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 1 signature score in AD")+ 
  geom_jitter_rast(col="#00000033",pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较
############################# type 1 for BP
genes <- list(c("IFNG",	"IL12A", "TNF"))
names(genes) <- 'Type_1_signature'
gene_score <- AddModuleScore_UCell(BP,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_1_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_1_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 1 signature score in BP")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

############################# type 1 for Pso
genes <- list(c("IFNG",	"IL12A", "TNF"))
names(genes) <- 'Type_1_signature'
gene_score <- AddModuleScore_UCell(Pso,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_1_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_1_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 1 signature score in Pso")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

############################# type 2
########HC
genes <- list(c("IL4",	"IL5", "IL13", "IL31"))
names(genes) <- 'Type_2_signature'
gene_score <- AddModuleScore_UCell(HC,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_2_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_2_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 2 signature score in HC")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较
########BP
genes <- list(c("IL4",	"IL5", "IL13", "IL31"))
names(genes) <- 'Type_2_signature'
gene_score <- AddModuleScore_UCell(BP,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_2_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_2_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 2 signature score in BP")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

########AD
genes <- list(c("IL4",	"IL5", "IL13", "IL31"))
names(genes) <- 'Type_2_signature'
gene_score <- AddModuleScore_UCell(AD,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_2_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_2_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 2 signature score in AD")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

########Pso
genes <- list(c("IL4",	"IL5", "IL13", "IL31"))
names(genes) <- 'Type_2_signature'
gene_score <- AddModuleScore_UCell(Pso, features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_2_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_2_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 2 signature score in Pso")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

####################type 3
##########HC
genes <- list(c("IL17A",	"IL17F",	"IL22"))
names(genes) <- 'Type_3_signature'
gene_score <- AddModuleScore_UCell(HC,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_3_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_3_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 3 signature score in HC")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

##########BP
genes <- list(c("IL17A",	"IL17F",	"IL22"))
names(genes) <- 'Type_3_signature'
gene_score <- AddModuleScore_UCell(BP,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_3_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_3_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 3 signature score in BP")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1,position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

##########AD
genes <- list(c("IL17A",	"IL17F",	"IL22"))
names(genes) <- 'Type_3_signature'
gene_score <- AddModuleScore_UCell(AD,features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_3_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_3_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 3 signature score in AD")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1,position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

##########Pso
genes <- list(c("IL17A",	"IL17F",	"IL22"))
names(genes) <- 'Type_3_signature'
gene_score <- AddModuleScore_UCell(Pso, features=genes,name="_score")

df<- FetchData(gene_score,vars = c("celltype","Type_3_signature_score"))
df$celltype <- factor(df$celltype,levels = c("Trm", "Treg", "Tcm","Tmm", "Tex", "CTL_1",  "CTL_2", "CTL_3", "NK", "ILC", "Cycling_1"))#设置顺序


ggplot(df,aes(x=celltype,y=Type_3_signature_score,fill=celltype, colour = celltype))+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Type 3 signature score in Pso")+ 
  geom_jitter_rast(col="#00000033", pch=20, cex=1,position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(11))+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  stat_compare_means(method="t.test",hide.ns = T,
                     comparisons =my_comparisons,
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较


###########################
Idents(Tsubset) <- Tsubset$celltype
Tex  <- subset(Tsubset, idents= "Tex")
Tex$dis <- factor(Tex$dis, levels = c("HC", "AD", "Pso", "BP"))


VlnPlot(Tex, features = "IFNG",pt.size=0, group.by = "dis")+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))



VlnPlot(Tex, features = "TNF",pt.size=0, group.by = "dis")+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


VlnPlot(Tex, features = "IL4",pt.size=0, group.by = "dis")+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

VlnPlot(Tex, features = "IL13",pt.size=0, group.by = "dis")+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


VlnPlot(Tex, features = "IL17A",pt.size=0, group.by = "dis")+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


VlnPlot(Tex, features = "IL17F",pt.size=0, group.by = "dis")+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black", face = "italic",
                                  size = 13, margin = margin(t = 1, b = 12)),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

