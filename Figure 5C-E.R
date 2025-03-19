##############################################Figure 2C################################
library(VennDiagram)
#最后使用的网页作图软件
#https://bioinformatics.psb.ugent.be/cgi-bin/liste/Venn/calculate_venn.htpl


#读入作图文件，all.txt即上述提到的记录group1-4的元素名称的文件
dat <- read.xlsx('H:/BP data analysis/Integrated with SI data/Mast/AD_Pso_BP_specific_gene.xlsx', colNames = T, rowNames = F)
dat
#以2个分组为例
#指定统计的分组列，并设置作图颜色、字体样式等
group1 = na.omit(dat$AD)
group2 = na.omit(dat$Pso)
group3 = na.omit(dat$BP)

venn_list <- list(group1, group2, group3)
names(venn_list) <- c("AD","Pso","BP")

venn.diagram(venn_list, filename = 'H:/BP data analysis/Integrated with SI data/Trm_ABP.png', imagetype = 'png', 
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),alpha = 0.50, 
             cat.col = c("#440154ff", '#21908dff', '#fde725ff'), cat.cex = 1.5, cat.fontfamily = "sans", 
             cat.default.pos = "outer",
             cat.pos = c(-27, 0, 27),
             cat.dist = c(0.055, 0.055, 0.055),
             col = c("#440154ff", '#21908dff', '#fde725ff'), cex = 1.5, fontfamily = "sans")


##############################################Figure 2D#########################################################
#包下载与加载#
#install.packages("stringi")
#install.packages("GOplot")
library(stringi)       #处理表格数据的包
library(GOplot)
#数据载入#
#气泡图——多组数据#
go_enrich = read.xlsx("H:/BP data analysis/Integrated with SI data/80 percentage DEG/enrichment.all.xlsx",sheet= 2,sep=',') 
head(go_enrich)

#数据处理#
go_enrich$term <- paste(go_enrich$category, go_enrich$`term,description`, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term, ordered = T) #转成因子，防止重新排列
head(go_enrich)
#绘图#
#纵向柱状图-根据ONTOLOGY类型绘制#
p1 <- ggplot(go_enrich,
             aes(x=`term,description`,y=`false,discovery,rate`, fill=category)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#9C6BA3", "#98C897", "#FE9C9D",  "#9DBAD2") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab(" ") +  #x轴标签
  ylab("FDR") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(),
        axis.text=element_text(size=12, colour = "black"),
        axis.title=element_text(size=12, colour = "black"))
p1
#根据ONTOLOGY分类信息添加分组框#
p1+facet_grid(category~., scale = 'free_y', space = 'free_y')


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