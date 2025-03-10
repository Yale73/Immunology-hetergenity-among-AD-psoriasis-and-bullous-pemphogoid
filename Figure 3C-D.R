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
  theme_bw()
p1
#根据ONTOLOGY分类信息添加分组框#
p1+facet_grid(category~., scale = 'free_y', space = 'free_y')

