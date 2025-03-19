####################################Figure 6A-C################################

######################################BP的弦图##################################
# 1.读取文件
df = read.csv('H:/BP data analysis/Integrated with SI data/CellChat/Ligand-receptor/BP.csv',header = T)

#---------------------
# 2. 整体布局与初始化
library(circlize)
circos.clear()  #这个命令用于清空画布，画错时要运行此命令重新再画
#整体布局
circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.02,0,0.02,0))
#画第一圈
fa = df$gene_id
fa = factor(fa,levels = fa)
circos.initialize(factors = fa, xlim = c(0,1))# 初始化

max(df$fc) #[1] 10.57101  
min(df$fc) #[1] -8.30034

#---------------------
# 设置fc的大小对应的颜色，随着颜色从黑到黄到红过渡，fc值从-10至0至10
col_fun = colorRamp2(c(-9, 0, 10), c("black", "yellow", "red"))
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.15, bg.border = 'black', bg.col = col_fun(df$fc),
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

# 标注基因
for(i in 1:nrow(df)){
  circos.axis(sector.index= df[i,4], direction = "outside", labels=df[i,3], 
              labels.facing = "clockwise",labels.cex=0.6, col = 'black',minor.ticks=0, major.at=seq(1, length(df$gene)))
}
# labels.away.percentage=0.1, 

#---------------------
# 画第二圈 细胞类型
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.15, bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

# cell
highlight.sector(as.character(df$gene_id[1:15]), track.index = 2,
                 text = 'Tex', niceFacing = F, font = 2, col = '#CCEBC5')
highlight.sector(as.character(df$gene_id[16:26]), track.index = 2,
                 text = 'Baso', niceFacing = F, font = 2, col ='#FFFFB3')
highlight.sector(as.character(df$gene_id[27:110]), track.index = 2,
                 text = 'Mac1', niceFacing = F, font = 2)
highlight.sector(as.character(df$gene_id[111:174]), track.index = 2,
                 text = 'Mac2', niceFacing = F, font = 2, col = '#ff9a9b')
highlight.sector(as.character(df$gene_id[175:204]), track.index = 2,
                 text = 'migDC', niceFacing = F, font = 2, col = '#FFCC99')
highlight.sector(as.character(df$gene_id[205:306]), track.index = 2,
                 text = 'moDC', niceFacing = F, font = 2, col = '#ccccff')
highlight.sector(as.character(df$gene_id[307:353]), track.index = 2,
                 text = 'Trm', niceFacing = F, font = 2, col = '#c94e65')

#---------------------
#画第三圈 配体与受体
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.08, bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

#红蓝配
highlight.sector(as.character(df$gene_id[1:11]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[12:15]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[16:20]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[21:26]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[27:68]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[69:110]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[111:147]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[148:174]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[175:193]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[194:204]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[205:257]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[258:306]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[307:332]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[333:353]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')


#画互作关系-箭头
## 读取数据库中配体-受体对关系文件
lr = read.csv('H:/BP data analysis/Integrated with SI data/CellChat/Ligand-receptor/ligand_receptors.csv',header = TRUE)
## 将配体和受体信息分开
ldf = subset(df,lr == 'ligand')
rdf = subset(df,lr == 'receptor')
##创建新的数据框
lrid = lr
lrid = within(lrid,{geneid1 = NA})
lrid = within(lrid,{geneid2 = NA})
lrid = lrid[,c('Ligand','geneid1','Receptor','geneid2')]
## 添加geneid，用于后面画箭头的索引使用。
for(i in 1:nrow(rdf)){
  index = which(lrid$Receptor %in% rdf[i,3])
  for(ind in index){lrid[ind,4] = as.character(rdf[i,4])}
}

for(i in 1:nrow(ldf)){
  index = which(lrid$Ligand %in% ldf[i,3])
  for(ind in index){lrid[ind,2] = as.character(ldf[i,4])}
}
##去除NA
lrid = na.omit(lrid)

##画图
for(i in 1:nrow(lrid)){
  circos.link(sector.index1 = lrid[i,2],point1 = 0, sector.index2 = lrid[i,4], point2 = 0 ,directional = 1,
              h=1, lwd=1, col="black",lty=1)
}



######################################AD的弦图##################################
# 1.读取文件
df = read.csv('H:/BP data analysis/Integrated with SI data/CellChat/Ligand-receptor/AD.csv',header = T)

#---------------------
# 2. 整体布局与初始化
library(circlize)
circos.clear()  #这个命令用于清空画布，画错时要运行此命令重新再画
#整体布局
circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.02,0,0.02,0))
#画第一圈
fa = df$gene_id
fa = factor(fa,levels = fa)
circos.initialize(factors = fa, xlim = c(0,1))# 初始化

max(df$fc) #[1] 10.57101  
min(df$fc) #[1] -8.30034

#---------------------
# 设置fc的大小对应的颜色，随着颜色从黑到黄到红过渡，fc值从-10至0至10
col_fun = colorRamp2(c(-10, 0, 10), c("black", "yellow", "red"))
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.15, bg.border = 'black', bg.col = col_fun(df$fc),
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

# 标注基因
for(i in 1:nrow(df)){
  circos.axis(sector.index= df[i,4], direction = "outside", labels=df[i,3], 
              labels.facing = "clockwise",labels.cex=0.7, col = 'black',minor.ticks=0, major.at=seq(1, length(df$gene)))
}
# labels.away.percentage=0.1, 

#---------------------
# 画第二圈 细胞类型
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.15, bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

# cell
highlight.sector(as.character(df$gene_id[1:2]), track.index = 2,
                 text = 'Mac1', niceFacing = F, font = 2)
highlight.sector(as.character(df$gene_id[3:15]), track.index = 2,
                 text = 'Mac2', niceFacing = F, font = 2, col = '#ff9a9b')
highlight.sector(as.character(df$gene_id[16:36]), track.index = 2,
                 text = 'migDC', niceFacing = F, font = 2, col = '#FFCC99')
highlight.sector(as.character(df$gene_id[37:54]), track.index = 2,
                 text = 'moDC', niceFacing = F, font = 2, col = '#ccccff')
highlight.sector(as.character(df$gene_id[55:82]), track.index = 2,
                 text = 'Trm', niceFacing = F, font = 2, col = '#c94e65')

#---------------------
#画第三圈 配体与受体
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.08, bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

#红蓝配
highlight.sector(as.character(df$gene_id[1:2]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[3:12]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[13:15]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[16:33]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[34:36]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[37:51]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[52:54]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[55:73]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[74:82]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')

#画互作关系-箭头
## 读取数据库中配体-受体对关系文件
lr = read.csv('H:/BP data analysis/Integrated with SI data/CellChat/Ligand-receptor/ligand_receptors.csv',header = TRUE)
## 将配体和受体信息分开
ldf = subset(df,lr == 'ligand')
rdf = subset(df,lr == 'receptor')
##创建新的数据框
lrid = lr
lrid = within(lrid,{geneid1 = NA})
lrid = within(lrid,{geneid2 = NA})
lrid = lrid[,c('Ligand','geneid1','Receptor','geneid2')]
## 添加geneid，用于后面画箭头的索引使用。
for(i in 1:nrow(rdf)){
  index = which(lrid$Receptor %in% rdf[i,3])
  for(ind in index){lrid[ind,4] = as.character(rdf[i,4])}
}

for(i in 1:nrow(ldf)){
  index = which(lrid$Ligand %in% ldf[i,3])
  for(ind in index){lrid[ind,2] = as.character(ldf[i,4])}
}
##去除NA
lrid = na.omit(lrid)

##画图
for(i in 1:nrow(lrid)){
  circos.link(sector.index1 = lrid[i,2],point1 = 0, sector.index2 = lrid[i,4], point2 = 0 ,directional = 1,
              h=1, lwd=1, col="black",lty=1)
}

######################################Pso的弦图##################################
# 1.读取文件
df = read.csv('H:/BP data analysis/Integrated with SI data/CellChat/Ligand-receptor/Pso.csv',header = T)

#---------------------
# 2. 整体布局与初始化
library(circlize)
circos.clear()  #这个命令用于清空画布，画错时要运行此命令重新再画
#整体布局
circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.02,0,0.02,0))
#画第一圈
fa = df$gene_id
fa = factor(fa,levels = fa)
circos.initialize(factors = fa, xlim = c(0,1))# 初始化

max(df$fc) #[1] 10.57101  
min(df$fc) #[1] -8.30034

#---------------------
# 设置fc的大小对应的颜色，随着颜色从黑到黄到红过渡，fc值从-10至0至10
col_fun = colorRamp2(c(-10, 0, 10), c("black", "yellow", "red"))
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.15, bg.border = 'black', bg.col = col_fun(df$fc),
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

# 标注基因
for(i in 1:nrow(df)){
  circos.axis(sector.index= df[i,4], direction = "outside", labels=df[i,3], 
              labels.facing = "clockwise",labels.cex=0.7, col = 'black',minor.ticks=0, major.at=seq(1, length(df$gene)))
}
# labels.away.percentage=0.1, 

#---------------------
# 画第二圈 细胞类型
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.15, bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

# cell
highlight.sector(as.character(df$gene_id[1:3]), track.index = 2,
                 text = 'Mac1', niceFacing = F, font = 2)
highlight.sector(as.character(df$gene_id[4:9]), track.index = 2,
                 text = 'Mac2', niceFacing = F, font = 2, col = '#ff9a9b')
highlight.sector(as.character(df$gene_id[10:21]), track.index = 2,
                 text = 'migDC', niceFacing = F, font = 2, col = '#FFCC99')
highlight.sector(as.character(df$gene_id[22:29]), track.index = 2,
                 text = 'moDC', niceFacing = F, font = 2, col = '#ccccff')
highlight.sector(as.character(df$gene_id[30:57]), track.index = 2,
                 text = 'Trm', niceFacing = F, font = 2, col = '#c94e65')

#---------------------
#画第三圈 配体与受体
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.08, bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

#红蓝配
highlight.sector(as.character(df$gene_id[1:1]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[2:3]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[4:8]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[9:9]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[10:21]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[22:26]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[27:29]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')
highlight.sector(as.character(df$gene_id[30:48]), track.index = 3,
                 text = 'ligand', niceFacing = F,col = '#FB8072',text.col = 'white')
highlight.sector(as.character(df$gene_id[49:57]), track.index = 3,
                 text = 'receptor', niceFacing = F,col = '#80B1D3',text.col = 'white')

#画互作关系-箭头
## 读取数据库中配体-受体对关系文件
lr = read.csv('H:/BP data analysis/Integrated with SI data/CellChat/Ligand-receptor/ligand_receptors.csv',header = TRUE)
## 将配体和受体信息分开
ldf = subset(df,lr == 'ligand')
rdf = subset(df,lr == 'receptor')
##创建新的数据框
lrid = lr
lrid = within(lrid,{geneid1 = NA})
lrid = within(lrid,{geneid2 = NA})
lrid = lrid[,c('Ligand','geneid1','Receptor','geneid2')]
## 添加geneid，用于后面画箭头的索引使用。
for(i in 1:nrow(rdf)){
  index = which(lrid$Receptor %in% rdf[i,3])
  for(ind in index){lrid[ind,4] = as.character(rdf[i,4])}
}

for(i in 1:nrow(ldf)){
  index = which(lrid$Ligand %in% ldf[i,3])
  for(ind in index){lrid[ind,2] = as.character(ldf[i,4])}
}
##去除NA
lrid = na.omit(lrid)

##画图
for(i in 1:nrow(lrid)){
  circos.link(sector.index1 = lrid[i,2],point1 = 0, sector.index2 = lrid[i,4], point2 = 0 ,directional = 1,
              h=1, lwd=1, col="black",lty=1)
}




###################################Figure S6####################################
#######HC
p1 <- plotGeneExpression(cellchat.HC, signaling = 'IFN-II', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IFN-II signaling in HC")
p2 <- plotGeneExpression(cellchat.HC, signaling = 'TNF', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("TNF signaling in HC")
p3 <- plotGeneExpression(cellchat.HC, signaling = 'IL4', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL4 signaling in HC")
p4 <- plotGeneExpression(cellchat.HC, signaling = 'IL17', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL17 signaling in HC")


p5 <- plotGeneExpression(cellchat.BP, signaling = 'IFN-II', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IFN-II signaling in BP")
p6 <- plotGeneExpression(cellchat.BP, signaling = 'TNF', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("TNF signaling in BP")
p7 <- plotGeneExpression(cellchat.BP, signaling = 'IL4', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL4 signaling in BP")
p8 <- plotGeneExpression(cellchat.BP, signaling = 'IL17', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL17 signaling in BP")


p9 <- plotGeneExpression(cellchat.AD, signaling = 'IFN-II', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IFN-II signaling in AD")
p10 <- plotGeneExpression(cellchat.AD, signaling = 'TNF', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("TNF signaling in AD")
p11 <- plotGeneExpression(cellchat.AD, signaling = 'IL4', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL4 signaling in AD")
p12 <- plotGeneExpression(cellchat.AD, signaling = 'IL17', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL17 signaling in AD")

p13 <- plotGeneExpression(cellchat.Pso, signaling = 'IFN-II', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IFN-II signaling in Pso")
p14 <- plotGeneExpression(cellchat.Pso, signaling = 'TNF', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("TNF signaling in Pso")
p15 <- plotGeneExpression(cellchat.Pso, signaling = 'IL4', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL4 signaling in Pso")
p16 <- plotGeneExpression(cellchat.Pso, signaling = 'IL17', type = 'dot', color.use = c("white", "#b2182b"))+ggtitle("IL17 signaling in Pso")

plot_grid(p1, p2, p3, p4, ncol =2)
plot_grid(p5, p6, p7, p8, ncol = 2)
plot_grid(p9, p10, p11, p12, ncol = 2)
plot_grid(p13, p14, p15, p16, ncol = 2)
