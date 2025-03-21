library(Seurat)
library(ComplexHeatmap)
library(RColorBrewer)

Merge <- readRDS("H:/BP data analysis/Integrated with SI data/Final objects/BP_AD_Pso.rds")

donor <- c("BP1", "BP2", "BP3", "BP4","BP5","BP6","BP7","BP8")
genes <- c("TWIST1",	"LGALS1",	"IL32",	"CAPG",	"ITM2C",	"MFHAS1",	"ANXA1",	"SOS1",	"CSGALNACT1",	"LMO4",	"IFITM2", "S100A10",	"PLA2G16",	"SYNE2",	"THADA",	"NEAT1",	"IL17RB",	"RPL36A",	"ARHGAP21",	"NBAS",	"ACTG1",	"PRKX", "TGFBR3",	"TNFSF10",	"AHNAK",	"ISG15",	"RPL17",	"CD99",	"TSHZ2",	"MMP25",	"IFITM1",		"BIRC3",
           "LPCAT2",	"CRIP1",	"CLU",	"PLP2",	"ZFP36",	"ZFP36L2",	"TUBA1B",	"SLC5A3",	"TAGLN2",	 "CXCL13",	"MTRNR2L12",	"CD7",	"MGAT4A",	"FTH1",	"LAYN",	"IL17F",	"KLRB1",	"GNLY",	"CPM",	"CTSH",	"GBP5", "SOX4",	"CLEC2B",	"GZMB",	"CD2",	"ODF2L",	"LAG3",	"LRRN3",	"ARHGEF12",	"PTPN13",	"TNFAIP3",	"TRPS1",	"METRNL",
           "BTG1",	"JUN",	"SPOCK2",	"HIST1H1E",	"RBPJ",	"MAP3K4",	"H1FX",	"UBC",	"GALNT1",	"PNRC1",	"GABPB1-AS1",	"RPS26", "MUC20-OT1",	"CHN1",	"NAP1L4",	"PTMS",	"CTLA4",	"DAPK2",	"RAP1B",	"YPEL2",	"SLA2",	"CBLB",	"ADGRG1")

## donor DEGs
set1 = list()
for (j in donor){
  tryCatch({
    gc() 
    set1[[paste0(j)]] <- FindMarkers(Merge, ident.1 = j, ident.2 ="NML", group.by = "donor", verbose =TRUE, assay = "RNA", 
                                     slot = "data", subset.ident = "Trm", test.use = "MAST",
                                     min.cells.feature = 0, min.cells.group = 0, logfc.threshold = 0, 
                                     min.pct = 0, min.diff.pct = 0, features = genes)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

### ref DEGs
marker <- FindMarkers(Merge, ident.1 = "AD", ident.2 ="Pso", group.by = "STATUS", verbose =TRUE, assay = "RNA", 
                      slot = "data", subset.ident = "Trm", test.use = "MAST", min.cells.feature = 0,
                      min.cells.group = 0, logfc.threshold = 0,  min.pct = 0, min.diff.pct = 0, features = genes)


##########Add a column to split AD and PV
library(stringr)
ad_genes = c("TWIST1",	"LGALS1",	"IL32",	"CAPG",	"ITM2C",	"MFHAS1",	"ANXA1",	"SOS1",	"CSGALNACT1",	"LMO4",	"IFITM2","S100A10",	"PLA2G16",	"SYNE2",	"THADA",	"NEAT1",	"IL17RB",	"RPL36A",	"ARHGAP21",	"NBAS",	"ACTG1",	"PRKX", "TGFBR3",	"TNFSF10",	"AHNAK", "ISG15",	"RPL17",	"CD99",	"TSHZ2",	"MMP25",	"IFITM1",		"BIRC3", "LPCAT2",	"CRIP1",	"CLU",	"PLP2",	"ZFP36",	"ZFP36L2",	"TUBA1B",	"SLC5A3",	"TAGLN2")

pv_genes = c("CXCL13",	"MTRNR2L12",	"CD7",	"MGAT4A",	"FTH1",	"LAYN",	"IL17F",	"KLRB1",	"GNLY",	"CPM",	"CTSH",	"GBP5", "SOX4",	"CLEC2B",	"GZMB",	"CD2",	"ODF2L",	"LAG3",	"LRRN3",	"ARHGEF12",	"PTPN13",	"TNFAIP3",	"TRPS1",	"METRNL", "BTG1",	"JUN",	"SPOCK2",	"HIST1H1E", "RBPJ",	"MAP3K4",	"H1FX",	"UBC",	"GALNT1",	"PNRC1",	"GABPB1-AS1",	"RPS26",
             "MUC20-OT1",	"CHN1",	"NAP1L4",	"PTMS",	"CTLA4",	"DAPK2",	"RAP1B",	"YPEL2",	"SLA2", "CBLB",	"ADGRG1")
marker$group <- ifelse(rownames(marker) %in% ad_genes, "AD", "PV")


###############Integrate two DEG lists for heatmap
set1[["ADvsPso"]] <- marker

library(plyr)
for(i in 1:length(set1)){
  colnames(set1[[i]]) <- paste0( names(set1)[i], "_", colnames(set1[[i]]) )
  set1[[i]]$ROWNAMES  <- rownames(set1[[i]])
}

data <- join_all(set1, by="ROWNAMES", type="full" )
rownames(data) <- data$ROWNAMES
data$ROWNAMES <- NULL

RashX1 <- subset(data, ADvsPso_group=="AD")
RashX1 <- RashX1[, which(str_detect(colnames(RashX1), "_avg_log2FC"))]
colnames(RashX1)
colnames(RashX1) <- c("BP1","BP2","BP3","BP4", "BP5","BP6","BP7", "BP8", "ADvsPso")
RashX2 <- subset(data, ADvsPso_group=="PV")
RashX2 <- RashX2[, which(str_detect(colnames(RashX2), "_avg_log2FC"))]
colnames(RashX2)
colnames(RashX2) <- c("BP1","BP2","BP3","BP4", "BP5","BP6","BP7", "BP8", "ADvsPso")
RashX2$ADvsPso = RashX2$ADvsPso*(-1)

########################## make figures ###########################
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
Heatmap(as.matrix(RashX1),
        name = "AD", #### annotate legend
        col = colorRamp2(c(-2,-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2), c("#2D004B","#542788","#8073AC","#B2ABD2","#D8DAEB","#F7F7F7","#FEE0B6","#FDB863","#E08214","#B35806","#7F3B08")), #### set the color scales
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9),
        #cluster_columns = F,
        clustering_distance_columns = "canberra", #canberra, euclidean
        column_title = "AD-specific genes",
        show_column_names = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12))

############################# Pso genes
Heatmap(as.matrix(RashX2),
        name = "Pso", #### annotate legend
        col = colorRamp2(c(-2,-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2), c("#7F3B08","#B35806","#E08214","#FDB863","#FEE0B6","#F7F7F7","#D8DAEB","#B2ABD2","#8073AC","#542788","#2D004B")), #### set the color scales
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9),
        #cluster_columns = F,
        clustering_distance_columns = "canberra",
        column_title = "Pso-specific genes",
        show_column_names = T,
        show_row_names = T,
        row_dend_side = "left", column_dend_side = "top",
        column_title_gp = gpar(fontsize = 12))

###############################Make the hyperdimensional plot and statistic analysis
library(concaveman)
library(ggforce)
library(Seurat)
library(monocle3)
library(dplyr)
library(Rmisc)
library(ggrepel)

exp_mat <- Merge@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
cell_metadat <- Merge@meta.data #pull out cell meta-data from Seurat object
gene_annot = data.frame(Merge@assays[["RNA"]]@counts@Dimnames[[1]])#pull out gene names from Seurat object
names(gene_annot) = "gene_short_name"
row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

data.frame(names(cell_metadat))

human_human_big_cds <- new_cell_data_set(exp_mat,
                                         cell_metadata = cell_metadat,
                                         gene_metadata = gene_annot)

#create 1 giant ass cds object from all the cds's created above
rm(exp_mat)

#Subset to patients of interest
pts = c("Skin170","Skin198","Skin230","Skin231","Skin232","Skin233","Skin236","Skin165","Skin173","Skin194","Skin199","Skin211","Skin222","Skin234","Skin235","Skin175","BP1","BP2","BP3","BP4","BP5","BP6","BP7")
human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$orig.ident %in% pts & colData(human_human_big_cds)$Ident %in% "Trm"] #subset to all normal patients and a diseased patient
colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$dis
colData(human_human_big_cds)$dis_updated = gsub("Ski170","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski198","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski230","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski231","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski232","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski233","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski236","AD",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("Ski165","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski173","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski194","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski199","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski211","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski222","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski234","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("Ski235","PV",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("Skin175","BP",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("BP1","BP",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("BP2","BP",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("BP3","BP",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("BP4","BP",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("BP5","BP",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("BP6","BP",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("BP7","BP",colData(human_human_big_cds)$dis_updated)

unique(colData(human_human_big_cds)$dis_updated)

#--------------------------------------------------------------------#
#Get the genes of interest
human_human_big_cds@assays@data$counts[1:10,1:10]

#--------------------------------------------------------------------#
#AD gene data_matrix
ad_mat = t(human_human_big_cds@assays@data$counts[ad_genes,])
ad_mat[1:10,1:10]
ad_int = data.frame(colData(human_human_big_cds)$dis_updated, colData(human_human_big_cds)$donor, rowSums(ad_mat), colData(human_human_big_cds)$STATUS)
names(ad_int) = c("dis","sample","gene_sig", "status")
ad_dfc = summarySE(ad_int, measurevar='gene_sig', groupvars=c("dis", "sample", "status"))
head(ad_dfc)
names(ad_dfc)[c(5,7)] = c("ad_gene_sig","ad_se")

#--------------------------------------------------------------------#
#PV gene data_matrix
pv_mat = t(human_human_big_cds@assays@data$counts[pv_genes,])
pv_mat[1:10,1:10]
pv_int = data.frame(colData(human_human_big_cds)$dis_updated, colData(human_human_big_cds)$donor, rowSums(pv_mat), colData(human_human_big_cds)$STATUS)
names(pv_int) = c("dis","sample","gene_sig", "status")
pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis", "sample", "status"))
head(pv_dfc)
names(pv_dfc)[c(5,7)] = c("pv_gene_sig","pv_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(ad_dfc[,c(1,2,3, 5,7)], pv_dfc[,c(5,7)]) #combine the gene signatures

ggplot(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig)) +
  geom_point(alpha=1,aes(color=dis), size=2) +
  geom_text_repel(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig,label=status),box.padding = 0.4) + #add labels
  geom_errorbarh(aes(xmax = ad_gene_sig + ad_se, xmin = ad_gene_sig - ad_se, color=dis)) +
  geom_errorbar(aes(ymax = pv_gene_sig + pv_se, ymin = pv_gene_sig - pv_se, color=dis)) +
  #stat_ellipse(aes(color=dis)) +
  theme_classic() +
  ggtitle("BP Mapping") + 
  xlab("AD-specific genes") +
  ylab("PV-specific genes") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkblue",
                                        "Pso" = "darkred",
                                        "BP" = "darkorange")) +
  scale_fill_manual(name="",values = c("AD" = "darkblue",
                                       "Pso" = "darkred",
                                       "BP" = "darkorange")) +
  geom_mark_hull(aes(fill=dis, label=dis), concavity=5) +
  theme(axis.text=element_text(size=15,color='black'),
        axis.title=element_text(size=15,color='black'),
        plot.title=element_text(size=15,color='black'))


colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor

###########################################Statistic analysis ###########################################
library(raster)
#pointDistance

all_coords_mat = as.matrix(data.frame(re_int$ad_gene_sig,re_int$pv_gene_sig))
row.names(all_coords_mat) = re_int$sample

library(vegan)
dist_mat2 =vegdist(all_coords_mat, method = "canberra", diag = FALSE, upper = TRUE, p = 2)
dist_mat2 <- as.matrix(dist_mat2)
rownames(dist_mat2) = re_int$sample
colnames(dist_mat2) = re_int$sample

ind_samps = subset(re_int, dis=="BP")$sample
ind_samps
ind_samps <-as.character(ind_samps)
ind_samps