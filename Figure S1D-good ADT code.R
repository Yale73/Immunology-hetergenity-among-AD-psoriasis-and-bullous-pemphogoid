########ADT CALCULATION FUNCTION
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['ADT']]@scale.data
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>1.0)/ncells
  }else{return(NA)}
}

genes <- c("CD19","CD25","CD45RO","CD80","CD123","HLA-DR",
           "CD1c","CD161","CD141","CD294","CD127","CD207",
           "CD197","CD196","CD21","CD195","CD56","CD14", 
           "CD11c", "CD3", "CD4", "CD45RA","CD69","CD8")

###############few samples
Idents(Merge) <- Merge$donor

GoodADT <- subset(Merge, idents = c("154", "155", "194", "203", "222", "234"), invert=T)

Idents(GoodADT) <- GoodADT$Ident1

A <- PrctCellExpringGene(GoodADT, genes=genes, group.by="Ident")
library(openxlsx)
write.xlsx(A, "H:/BP data analysis/Integrated with SI data/Our with 1.0 cutoff-good-donor-Ident.xlsx", colNames=T, rowNames=T)

library(openxlsx)
ADT <- read.xlsx("H:/BP data analysis/Integrated with SI data/Our with 1.0 cutoff-good-donor-Ident.xlsx", sheet = 1, colNames = T, rowNames = T)
#colnames(ADT)#check colnames
#rownames(ADT) <- ADT[,1]# change first column to rownames
#ADT[,1] <- NULL# remove the first column
#colnames(ADT)#double check colnames

breaksList = as.numeric(seq(0.1, 0.6, by = .01))
library(viridis)
pheatmap(ADT, color = inferno(length(breaksList)),
         show_rownames = T, 
         show_colnames = T,
         cluster_cols =F,
         cluster_rows =T,
         fontsize_row = 10, 
         fontsize_col = 10, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "CITE-seq")
