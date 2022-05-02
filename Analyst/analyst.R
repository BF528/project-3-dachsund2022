# library
library(dplyr)
library(Seurat)
library(patchwork)
library(ggrastr)
# read file
df <- readRDS("GSM2230760_seurat.rda")
# Identify marker genes for each cluster
cluster1.markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25 )
signif <- cluster1.markers[cluster1.markers$p_val_adj<0.05,]
write.csv(signif, file="marker_genes.csv")

# Label clusters as a cell type based on marker genes
df_labels<-signif%>%mutate(celltype =
                             case_when(signif$cluster == "0" ~ "Delta",
                                       signif$cluster == "1" ~ "Beta_1",
                                       signif$cluster == "2" ~ "Gamma",
                                       signif$cluster == "3" ~ "Acinar",
                                       signif$cluster == "4" ~ "Alpha",
                                       signif$cluster == "5" ~ "Ductal",
                                       signif$cluster == "6" ~ "Beta_2",
                                       signif$cluster == "7" ~ "Unknown_1",
                                       signif$cluster== "8" ~ "Unknow_2",
                                       signif$cluster == "9" ~ "Endothelial",
                                       signif$cluster == "10" ~ "Macrophage_1",
                                       signif$cluster == "11" ~ "Exocrine Glandular",
                                       signif$cluster == "12" ~ "Macrophages_2"))

# 3.	Visualize the clustered cells using a projection method
png("umap_labels_num.png")
df <- RenameIdents(object = df, '0' = "0",'1' = "1",'2' = "2",'3' = "3",'4' = "4",'5' = "5",'6' = "6",'7' = "7",'8' = "8",'9' = '9','10' = '10','11' = '11','12' = '12')
DimPlot(df, reduction = "umap", label = 'true')
dev.off()          

# 4. Visualize the top marker genes per cluster
png("heatmap_clusters.png", height = 800, width = 1500)
top10 <- signif %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pl <- DoHeatmap(df, features = top10$gene) 
pl
dev.off()

# 5. novel markers: top 2 of each cluster
png("vln_jit_2.png", width=1500, height=1500)
# a <- VlnPlot(df, features = c("SST","SP100", "INS", "GNAS", "FN1", "COL1A1", "REG1B", "REG1A", "TTR", "GCG", "KRT19", "MMP7", "EEF1A2", "EDN3", "ACER3", "SPATS1","GC", "PLCE1", "COL6A2", "COL3A1", "CRP", "KRT18", "ALDOB", "LCN2", "ACP5", "APOE"))
a <- VlnPlot(df, features = c("SP100", "PRRG3", "DLK1","INS","FN1","COL1A1","REG1B","REG1A","TTR","GCG","CXCL1","KRT19",'EEF1A2','EDN3','ACER3','AL022322.2','GC','PLCE1', 'COL1A2','SPARC','CRP','KRT18','ALDOB','PRSS21','ACP5','APOE'))
a
dev.off()