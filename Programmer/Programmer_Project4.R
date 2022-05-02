#Programmer -Shreenk.Project BF528: Applications in Translational Bioinformatics, Project 4 

library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
# use the count matrix to create a Seurat object.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

x <- file.path("/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant/alevin/quants_mat.gz")  
txi <- tximport(x, type="alevin")
v<- txi$counts

# Creating Seurat Object
pbmc <- CreateSeuratObject(counts = v , min.cells = 3, min.features = 200, project = "10X_PBMC")

# Data QC - Create Mitochondrial Percent
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Data QC - nFeatures, nCount and PercentMT
png("dataqc_vln1.png")
a<-VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
a


png("dataqc_scatter.png", width = 800)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
b<-plot1 + plot2
dev.off()
b


cells_subset <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)

#Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#Identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labelsVizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE,xnudge =0,ynudge= 0)
plot1 + plot2


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#VISULISATION OF DATA USING VIMDIM
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")


#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)



pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# Save RDS
saveRDS(cells_subset, file = "programmer_output.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers)


# Cells Per Cluster
cluster_numbers <- table(pbmc@meta.data$seurat_clusters)
cluster_numbers <- as.data.frame(cluster_numbers)
cluster_numbers %>% ggplot(mapping = aes(x = Var1, y = Freq)) +
  geom_bar(stat="identity") +
  ggtitle("Distribution of Cells Across Clusters") +
  xlab("Cluster ID") +
  ylab("Number of Cells")