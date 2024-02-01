suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(stringr)
  library(Matrix)
  library(scales)
  library(cowplot)
  library(RCurl)
  library(harmony)
  library(DoubletFinder)
})

rds = 'RDS/All_Merged.rds'
dt <- readRDS(file = rds)

count_file = dt@assays$RNA@counts
dt<-CreateSeuratObject(count=count_file)

dt$sample = dt$orig.ident
# Add number of genes per UMI for each cell to metadata
dt$log10GenesPerUMI <- log10(dt$nFeature_RNA) / log10(dt$nCount_RNA)
dt$mitoRatio <- PercentageFeatureSet(object = dt, pattern = "^MT-")
dt$mitoRatio <- dt@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- dt@meta.data
head(metadata)
metadata$type <- 'tumor'
for (i in c("^S01",'^S11','^S21','^S23','^S25','^S31','^S35','^S04','^S06','^S09')){
  metadata$type[which(str_detect(metadata$sample, i))] <- "normal"}
for (i in c('^S38','^S19','^S20')){
  metadata$type[which(str_detect(metadata$sample, i))] <- "metastasis"}
for (i in c('^S37')){
  metadata$type[which(str_detect(metadata$sample, i))] <- "peritoneum"}


dt@meta.data = metadata

# Visualize the number of cell counts per sample
VlnPlot(dt, features = c("nFeature_RNA"))
VlnPlot(dt, features = c('nCount_RNA'))
FeatureScatter(dt, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

metadata %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)+
  scale_x_continuous(limits = c(0,0.3))

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Extract counts
counts <- GetAssayData(object = dt, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

## QC DONE
dt <- CreateSeuratObject(filtered_counts, meta.data = dt@meta.data)
filtered_counts=NA
counts <-NA
nonzero <- NA


#### Analysis ####
dt <- NormalizeData(dt, normalization.method = "LogNormalize", scale.factor = 10000)
dt <- FindVariableFeatures(dt, selection.method = "vst", nfeatures = 3000)



# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dt), 10)
plot1=VariableFeaturePlot(dt)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

## PCA
all.genes <- rownames(dt)
dt <- ScaleData((dt), features = all.genes)
dt <- RunPCA(dt, features = VariableFeatures(object = dt))


## HARMONY
options(repr.plot.height = 2.5, repr.plot.width = 6)
dt <- RunHarmony(dt, c("sample"))

harmony_embeddings <- Embeddings(dt, 'harmony')



DimPlot(object = dt, reduction = "harmony", pt.size = .1, group.by = "sample")
ElbowPlot(object = dt,  ndims = 50,reduction='harmony')


## Clustering
dim = 40
dt <- RunUMAP(dt, dims = 1:dim,  reduction = "harmony",n.neighbors =30,min.dist = 0.1,n.epochs=300)
DimPlot(dt, reduction = "umap",label=TRUE,raster=FALSE)

dt <- FindNeighbors(dt,reduction = "harmony", dims = 1:dim)
dt <- FindClusters(dt, resolution = 0.1)

DimPlot(dt, reduction = "umap",label=TRUE,raster=FALSE)
#DimPlot(object = dt,  reduction = "umap",group.by = 'RNA_snn_res.0.3',label=TRUE)


gene = c('THY1') #JCHAIN, PTPRC
FeaturePlot(dt, reduction = "umap", 
            features = c(gene), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,raster=FALSE) 

## DoubletFinder
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- dt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(dt@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
dt<- doubletFinder_v3(dt, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(dt)
Idents(object = dt) <- 'DF.classifications_0.25_0.09_11898'
DimPlot(dt, reduction = "umap",label=TRUE,raster=FALSE)

dt <- subset(dt, ident = c('Singlet'))

sub <- subset(dt, ident = c('4','7'))



