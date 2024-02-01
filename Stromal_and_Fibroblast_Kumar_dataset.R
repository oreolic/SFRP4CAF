suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(harmony)
  library(ggplot2)
  library(stringr)
  library(Matrix)
  library(scales)
  library(cowplot)
  library(RCurl)
  library(future)
  library(harmony)
  library(DUBStepR)
  library(DoubletFinder)
  library(Rcpp)
})


sub <- readRDS('Kumar_Clustering.rds')
DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10)#,group.by='type')
Idents(object = sub) <- "type"
sub <- subset(sub, ident = c('tumor','normal'))

#Filtering out low complexity (log10GenesPerUMI<0.8) and high mitochondria
sub <- subset(x = sub,
              subset= 
                (log10GenesPerUMI > 0.80) & 
                (mitoRatio < 0.2))

cell_number = 0.01
k_number = 10
pcs =30
sub<-CreateSeuratObject(count=sub@assays$RNA@counts,meta.data = sub@meta.data)
sub <- NormalizeData(object = sub, normalization.method = "LogNormalize")

#DubStepR for variant feature
dubstepR.out <- DUBStepR(input.data = sub@assays$RNA@data, min.cells = cell_number*ncol(sub), optimise.features = T, k = k_number, num.pcs = pcs, error = 0)

sub@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
summary(dubstepR.out$optimal.feature.genes)

sub <- ScaleData(sub, features = rownames(sub))
sub <- RunPCA(sub, features = VariableFeatures(object = sub))

options(repr.plot.height = 2.5, repr.plot.width = 6)
sub <- sub %>% 
  RunHarmony(c("sample"), plot_convergence = TRUE)


harmony_embeddings <- Embeddings(sub, 'harmony')
ElbowPlot(object = sub,  ndims = 50,reduction ='harmony')


dim  =50
sub <- RunUMAP(sub, dims = 1:dim, reduction='harmony',n.components = 2 ,n.neighbors = 50 ,min.dist =0.01,n.epochs=300)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 5)

sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:dim)
sub <- FindClusters(sub, resolution =0.2)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10,group.by='type')


gene = c('JCHAIN') #JCHAIN, PTPRC
FeaturePlot(sub, reduction = "umap",
            features = c(gene), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

VlnPlot(sub, features = c(gene),pt.size = 0)#

sub <- subset(sub, subset = seurat_clusters != "4") #filtering out JCHAIN (immune)
sub <- subset(sub, subset = seurat_clusters != "6") # filtering out PTPRC (immune)
sub <- subset(sub, subset = seurat_clusters != "7") # filtering out PTPRC (immune)
# REHARMONY
options(repr.plot.height = 2.5, repr.plot.width = 6)
sub <- sub %>% 
  RunHarmony(c("sample"), plot_convergence = TRUE)

head(sub)

harmony_embeddings <- Embeddings(sub, 'harmony')
ElbowPlot(object = sub,  ndims = 50,reduction ='harmony')


dim  = 20
sub <- RunUMAP(sub, dims = 1:dim, reduction='harmony',n.components = 2 ,n.neighbors = 50 ,min.dist =0.1,n.epochs=300)
#sub<-RunTSNE()

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 5)

sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:dim)
sub <- FindClusters(sub, resolution =0.1)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10)


gene = c('CXCL14') #JCHAIN, PTPRC
FeaturePlot(sub, reduction = "umap", 
            features = c(gene), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
VlnPlot(sub, features = c(gene),pt.size = 0)#


sub <- RenameIdents(object = sub, 
                    '0'='Fibroblast',
                    '1'= 'Mesothelial',
                    '2'= 'Pericyte',
                    '3'= 'SMC',
                    '4'='Proliferation',
                    '5'= 'apCAF'
)

umap = DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10)


## Fibroblast Analysis
DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10)


### SELECT FIBROBLAST ONLY
sub <- subset(sub, ident = c('Fibroblast'))

cell_number = 0.01
k_number = 10
pcs =30
sub<-CreateSeuratObject(count=sub@assays$RNA@counts,meta.data = sub@meta.data)
sub <- NormalizeData(object = sub, normalization.method = "LogNormalize")


dubstepR.out <- DUBStepR(input.data = sub@assays$RNA@data, min.cells = cell_number*ncol(sub), optimise.features = T, k = k_number, num.pcs = pcs, error = 0)

#dubstepR.out = readRDS('All_Pipeline/CancerNormal_Fibroblast_Only_DubStepR_0921.rds')
saveRDS(dubstepR.out,'All_Pipeline/CancerNormal_Fibroblast_Only_DubStepR_0925.rds')



sub@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
summary(dubstepR.out$optimal.feature.genes)


sub <- ScaleData(sub, features = rownames(sub))
sub <- RunPCA(sub, features = VariableFeatures(object = sub))

options(repr.plot.height = 2.5, repr.plot.width = 6)
sub <- sub %>% 
  RunHarmony(c("sample"), plot_convergence = TRUE)

head(sub)

harmony_embeddings <- Embeddings(sub, 'harmony')


ElbowPlot(object = sub,  ndims = 50,reduction ='harmony')


dim  = 20
sub <- RunUMAP(sub, dims = 1:dim, reduction='harmony',n.components = 2 ,n.neighbors = 12 ,min.dist =0.2,n.epochs=300)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 5)

sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:dim)
sub <- FindClusters(sub, resolution =0.3)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10)

DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 10,group.by='type')

gene = c('SFRP4') #JCHAIN, PTPRC
FeaturePlot(sub, reduction = "umap",
            features = c(gene), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


markers <- FindAllMarkers(object = sub, only.pos = TRUE, logfc.threshold = 0.5)    

VlnPlot(sub, features = c(gene),pt.size = 0)#
DotPlot(sub,features=gene)

sub <- RenameIdents(object = sub, 
                    '0'='CXCL14POSTN',
                    '1'= 'CXCL14APOE',
                    '2'= 'CXCL14APOE',
                    '3'= 'SFRPfour',
                    '4'= 'SFRPtwo')


DimPlot(object = sub, reduction = "umap") 
umap <- DimPlot(object = sub, reduction = "umap") 
