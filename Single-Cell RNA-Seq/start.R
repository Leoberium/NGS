library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(monocle3)
load(file = 'hse.object.Rdata')

# test
str(sc.object)
sc.object@assays$RNA@counts[1:5, 1:5]
sc.object@active.ident
# test

data <- sc.object
data <- RunPCA(object = data, verbose = FALSE)
data <- RunTSNE(object = data, dims = 1:30, verbose = FALSE)
data <- FindNeighbors(object = data, dims = 1:30, verbose = FALSE)
data <- FindClusters(object = data, resolution = 0.2, verbose = FALSE)
DimPlot(object = data, label = TRUE)
TSNEPlot(object = data)
UMAPPlot(object = data)
PCAPlot(object = data)

FeaturePlot()
Idents()
