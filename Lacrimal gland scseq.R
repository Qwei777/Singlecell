rm(list = ls())

install.packages('Seurat')
install.packages('dplyr')
update.packages('dplyr')
install.packages("mindr")
install.packages("tidyverse")
library(Seurat)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(mindr)
library(tidyverse)
list.files("filtered_feature_bc_matrix")
eye.data <- Read10X(data.dir = "filtered_feature_bc_matrix")
eye <- CreateSeuratObject(counts = eye.data, project = "lacrimal gland")
eye

ncol(eye)
eye[["percent.mt"]] <- PercentageFeatureSet(eye, pattern = "^mt-")
head(eye@meta.data,5)
VlnPlot(eye, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(eye, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(eye, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2)) 
plot2
eye<- subset(eye, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
VlnPlot(eye, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
##ncol(as.data.frame(eye[["RNA"]]@counts))用上一步方法检测即可，总是报错
memory.size()
##计算 分群
eye <- NormalizeData(eye, normalization.method = "LogNormalize", scale.factor = 10000)
eye <- FindVariableFeatures(eye, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(eye), 10)
plot1 <- VariableFeaturePlot(eye)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
###eye <- ScaleData(eye, features = rownames(eye))
eye <- ScaleData(eye)
eye <- RunPCA(eye, features = VariableFeatures(object = eye))
print(eye[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(eye, dims = 1:2, reduction = "pca")
DimPlot(eye, reduction = "pca")
DimHeatmap(eye, dims = 1:15, cells = 500, balanced = TRUE)
eye <- JackStraw(eye, num.replicate = 100)
eye <- ScoreJackStraw(eye, dims = 1:20)
JackStrawPlot(eye, dims = 1:20)  
ElbowPlot(eye) 
eye <- FindNeighbors(eye, dims = 1:14)
eye <- FindClusters(eye, resolution = 0.8)
eye <- RunUMAP(eye, dims = 1:10)
DimPlot(eye, reduction = "umap",label = T)

###eye <- RunTSNE(eye, dims = 1:10)
###DimPlot(eye, reduction = "tsne")
eyemarkers <- FindAllMarkers(eye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#if(!require(dplyr))install.packages("dplyr")
eyemarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
VlnPlot(eye, features = c("Cd3e"), slot = "counts", log = TRUE)
FeaturePlot(eye, features = c("Eomes"))

top10 <- eyemarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(eye, features = top10$gene) + NoLegend()

write.table(eye.markers,  # 要导出的数据
            file = "mtcars.txt", # 指定导出数据的文件名称和格式
            sep = "\t", # 字段分隔符
            row.names = TRUE, # 逻辑词，是否将行名一起导出
            col.names = NA) # 逻辑词，是否将列名一起导出



saveRDS(eye,'Desktop/eye.rds')
eye<- readRDS('eye.rds')
##png(filename="p.png",width=1200,height=900,res=300)
##print(p)
##dev.off()
##三个文件eye,eye.markers,top10
#http://biocc.hrbmu.edu.cn/CellMarker/