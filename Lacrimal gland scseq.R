if(!require(multtest))install.packages("multtest")
## Loading required package: multtest
## Warning: package 'multtest' was built under R version 4.0.3
## Loading required package: BiocGenerics
## Warning: package 'BiocGenerics' was built under R version 4.0.5
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## The following object is masked from 'package:imager':
## 
##     width
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
## Loading required package: Biobase
## Warning: package 'Biobase' was built under R version 4.0.3
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Attaching package: 'Biobase'
## The following object is masked from 'package:imager':
## 
##     channel
if(!require(Seurat))install.packages("Seurat")
## Loading required package: Seurat
## Warning: package 'Seurat' was built under R version 4.0.5
## Registered S3 method overwritten by 'spatstat.geom':
##   method      from  
##   plot.imlist imager
## Attaching SeuratObject
if(!require(dplyr))install.packages("dplyr")
## Loading required package: dplyr
## Warning: package 'dplyr' was built under R version 4.0.5
## 
## Attaching package: 'dplyr'
## The following object is masked from 'package:Biobase':
## 
##     combine
## The following objects are masked from 'package:BiocGenerics':
## 
##     combine, intersect, setdiff, union
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
if(!require(patchwork))install.packages("patchwork")
## Loading required package: patchwork
## Warning: package 'patchwork' was built under R version 4.0.5
if(!require(R.utils))install.packages("R.utils")
## Loading required package: R.utils
## Warning: package 'R.utils' was built under R version 4.0.5
## Loading required package: R.oo
## Warning: package 'R.oo' was built under R version 4.0.3
## Loading required package: R.methodsS3
## Warning: package 'R.methodsS3' was built under R version 4.0.3
## R.methodsS3 v1.8.1 (2020-08-26 16:20:06 UTC) successfully loaded. See ?R.methodsS3 for help.
## R.oo v1.24.0 (2020-08-26 16:11:58 UTC) successfully loaded. See ?R.oo for help.
## 
## Attaching package: 'R.oo'
## The following object is masked from 'package:R.methodsS3':
## 
##     throw
## The following object is masked from 'package:magrittr':
## 
##     equals
## The following objects are masked from 'package:methods':
## 
##     getClasses, getMethods
## The following objects are masked from 'package:base':
## 
##     attach, detach, load, save
## R.utils v2.11.0 (2021-09-26 08:30:02 UTC) successfully loaded. See ?R.utils for help.
## 
## Attaching package: 'R.utils'
## The following object is masked from 'package:magrittr':
## 
##     extract
## The following object is masked from 'package:utils':
## 
##     timestamp
## The following objects are masked from 'package:base':
## 
##     cat, commandArgs, getOption, inherits, isOpen, nullfile, parse,
##     warnings
rm(list = ls())


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
#人源的需要换成MT
head(eye@meta.data,5)
VlnPlot(eye, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(eye, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(eye, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2)) 
## Warning: CombinePlots is being deprecated. Plots should now be combined using
## the patchwork system.
plot2
eye<- subset(eye, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
VlnPlot(eye, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
##ncol(as.data.frame(eye[["RNA"]]@counts))
memory.size()
##计算 分群
eye <- NormalizeData(eye, normalization.method = "LogNormalize", scale.factor = 10000)
#CLR、RC
#normalizes the feature expression measurements for each cell by the total expression
#eye[["RNA"]]@data
eye <- FindVariableFeatures(eye, selection.method = "vst", nfeatures = 2500)
#for PCA DoHeatmap
top10 <- head(VariableFeatures(eye), 10)
plot1 <- VariableFeaturePlot(eye)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
## When using repel, set xnudge and ynudge to 0 for optimal results
plot1 + plot2
# plot variable features with and without labels
## Warning: Transformation introduced infinite values in continuous x-axis
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Transformation introduced infinite values in continuous x-axis
## Warning: Removed 1 rows containing missing values (geom_point).
eye <- ScaleData(eye, features = rownames(eye))
##eye <- ScaleData(eye)
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
#http://biocc.hrbmu.edu.cn/CellMarker/
