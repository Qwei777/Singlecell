if(T){rm(list = ls())

############1.Database clustering
library(dplyr)
top10 <- eye.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(eye, features = top10$gene) + NoLegend()
VlnPlot(eye, features = top10$gene[1:20])
VlnPlot(eye, features = top10$gene[1:20],pt.size=0)
DimPlot(eye,label = T)
bfreaname.eye <- eye
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(eye)
eye <- RenameIdents(eye, new.cluster.ids)
DimPlot(eye, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
###2. SingleR clustering
if(!require(singleR))BiocManager::install(singleR)
BiocManager::install('SingleR')
library(SingleR)

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("SingleR")

library(SingleR)

## Loading required package: SingleR
## Loading required package: SummarizedExperiment
## Loading required package: MatrixGenerics
## Loading required package: matrixStats
## Warning: package 'matrixStats' was built under R version 4.0.5
## 
## Attaching package: 'matrixStats'
## The following object is masked from 'package:dplyr':
## 
##     count
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
## 
## Attaching package: 'MatrixGenerics'
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
## The following object is masked from 'package:Biobase':
## 
##     rowMedians
## Loading required package: GenomicRanges
## Loading required package: stats4
## Loading required package: S4Vectors
## 
## Attaching package: 'S4Vectors'
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
## The following object is masked from 'package:base':
## 
##     expand.grid
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## The following object is masked from 'package:R.oo':
## 
##     trim
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
## The following object is masked from 'package:grDevices':
## 
##     windows
## Loading required package: GenomeInfoDb
## Warning: package 'GenomeInfoDb' was built under R version 4.0.5
## 
## Attaching package: 'SummarizedExperiment'
## The following object is masked from 'package:SeuratObject':
## 
##     Assays
## The following object is masked from 'package:Seurat':
## 
##     Assays
rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/")) 
options(download.file.method = 'libcurl')
options(url.method='libcurl')
BiocManager::install("SingleR")
###
remove.packages('matrixStats')#移除后可能需要重新打开
if(!require(matrixStats))BiocManager::install('matrixStats')
# remove.packages(c('dplyr','ellipsis'))
# install.packages(c('dplyr','ellipsis'))
install.packages("SingleR")
# remove.packages(c('vctrs'))
# install.packages(c('vctrs'))
if(!require(celldex))
BiocManager::install('celldex')#There are some version conflicts and some packages need to be reinstalled。
## Loading required package: celldex
## 
## Attaching package: 'celldex'
## The following objects are masked from 'package:SingleR':
## 
##     BlueprintEncodeData, DatabaseImmuneCellExpressionData,
##     HumanPrimaryCellAtlasData, ImmGenData, MonacoImmuneData,
##     MouseRNAseqData, NovershternHematopoieticData
###Built-in data sets：
# cg=BlueprintEncodeData()
# cg=DatabaseImmuneCellExpressionData()
# cg=NovershternHematopoieticData()
# cg=MonacoImmuneData()
# cg=ImmGenData()
# cg=MouseRNAseqData()
# cg=HumanPrimaryCellAtlasData()

install.packages('celldex')

BiocManager::install("Biostrings")

install.packages("remotes")

install.packages("Biostrings")
library(Biostrings)
remotes::install_github("LTLA/celldex")

install.packages('SeuratObject')
library(SingleR)
library(celldex)
cg<-ImmGenData()

celldex::ImmGenData()

#Select the immune cell reference data set we want to utilize
## snapshotDate(): 2020-10-27
## see ?celldex and browseVignettes('celldex') for documentation
## Could not check id: EH3494 for updates.
##   Using previously cached version.
## loading from cache
## Could not check id: EH3494 for updates.
##   Using previously cached version.
## see ?celldex and browseVignettes('celldex') for documentation
## Could not check id: EH3495 for updates.
##   Using previously cached version.
## loading from cache
## Could not check id: EH3495 for updates.
##   Using previously cached version.

assay_for_SingleR <- GetAssayData(tcell, slot="data")
#Remove the expressed sequence from the sample
predictions <- SingleR(test=assay_for_SingleR, 
                       ref=cg, labels=cg$label.main)

table(predictions$labels)#Let's see what cells were annotated.

cellType=data.frame(seurat=tcell@meta.data$seurat_clusters,
                    predict=predictions$labels)#Get the relationship between numbering and predicted labels in seurat
sort(table(cellType[,1]))


table(cellType[,1:2])
#Accessing 1~2 columns of celltyple

#singleR without a suitable dataset, the results obtained can only be used as an auxiliary reference
