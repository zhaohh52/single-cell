
#安装和加载所需包
BiocManager::install("Seurat") 
BiocManager::install("dplyr")
BiocManager::install("patchwork")
library(dplyr)
library(Seurat)
library(patchwork) 

#导入数据 
setwd("F:/myworkspace/single") #设置工作环境到数据所在文件夹
counts.data <- Read10X(data.dir = "pbmc3k") # 设置项目名称
#创建Seurat对象,#过滤检测少于200个基因的细胞（min.features = 200）和少于3个细胞检测出的基因（min.cells = 3）
srat <- CreateSeuratObject(counts = counts.data, project = "pbmc3k", min.cells = 3, min.features = 200)


#QC,每个细胞基因数和总分子数在建立seurat对象时就已经自动计算好了
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
head(srat@meta.data, 5)
#nFeature_RNA代表每个细胞测到的基因数目。nCount_RNA代表每个细胞测到所有基因的表达量之和。percent.mt代表测到的线粒体基因的比例。
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
#nCount_RNA与nFeature_RNA的相关性
plot2 <- FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 #合并两图

srat <- subset(srat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#选取 2500 > nFeature_RNA >200 和percent.mt < 5的数据



#数据标准化,若要使用默认值srat <- NormalizeData(srat)
srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)


#鉴定高变基因,每个数据集返回2000个features 。这些将用于下游分析，如PCA
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(srat), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(srat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) # 显示top10的基因
plot1 + plot2


#数据缩放，在PCA降维之前的一个标准预处理步骤
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)


#线性降维
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
#print(srat[["pca"]], dims = 1:5, nfeatures = 5)
#查看PCA结果
VizDimLoadings(srat, dims = 1:2, reduction = "pca")
DimPlot(srat, reduction = "pca")
DimHeatmap(srat, dims = 1, cells = 500, balanced = TRUE) #1个PC 500个细胞
DimHeatmap(srat, dims = 1:15, cells = 500, balanced = TRUE)#15个PC 500个细胞


#确定数据的维度
#1. JackStraw方式
srat <- JackStraw(srat, num.replicate = 100)
srat <- ScoreJackStraw(srat, dims = 1:20)
JackStrawPlot(srat, dims = 1:15)#在10-12个PC之后，显著性大幅下降，也就是前10-12个维度包含了大部分的样本信息。
#2. Elbow plot
ElbowPlot(srat)#PC9-10附近有一个拐点（“elbow”），这表明大部分真实信号是在前10个pc中捕获的。
#综上，选择10个主成分为参数作为后续分析



#细胞聚类,Seurat使用KNN算法进行聚类。
srat <- FindNeighbors(srat, dims = 1:10)#dims = 1:10 即选取前10个主成分来分类细胞。
srat <- FindClusters(srat, resolution = 0.5)
head(Idents(srat), 5)#查看前5个细胞的分类ID


#非线性降维（UMAP/tSNE）
srat <- RunUMAP(srat, dims = 1:10)
DimPlot(srat, reduction = "umap")
DimPlot(srat, reduction = "umap", label = TRUE)# 显示在聚类标签
srat <- RunTSNE(srat, dims = 1:10)# 使用TSNE聚类
DimPlot(srat, reduction = "tsne")
DimPlot(srat, reduction = "tsne", label = TRUE)# 显示在聚类标签
saveRDS(srat, file = "pbmc_tutorial.rds")#保存rds，用于后续分析


#找差异表达基因(聚类标志cluster biomarkers)FindMarkers函数
#cluster1.markers <- FindMarkers(srat, ident.1 = 1, min.pct = 0.25)# cluster 1的标记基因
#head(cluster1.markers, n = 5)
#cluster5.markers <- FindMarkers(srat, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)#找出区分cluster 5与cluster 0和cluster 3的所有标记
#head(cluster5.markers, n = 5)
srat.markers <- FindAllMarkers(srat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
srat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#VlnPlot(srat, features = c("MS4A1", "CD79A")) #可视化
#VlnPlot(srat, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(srat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))#可视化
top10 <- srat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#每个聚类前10个差异基因表达热图(如果小于10，则绘制所有标记)
DoHeatmap(srat, features = top10$gene) + NoLegend()

#鉴定细胞类型
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(srat)
srat <- RenameIdents(srat, new.cluster.ids)
DimPlot(srat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(srat, file = "pbmc3k_final.rds")

#计算细胞周期，使用CellCycleScoring函数计算细胞周期评分 
s.genes <- cc.genes$s.genes # 获取S期marker基因 
g2m.genes <- cc.genes$g2m.genes # 获取G2M期marker基因 
srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 
head(srat) 
RidgePlot(srat, features = c("GMNN", "HMGB2", "NKG7", "NUSAP1"), ncol = 2) #可视化细胞周期标志物在各个地区的分布
