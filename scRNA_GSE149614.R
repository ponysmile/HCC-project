#######01.data load ######
ct=fread("GSE149614_RAW/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",data.table = F)
ct[1:4,1:4]
rownames(ct)=ct[,1]
ct=ct[,-1]
sce=CreateSeuratObject(counts =  ct ,
                       min.cells = 5,
                       min.features = 300)
sce$group = ifelse(grepl("T",sce$orig.ident),'T',
                   ifelse(grepl("N",sce$orig.ident),'N',
                          ifelse(grepl("P",sce$orig.ident),'P','L')))

table(sce$group)
table(sce$orig.ident)
scemhh <- sce
library(stringr)
Idents(scemhh) <- scemhh$group
scemhh[["percent.mt"]] <- PercentageFeatureSet(scemhh, pattern = "^MT-")
scemhh[['percent.ribo']] <- PercentageFeatureSet(scemhh, pattern = "^RP[SL]") 
scemhh[['percent.hb']] <- PercentageFeatureSet(scemhh, pattern = '^HB[^(P)]') 
scemhh <- subset(scemhh,idents = c('N','T'))
library(dplyr)
scemhh <- scemhh  %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 1e4)  
scemhh <- scemhh  %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = 2000)  
scemhh <- scemhh  %>%
  ScaleData(vars.to.regress = "percent.mt")  
scemhh <- scemhh  %>%
  RunPCA(features = VariableFeatures(object = scemhh))   
ElbowPlot(scemhh)
scemhh <- scemhh  %>%
  RunUMAP(dims = 1:20)   
scemhh <- scemhh  %>%
  FindNeighbors(dims = 1:20)
temp <- scemhh  %>%
  FindClusters(resolution=seq(0.1,1,0.1))  
library(clustree)
clustree(temp)
scemhh <- scemhh %>% 
  FindClusters(resolution = 0.1)
scemhh <- RunUMAP(scemhh, dims = 1:10)
DimPlot(scemhh, reduction = 'umap',label=T)
DimPlot(scemhh, reduction = 'umap',group.by = 'orig.ident')

#######02.harmony#######
scemhh <- RunHarmony(scemhh,reduction = "pca",group.by.vars = "group",reduction.save = "harmony")
scemhh <- RunUMAP(scemhh, reduction = "harmony", dims = 1:30,reduction.name = "harmony")
scemhh <- FindNeighbors(scemhh
                        # , reduction = "harmony", dims = 1:5
)
scemhh@meta.data$seurat_clusters <- scemhh@meta.data$RNA_snn_res.0.1
Idents(scemhh) <- "seurat_clusters"

#######03.annotation ######
all_markers <- FindAllMarkers(object = scemhh,only.pos = T)
all_markers$rank <- all_markers$pct.1/all_markers$pct.2
library(dplyr)
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:10) %>%
  ungroup() 

scemhh <- RenameIdents(scemhh,
                       "0"= "Tcell","10"="Fibro",
                       "1"= "Mac",'11'='Hep5',
                       "2"= "Endo",'12'='Hep6',
                       "3"= "Hep1", '13'='Bcell',
                       "4"= "Mac",'14'='Tpro',
                       "5"= "Epi",'15'='Hep7',
                       "6"= "Hep2",'16'='Mac',
                       "7"= "Hep3",'17'='Mac',
                       "8"= "Hep4", '18'='Mac',
                       "9"="Bcell", '19'='Epi'
)
