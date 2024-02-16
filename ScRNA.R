######PACKAGES###################
library(Seurat)
library(patchwork)
library(hdf5r)
library(tidyverse)
library(Neighborseq)
library(RColorBrewer)
library(Neighborseq)
library(RColorBrewer)
library('celldex')
library('SingleR')

pbmc.cancerData <- Read10X("hg19")
pbmc<- CreateSeuratObject(pbmc.cancerData)
#str(pbmc)
#pbmc.cancerData[1:10,1:10]
pbmc[["percent.mt"]] <-PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc,features = c("nFeature_RNA",
                          "nCount_RNA",
                          "percent.mt"),ncol = 3)
#plot1<- FeatureScatter(pbmc,feature1 ="nCount_RNA",
                        #feature2 = "percent.mt")
#plot2<- FeatureScatter(pbmc,feature1 ="nCount_RNA", 
                        #feature2 = "nFeature_RNA")
#plot1 + plot2        
pbmc<-subset(pbmc,subset = nFeature_RNA >200 & nFeature_RNA <2500
             & percent.mt < 5)
pbmc<-NormalizeData(pbmc, normalization.method = "LogNormalize",
                    scale.factor = 10000)
pbmc<-FindVariableFeatures(pbmc,selection.method = "vst",
                           nfeatures = 2000)
#top20<-head(VariableFeatures(pbmc),10)
#plot3<- VariableFeaturePlot(pbmc)
#plot4<- LabelPoints(plot = plot3,points = top20)
#plot3 + plot4
pbmc<-ScaleData(pbmc, features = rownames(pbmc))
pbmc<-RunPCA(pbmc)
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#DimPlot(pbmc, reduction = "pca")
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
#pbmc<-JackStraw(pbmc, num.replicate = 100)
#pbmc<-ScoreJackStraw(pbmc, dims = 1:20)
#JackStrawPlot(pbmc, dims=1:15)
#ElbowPlot(pbmc)
pbmc<-FindNeighbors(pbmc,dims = 1:10)
pbmc<-FindClusters(pbmc,resolution = 0.5)
#head(Idents(pbmc), 5)
#sid_cluster<-Idents(pbmc)
#summary(sid_cluster)
pbmc<-RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc,reduction = "umap")
#saveRDS(pbmc, file = "../pbmc_tutorial.rds")
#cluster2<-FindMarkers(pbmc,ident.1=2, min.pct = 0.25)
#head(cluster2, 10)
# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5,ident.2 = c(0,3),min.pct = 0.25)
#head(cluster5.markers, n = 5)
#allMarkers<- FindAllMarkers(pbmc, min.pct = 0.25)
#library(tidyverse)
#allMarkers %>% group_by(cluster) %>% head(30) %>% print(n=25)
#allMarkers %>% group_by(cluster) %>% slice_max(n=5, 
                                               #order_by =avg_log2FC) %>% print(n=45)
#VlnPlot(pbmc, features = c("CD79A", "CD79B", "MS4A1", "CD19"), log = T)
#FeaturePlot(pbmc, features =c("CD79A", "CD79B", "MS4A1", "CD19", "CD3D", "NKG7","PRF1", "GZMA", 
                              #"CD14", "CD68", "CD8A") )
#allMarkers %>% group_by(cluster) %>% top_n(n=10, wt= avg_log2FC) -> top10
#DoHeatmap(pbmc, features = top10$gene) + NoLegend()
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc$cell_type<- RenameIdents(pbmc, new.cluster.ids)
pbmc$singleR<- Testing_main
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = Testing_main$labels) + NoLegend()
saveRDS(pbmc, file = "../pbmc3k_final.rds")
pbmc$cell_identity<-Idents(pbmc)
pbmc@meta.data

###cancer cancerDataa

cancerData<- readRDS('~/Downloads/breast_cancer_scRNA')
#str(cancerData)
cancerData <- Read10X_h5('~/Downloads/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
cancerData1<- cancerData$`Gene Expression`
cancerData<- CreateSeuratObject(counts = cancerData1, min.cells = 3, min.features = 100)
str(cancerData)
#cancerData@assays$RNA
#cancerData@meta.data
#cancerData[1:10,1:10]
#table(cancerData$sample)
#table(cancerData$patient, cancerData$tissue)
cancerData[["percent.mt"]] <-PercentageFeatureSet(cancerData, pattern = "^MT-")
VlnPlot(cancerData,features = c("nFeature_RNA",
                          "nCount_RNA",
                          "percent.mt"),ncol = 3)
plot1<- FeatureScatter(cancerData,feature1 ="nCount_RNA",
                       feature2 = "percent.mt")
plot2<- FeatureScatter(cancerData,feature1 ="nCount_RNA", 
                       feature2 = "nFeature_RNA")
view(cancerData@meta.data)
plot1 + plot2        
#VlnPlot(cancerData, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = 'tissue')
ggplot(cancerData@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
  geom_point(size = 0.8) +
  scale_colour_viridis_c(option = 'F') + 
  lims(x = c(0, NA), y = c(0, NA)) +
  facet_wrap(~sample, nrow = 5, scales = 'free') +
  theme_minimal()
cancerData<-subset(cancerData,subset = nFeature_RNA >300 & nFeature_RNA <2500
             & percent.mt < 5)
cancerData<-NormalizeData(cancerData, normalization.method = "LogNormalize",
                    scale.factor = 10000)
cancerData<-FindVariableFeatures(cancerData,selection.method = "vst",
                           nfeatures = 2000)
#top20<-head(VariableFeatures(cancerData),20)
#plot3<- VariableFeaturePlot(cancerData)
#plot4<- LabelPoints(plot = plot3,points = top20)
#plot3 + plot4
cancerData<-ScaleData(cancerData, features = rownames(cancerData))
cancerData<-RunPCA(cancerData, features = VariableFeatures(object = cancerData))
#print(cancerData[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(cancerData, dims = 1:2, reduction = "pca")
DimPlot(cancerData, reduction = "pca")
#DimHeatmap(cancerData, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(cancerData, dims = 1:15, cells = 500, balanced = TRUE)
#cancerData<-JackStraw(cancerData, num.replicate = 100)
#cancerData<-ScoreJackStraw(cancerData, dims = 1:20)
#JackStrawPlot(cancerData, dims=1:20)
#ElbowPlot(cancerData)
cancerData<-FindNeighbors(cancerData,dims = 1:10)
cancerData<-FindClusters(cancerData,resolution= 0.3)
#cancerData<-FindClusters(cancerData,resolution =c(0.3,0.5,0.7,1))
#head(Idents(cancerData), 5)
#sid<-Idents(cancerData)
#summary(sid)
cancerData<-RunUMAP(cancerData, dims = 1:10)
DimPlot(cancerData,reduction = "umap", label=T)
#DimPlot(cancerData,group.by = 'RNA_snn_res.0.3')
#DimPlot(cancerData, reduction = 'umap', group.by = c('sample', 'patient', 'tissue', 'cell_type_major'), ncol = 2)
#saveRDS(cancerData, file = "../cancerData_tutorial.rds")
cluster2<-FindMarkers(cancerData,ident.1=2, min.pct = 0.25)
head(cluster2, 30)
# find all markers distinguishing cluster 5 from clusters 0 and 3
#clusterDIFF.markers <- FindMarkers(cancerData, ident.1 = 8,ident.2 = 2,min.pct = 0.25)
#head(clusterDIFF.markers, n = 30)
#allMarkers<- FindAllMarkers(cancerData, min.pct = 0.25)
#library(tidyverse)
#allMarkers %>% group_by(cluster) %>% head(30) %>% print(n=25)
#allMarkers %>% group_by(cluster) %>% slice_max(n=10, 
                                               #order_by =avg_log2FC) %>% print(n=140)
VlnPlot(cancerData, features = c('CD3D', 'CD79A','CD19', 'CD14'),group.by ='ScType_cell')
#VlnPlot(cancerData, features = c('CST3','CD1C','CD86','CD83','CD80','HLA-DRA','ITGAE','HLA-DMA','HLA-DMB'))
#VlnPlot(cancerData, features = c('MUC5B', 'DUSP4','MUC16','MUC4','MUC20'))
FeaturePlot(cancerData, features =c("CD79A", "CD79B", "MS4A1", "CD19", "CD3D", "NKG7","PRF1", "GZMA", 
                              "CD14", "CD68", "CD8A", "NCAM1", "FCGR3A", "CD4") )
#VlnPlot(cancerData, features =c("CD79A", "CD79B", "MS4A1", "CD19", "CD3D", "NKG7","PRF1", "GZMA", 
                                    #"CD14", "CD68", "CD8A", "NCAM1","FCGR3A","CD4") )
#allMarkers %>% group_by(cluster) %>% top_n(n=10, wt= avg_log2FC) -> top10
#DoHeatmap(cancerData, features = top10$gene) + NoLegend()
#cancerData<- RenameIdents(cancerData, "5" = "Mono-Macrophages")
#cancerData<- RenameIdents(cancerData, "5" = "Mono-Macrophages")
new.cluster.ids <- c('ENDOTHELIAL CELL', 'FIBROBLAST/SMOOTH', 'Mast cells', 'neurofibroblast', 'UNKNOWN2', 'cytotoxic B cells', 'Plasma cells', 'alveolar epi cells', 'Mono-Macrophage', 'CD4+ T-cell', 'CD8 T-cell', 'B cells', 'Unknown', 'DC')
names(new.cluster.ids) <- levels(cancerData)
cancerData <- RenameIdents(cancerData, new.cluster.ids)
cancerData$celltype_Sid<- Idents(cancerData)
DimPlot(cancerData, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'celltype_Sid') + NoLegend()
#saveRDS(cancerData, file = "../cancerData3k_final.rds")

########integration#######################################################
cancerData@assays$RNA@data <- cancerData@assays$RNA@counts
cancerData.list<- SplitObject(cancerData, split.by = 'patient')

cancerData.list<- lapply(X=cancerData.list, FUN = function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = cancerData.list)
anchors <- FindIntegrationAnchors(object.list = cancerData.list, anchor.features = features, verbose = F)
cancerData.combined<- IntegrateData(anchorset = anchors)
DefaultAssay(cancerData.combined) <- 'integrated'
cancerData.combined <- ScaleData(cancerData.combined, verbose = FALSE)
cancerData.combined <- RunPCA(cancerData.combined, npcs = 10, verbose = FALSE)
cancerData.combined <- FindNeighbors(cancerData.combined, dims = 1:10)
cancerData.combined <- RunUMAP(cancerData.combined, dims = 1:10)
DimPlot(cancerData.combined, reduction = 'umap', group.by = c('sample', 'patient', 'tissue', 'cell_type_major'), ncol = 2)
DimPlot(cancerData.combined, reduction = 'umap', split.by = 'tissue')
cancerData.combined <- FindClusters(cancerData.combined, resolution = 0.3) 
#DimPlot(cancerData.combined, reduction = 'umap', group.by = c('integrated_snn_res.0.2', 'integrated_snn_res.0.4', 'integrated_snn_res.0.6', 'integrated_snn_res.0.8', 'integrated_snn_res.1', 'integrated_snn_res.1.2'), label = T, repel = T, ncol = 2)
#X11()
install.packages("clustree")
library(clustree)
#clustree(cancerData.combined, prefix = "integrated_snn_res.")
#cancerData.combined@meta.data
#plot_df <- cancerData.combined@meta.data %>% pivot_longer(cols = starts_with('integrated_snn_res.'), names_to = 'res')
#ggplot(plot_df, aes(x = value, fill = cell_type_major)) +
 # geom_bar(position = 'fill') +
  #scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  #facet_wrap(~res, scales = 'free') +
  #theme_minimal()
p1 <- DimPlot(cancerData.combined, reduction = 'umap', label = T, repel = T) + NoLegend()
#p2 <- DimPlot(cancerData.combined, reduction = 'umap', group.by = 'cell_type_minor', label = T, repel = T) + NoLegend()
#p1+p2
markers <- list(`NK cells` = c('NCAM1', 'NCR1', 'KLRK1'),
                `Cytotoxic T/\nNK cells` = c('GNLY', 'PFN1', 'GZMA', 'GZMB', 'GZMM', 'GZMH'),
                `Exhausted T/\nTregs` = c('FOXP3', 'CTLA4', 'TIGIT', 'TNFRSF4', 'LAG3', 'PDCD1'),
                `T cells`  = c('CD8A', 'CD3E', 'CD4'),
                `Naive\nT cells` = c('IL7R'),
                `B\ncells` = c('CD19'),
                `Mast\ncells` = c('ENPP3', 'KIT'),
                `pDC` = c('IL3RA', 'LILRA4'),
                `Monocytic\nlineage` = c('HLA.DRA', 'FCGR3A', 'CD68', 'ANPEP', 'ITGAX', 'CD14', 'ITGAM', 'CD33'))
DotPlot(cancerData.combined, features = markers, group.by = 'cell_type_major', assay = 'RNA')+
 scale_colour_viridis_c(option = 'H') +
  theme(axis.text = element_text(size = 8),
       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = 'bottom',
        strip.text = element_text(size = 7)) +
  coord_cartesian(clip = 'off')

 DefaultAssay(cancerData.combined) <- "RNA"
marker_6 <- FindConservedMarkers(cancerData.combined, ident.1 = 6, grouping.var = "tissue", verbose = FALSE)
head(marker_6)
#install.packages('BiocManager')
#BiocManager::install('fftw3')
#install.packages('metap')
#install.packages("qqconf")
#install.packages("fftw")
library(multtest)
library(metap)
library(SeuratData)
#InstallData('ifnb')
#LoadData('ifnb')
Idents(cancerData.combined)<- cancerData.combined@meta.data$cell_type_minor
Idents(cancerData.combined)
DimPlot(cancerData.combined, reduction = 'umap', label = T)
cancerData.combined$celltyp_diff<- paste0(cancerData.combined$cell_type_minor,'_',cancerData.combined$tissue)
view(cancerData.combined@meta.data)
Idents(cancerData.combined)<- cancerData.combined$celltyp_diff
summary(Idents(cancerData.combined))
Finded1<-FindMarkers(cancerData.combined, ident.1 = 'T:CD8+EM_NORMAL', ident.2 = 'T:CD8+EM_TUMOR', pos.only=T, min.pct = 0.25)
FeaturePlot(cancerData.combined, features = c("IRF8", "CD96", "IFIT1"), split.by = "tissue", max.cutoff = 3,
            cols = c("grey", "red"), label = F, repel = T)
plots <- VlnPlot(cancerData.combined, features = c("IRF8", "PRF1", "GZMA"), split.by = "tissue", group.by = "cell_type_major",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)
#X11()
library(cowplot)
theme_set(theme_cowplot())
Idents(cancerData.combined)<-cancerData.combined$cell_type_minor 
TCELL <- subset(cancerData.combined, idents = 'T:CD8+EM')
Idents(TCELL) <- "tissue"
  avg.T_CELL <- as.data.frame(log1p(AverageExpression(TCELL, verbose = FALSE)$RNA))
 avg.T_CELL$gene <- rownames(avg.T_CELL)
 #points.to.label<- head(desc(avg.T_CELL$gene), 20)
 genes.to.label = c("S100A4","ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8", "IFNG", "GNLY", "GZMA")
 p1 <- ggplot(avg.T_CELL, aes(TUMOR,NORMAL)) + geom_point() + ggtitle("T CELLS")
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = T)
p1 
######SCTransform method########################
ifnb.list <- SplitObject(cancerData, split.by = "patient")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "tissue")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "cell_type_major", label = TRUE,
              repel = TRUE)
table(cancerData.combined$cell_type_minor)
x11()
p1 + p2
library(SeuratData)
InstallData("panc8")
LoadData("panc8")
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]
pancreas.list<- lapply(X=pancreas.list, FUN = function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
})
ref.list<- pancreas.list[c("celseq","celseq2","smartseq2", "fluidigmc1")]
features <- SelectIntegrationFeatures(object.list = ref.list)
anchors.pancreas <- FindIntegrationAnchors(object.list = ref.list, anchor.features = features, verbose = F)
ref.combined<- IntegrateData(anchorset = anchors.pancreas)
DefaultAssay(ref.combined) <- 'integrated'
ref.combined <- ScaleData(ref.combined, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, npcs = 30, verbose = FALSE)
ref.combined <- FindNeighbors(ref.combined, dims = 1:10)
ref.combined <- RunUMAP(ref.combined, dims = 1:30, reduction = "pca")
p1 <- DimPlot(ref.combined, reduction = "umap", group.by = "tech")
p2 <- DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  NoLegend()
view(ref.combined@meta.data)
x11()
p1 + p2
#####reference query#############
pancreas.query <- pancreas.list[["indrop"]]
DefaultAssay(ref.combined) <- 'integrated'
pancreas.anchors <- FindTransferAnchors(reference = ref.combined, query = pancreas.query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = ref.combined$celltype,
                            dims = 1:30)
view(predictions)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
table(pancreas.query$predicted.id)
table(pancreas.query$celltype)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")
ref.combined <- RunUMAP(ref.combined, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = ref.combined, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
###########automatic annotation####################################
#install.packages("scCATCH")
#library(scCATCH)
#sc_data<- cancerData.combined@assays$RNA@data
#sc_cluster<- factor(cancerData@meta.data$cluster)
#ncol(sc_data)
#length(sc_cluster)
#obj <- createscCATCH(data = sc_data, cluster = 'sc_cluster')
##########################singleR#######
#install.packages('celldex')
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

# The following initializes usage of Bioc devel
#BiocManager::install(version='3.16')

#BiocManager::install("celldex")
browseVignettes("celldex")
library('celldex')
library('SingleR')
sc_data<- cancerData@assays$RNA@data
ref1 <- HumanPrimaryCellAtlasData()
#BiocManager::install('scRNAseq')
#browseVignettes("SingleR")
#install.packages('SingleR')
#BiocManager::install('SingleR')

celltype_hpca_fine<- SingleR(test = sc_data, ref = ref1, assay.type.test=1,
                  labels = ref1$label.fine)
celltype_main<-SingleR(test = sc_data, ref = ref1, assay.type.test=1,
                      labels = ref1$label.main)
ref2<- celldex::ImmGenData()
celltype_imm_main<- SingleR(test = sc_data, ref = ref2, assay.type.test=1,
                  labels = ref2$label.main)
celltype_imm_fine<-SingleR(test = sc_data, ref = ref2, assay.type.test=1,
                      labels = ref2$label.fine)
ref_encode<-BlueprintEncodeData() 
celltype_encode_main<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                          labels = ref_encode$label.main)
celltype_encode_fine<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                           labels = ref_encode$label.fine)
cancerData@meta.data
view(Testing)
view(table(celltype_main$labels))
cancerData$celltype_encode_main<- celltype_encode_main$pruned.labels
cancerData$celltype_encode_fine<- celltype_encode_fine$pruned.labels
cancerData$celltype_imm_main<- celltype_imm_main$pruned.labels
view(cancerData@meta.data)
p1<-DimPlot(cancerData, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = c('celltype_encode_main', 'ScType_cell'), repel = T) + NoLegend()
p3<-DimPlot(cancerData, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = c('ScType_cell','celltype_encode_fine'), repel = T) + NoLegend()
#p4<-DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'blue_fine_cell', repel=T) + NoLegend()
x11()
p1
p3
myplot<- p1+p3

#######scType##############################
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
## DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = cancerData[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(cancerData@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(cancerData@meta.data[cancerData@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(cancerData@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
cancerData@meta.data$ScType_cell = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  cancerData@meta.data$ScType_cell[cancerData@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

p5<-DimPlot(cancerData, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'ScType_cell')     
x11()
p5
pbmc$ScType_cell
p6<-DimPlot(cancerData, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = c('cell_type_major', 'blue_main_cell', 'ScType_cell'), repel = T) + NoLegend()
p7<-DimPlot(cancerData, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = c('blue_fine_cell','ScType_cell'), repel = T) + NoLegend()
p6
p7
VlnPlot(cancerData, features = c('CD3D', 'CD79A'),group.by ='celltype_encode_main', ncol = 2)
VlnPlot(cancerData, features = c('CD3G', 'NKG7','IFNG', 'FCGR3A'),group.by ='ScType_cell', ncol = 2)
VlnPlot(cancerData, features = c('PDCD1', 'HAVCR2','NCAM1', 'JCHAIN'),group.by ='ScType_cell', ncol = 2)
FeaturePlot(cancerData, features =c("CD79A", "CD79B", "MS4A1", "CD19", "CD3D", "NKG7","PRF1", "GZMA", 
                                    "CD14", "CD68", "CD8A", "NCAM1", "FCGR3A", "CD4") )

# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
#install.packages("data.tree")
library(data.tree)
library(igraph)
#install.packages("igraph")
# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#BiocManager::install('scater')
library('SingleCellExperiment')
library('scuttle')
library('scater')
#install.packages('gggr')
#BiocManager::install('gggr')
scater::multiplot(DimPlot(cancerData, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)
# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

# guess a tissue type
tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = tubo_ovarian@assays$RNA@scale.data, scaled = TRUE, assay = "RNA")  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used         

#################CELL COMMUNICATION###########################
     ######## NICHE NET###############
devtools::install_github("saeyslab/nichenetr")
BiocManager::install('ComplexHeatmap')
library('nichenetr')
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5]
##############################################################################################################################
######NEIGHBOUR-SEQ############################
devtools::install_github('sjdlabgroup/Neighborseq')
install.packages("qqconf")
install.packages("fftw")
library(Neighborseq)
library(RColorBrewer)
library(Neighborseq)
library(RColorBrewer)

my_cell<-cancerData@meta.data$cell_type_minor
#data<-subset(pbmc,subset = nFeature_RNA >200)
ns.data = prep_cell_mat(ge = cancerData@assays$RNA,celltypes = my_cell, logfc.threshold = 0.5)
table(ns.data$celltypes)
ns.data$markers
set.seed(0) # for reproducibility
mt = multiplet_types(ns.data$celltypes)
am = artificial_multiplets(cell.mat = ns.data$cell.mat, 
                           celltypes = ns.data$celltypes, 
                           multiplet_classes = mt)
rf = multiplet_rf(am)
#view(rf$test)
#> [1]  train-mlogloss:3.545750
mroc = mroc_format(rf$test, rf$pred)
mroc_plot(mroc)
#p5
#view(cancerData)
#library(ggraph)
pred = xgpred(rf, ns.data$cell.mat)
result = multiplet_pop_sig(pred = pred)
view(result)
plot_interactions(result, 
                  legend.position = 'right', 
                  width_range = c(0.5,1))
view(table(pbmc$seurat_clusters, pbmc$ScType_cell))
#> Loading required package: ggraph
#> Loading required package: ggplot2
set.seed(0) # for reproducibility
ns.data = prep_cell_mat(ge = cancerData@assays$RNA,celltypes = cancerData$cell_type_major, logfc.threshold = 0.5)
haber.ns = neighborseq(ns.data$cell.mat, ns.data$celltypes,iter = 10, do.mroc = T, sample = cancerData$tissue)
ns.data$celltypes
# plot with colored nodes
#color = colorRampPalette((brewer.pal(9,"Greens")))(11)[6:12]; 
#names(color) = c('Paneth','Stem','TA','EP','Enterocyte'); color=c(color, Goblet='gold3')
x11()
result_normal<- haber.ns$combined_result %>% filter(sample== 'NORMAL')
result_tumor<-haber.ns$combined_result %>% filter(sample== 'TUMOR')
p1<-plot_interactions(result_normal, combined = T,
                  legend.position = 'left', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green')

p2<-plot_interactions(result_tumor, combined = T,
                      legend.position = 'right', 
                      width_range = c(0.5,1), min_color = 'red', max_color = 'green')

p1+p2
mroc_plot(haber.ns$mroc[[1]])
to_plot<- (1:10)
for (i in to_plot) {
  my_plot<- mroc_plot(haber.ns$mroc[[i]])
  print(my_plot)
}
neigbor_diff<- FindMarkers(ns.data, ident.1= pred$`B:_MONOCYTE:`)
view(pred$`B:`)
#######################OVARIAN ANALYSIS###############
message("Loading matrix...")
tuboData<- readRDS('~/Downloads/2095-Olbrecht_counts_matrix')
meta <- read_csv("~/Downloads/2093-Olbrecht_metadata.csv")
colnames(tuboData)
message("Metadata per cell...")
cell_names <- data.frame(cell_label = colnames(tuboData),
                         sample_name = sub(".*_(.*)",
                                           "\\1",
                                           colnames(tuboData)),
                         stringsAsFactors = FALSE,
                         row.names = colnames(tuboData))
stopifnot(all(cell_names$sample_name %in% meta$sample_name))
meta_seurat <- left_join(x = cell_names,
                         y = as.data.frame(meta),
                         by = "sample_name")
# Seurat needs rownames, add them manually
# left_join preserves the order of cells:
stopifnot(all.equal(meta_seurat$cell_label, cell_names$cell_label))
# add cell labels as rownames
rownames(meta_seurat) <- meta_seurat$cell_label
library(Seurat)
message("Creating Seurat object...")
tubo_ovarian<- CreateSeuratObject(tuboData,
                        project              = "OV",
                        min.cells            = 10,     # only genes > 10 cells
                        min.genes            = 200,    # only cells with > 200 genes
                        is.expr              = 0,      # expression threshold for 'detected' gene
                        normalization.method = NULL,   # no normalization for now
                        scale.factor         = 10000,  # scale factor in normalization process
                        # (not used given norm.method == NULL)
                        do.scale             = FALSE,  # gene-based z-score
                        do.center            = FALSE,  # gene-based centering
                        names.field          = 1,
                        names.delim          = "_",
                        meta.data            = meta_seurat,
                        display.progress     = TRUE)

tubo_ovarian[["percent.mt"]] <-PercentageFeatureSet(tubo_ovarian, pattern = "^MT-")
VlnPlot(tubo_ovarian,features = c("nFeature_RNA",
                                "nCount_RNA",
                                "percent.mt"),ncol = 3)
plot1<- FeatureScatter(tubo_ovarian,feature1 ="nCount_RNA",
                       feature2 = "percent.mt")
plot2<- FeatureScatter(tubo_ovarian,feature1 ="nCount_RNA", 
                       feature2 = "nFeature_RNA")
view(tubo_ovarian@meta.data)
plot1 + plot2  
tubo_ovarian$tissue<- tubo_ovarian$sample_type
VlnPlot(tubo_ovarian, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = 'tissue')
ggplot(tubo_ovarian@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
  geom_point(size = 0.8) +
  scale_colour_viridis_c(option = 'F') + 
  lims(x = c(0, NA), y = c(0, NA)) +
  facet_wrap(~patient_id, nrow = 5, scales = 'free') +
  theme_minimal()
tubo_ovarian<-subset(tubo_ovarian,subset = nFeature_RNA >300 & nFeature_RNA <5000
                   & percent.mt < 25)
tubo_ovarian<-NormalizeData(tubo_ovarian, normalization.method = "LogNormalize",
                          scale.factor = 10000)
tubo_ovarian<-FindVariableFeatures(tubo_ovarian,selection.method = "vst",
                                 nfeatures = 5000)
#top20<-head(VariableFeatures(tubo_ovarian),20)
#plot3<- VariableFeaturePlot(tubo_ovarian)
#plot4<- LabelPoints(plot = plot3,points = top20)
#plot3 + plot4
tubo_ovarian<-ScaleData(tubo_ovarian, features = rownames(tubo_ovarian))
tubo_ovarian<-RunPCA(tubo_ovarian, features = VariableFeatures(object = tubo_ovarian))
#print(tubo_ovarian[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(tubo_ovarian, dims = 1:2, reduction = "pca")
DimPlot(tubo_ovarian, reduction = "pca")
ElbowPlot(tubo_ovarian)
tubo_ovarian<-FindNeighbors(tubo_ovarian,dims = 1:20)
tubo_ovarian<-FindClusters(tubo_ovarian,resolution= 0.35)
tubo_ovarian<-RunUMAP(tubo_ovarian, dims = 1:20)
DimPlot(tubo_ovarian,reduction = "umap", label=T)
#DimPlot(tubo_ovarian,group.by = 'RNA_snn_res.0.3')
DimPlot(tubo_ovarian, reduction = 'umap', group.by = c('sample_name', 'patient_id', 'tissue', 'sample_site'), ncol = 2)
#saveRDS(tubo_ovarian, file = "./tubo_ovarian_tutorial.rds")
library('celldex')
library('SingleR')
sc_data<- tubo_ovarian@assays$RNA@data
ref1 <- HumanPrimaryCellAtlasData()

celltype_hpca_fine<- SingleR(test = sc_data, ref = ref1, assay.type.test=1,
                             labels = ref1$label.fine)
celltype_hpcamain<-SingleR(test = sc_data, ref = ref1, assay.type.test=1,
                       labels = ref1$label.main)
ref_encode<-BlueprintEncodeData() 
celltype_encode_main<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.main)
celltype_encode_fine<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.fine)
tubo_ovarian$celltype_main<- celltype_encode_main$pruned.labels
tubo_ovarian$celltype_fine<- celltype_encode_fine$pruned.labels
p2<-DimPlot(tubo_ovarian, reduction = 'umap', group.by = 'celltype_fine', label = T, repel = T)
p3<-DimPlot(tubo_ovarian, reduction = 'umap', group.by = 'celltype_main', label = T, repel = T)
x11()
view(table(tubo_ovarian$celltype_main))
view(table(tubo_ovarian$celltype_fine))

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
## DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = c("Immune system") # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = tubo_ovarian[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(tubo_ovarian@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(tubo_ovarian@meta.data[tubo_ovarian@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(tubo_ovarian@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])
tubo_ovarian@meta.data$ScType_cell = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  tubo_ovarian@meta.data$ScType_cell[tubo_ovarian@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
p5<-DimPlot(tubo_ovarian, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'ScType_cell')     
X11()
p3 + p5
saveRDS(tubo_ovarian, file = "./tubo_ovarian_annotated.rds")

VlnPlot(tubo_ovarian, features = c("CD79A", "CD79B", "MS4A1", "CD19", 'CD3D', 'KRT17'), log = T, group.by = 'celltype_main')
FeaturePlot(tubo_ovarian, features =c("CD79A", "CD79B", "MS4A1", "CD19", "CD3D", "NKG7","PRF1", "GZMA", 
"CD14", "CD68", "CD8A", "KRT17","KRT6A") )

my_cell<-tubo_ovarian@meta.data$ScType_cell
#data<-subset(pbmc,subset = nFeature_RNA >200)
ns.data = prep_cell_mat(ge = tubo_ovarian@assays$RNA, celltypes = tubo_ovarian$seurat_clusters, logfc.threshold = 0.5)
#table(ns.data$celltypes)
view(table(ns.data$celltypes))
set.seed(0) # for reproducibility
length(unique(ns.data$celltypes))

mt = multiplet_types(celltypes = ns.data$celltypes)
am = artificial_multiplets(cell.mat = ns.data$cell.mat, 
                           celltypes = ns.data$celltypes, 
                           multiplet_classes = mt)


multiplet_types = function(celltypes, n = 2, exclude = NULL){
  celltypes = celltypes %>% as.character() %>% unique()
  
  if(!is.null(exclude)){celltypes = setdiff(celltypes, exclude)}
  d = list()
  for(i in 2:n){
    d[[i]] = combinations(length(celltypes), i, celltypes, repeats.allowed = T) %>% data.frame()
  }
  rbindlist(d, fill = T) %>% data.frame()
}


rf = multiplet_rf(am)
#view(rf$test)
#> [1]  train-mlogloss:3.545750
mroc = mroc_format(rf$test, rf$pred)
mroc_plot(mroc)
cluster_dotplot()
#p5
#view(tubo_ovarian)
#library(ggraph)
pred = xgpred(rf, ns.data$cell.mat)
result = multiplet_pop_sig(pred = pred, sample = tubo_ovarian$tissue)
view(result)
plot_interactions(result, 
                  legend.position = 'right', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('combined tissue')
#> Loading required package: ggraph
#> Loading required package: ggplot2
set.seed(0) # for reproducibility
ns.data = prep_cell_mat(ge = tubo_ovarian@assays$RNA, logfc.threshold = 0.5, celltypes =tubo_ovarian$seurat_clusters )
ns_tubOvarian = neighborseq(ns.data$cell.mat, ns.data$celltypes, do.mroc = T, sample = tubo_ovarian$tissue)
ns.data$celltypes
library(ggraph)
plot_interactions(ns_tubovarian$result,
                  legend.position = 'left', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('Combined interaction')

# plot with colored nodes
#color = colorRampPalette((brewer.pal(9,"Greens")))(11)[6:12]; 
#names(color) = c('Paneth','Stem','TA','EP','Enterocyte'); color=c(color, Goblet='gold3')
x11()
result_normal<-result %>% filter(sample== 'Normal')
result_tumor<-result %>% filter(sample== 'Tumor')
p1<-plot_interactions(result_normal,
                      legend.position = 'left', 
                      width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('Normal Tissues')

p2<-plot_interactions(result_tumor,
                      legend.position = 'right', 
                      width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('Tumor Microenvironment')

p1+p2
cell_cluster<- data.frame(tubo_ovarian$seurat_clusters, tubo_ovarian$ScType_cell) %>% unique()
tubo_ovarian<- readRDS('~/tubo_ovarian_annotated')
tubo_ovarian$celltype_main
unique(result$Cell_1)
naming<- setNames(new_name, unique(result$Cell_1))
new_name= as.character(cell_cluster$tubo_ovarian.ScType_cell)
unique(mt)
########HGSOC#####
HGSOC_Data<- read_tsv('~/R/GSE146026_Izar_HGSOC_ascites_SS2_log.tsv')
view(head(HGSOC_Data,20))
column
HGSOC_count<-HGSOC_Data[-(1:5),]
izar_meta <- as.data.frame(t(HGSOC_Data[1:5,]))
izar_meta <- izar_meta %>%
  row_to_names(row_number = 1)

izar <- HGSOC_count %>%
  distinct(Cell_ID, .keep_all= TRUE) %>%
  remove_rownames %>%
  column_to_rownames(var="Cell_ID")
izar <- as.matrix(izar)
class(izar) <- "numeric"
izar <- as(izar, "dgCMatrix")
HGSOC<- CreateSeuratObject(counts = izar, project = "IZAR", meta.data =izar_meta)

HGSOC[["percent.mt"]] <-PercentageFeatureSet(HGSOC, pattern = "^MT-")
VlnPlot(HGSOC,features = c("nFeature_RNA",
                                  "nCount_RNA",
                                  "percent.mt"),ncol = 3)
plot1<- FeatureScatter(HGSOC,feature1 ="nCount_RNA",
                       feature2 = "percent.mt")
plot2<- FeatureScatter(HGSOC,feature1 ="nCount_RNA", 
                       feature2 = "nFeature_RNA")
view(HGSOC@meta.data)
plot1 + plot2  
HGSOC$tissue<- HGSOC$sample_type
VlnPlot(HGSOC, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = 'Patient')
ggplot(HGSOC@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
  geom_point(size = 0.8) +
  scale_colour_viridis_c(option = 'F') + 
  lims(x = c(0, NA), y = c(0, NA)) +
  facet_wrap(~Patient, nrow = 5, scales = 'free') +
  theme_minimal()
HGSOC<-subset(HGSOC,subset = nFeature_RNA >1000 & nFeature_RNA <9000
                     & percent.mt < 25)
HGSOC<-NormalizeData(HGSOC, normalization.method = "LogNormalize",
                            scale.factor = 10000)
HGSOC<-FindVariableFeatures(HGSOC,selection.method = "vst",
                                   nfeatures = 5000)

HGSOC<-ScaleData(HGSOC, features = rownames(HGSOC))
HGSOC<-RunPCA(HGSOC, features = VariableFeatures(object = HGSOC))
#print(HGSOC[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(HGSOC, dims = 1:2, reduction = "pca")
DimPlot(HGSOC, reduction = "pca")
ElbowPlot(HGSOC)
HGSOC<-FindNeighbors(HGSOC,dims = 1:20)
HGSOC<-FindClusters(HGSOC,resolution= 0.35)
HGSOC<-RunUMAP(HGSOC, dims = 1:20)
DimPlot(HGSOC,reduction = "umap", label=T)
#DimPlot(HGSOC,group.by = 'RNA_snn_res.0.3')
DimPlot(HGSOC, reduction = 'umap', group.by = c('Patient', 'Time'))
sc_data<- HGSOC@assays$RNA@data
ref_encode<-BlueprintEncodeData() 

celltype_encode_main<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.main)
celltype_encode_fine<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.fine)
HGSOC$celltype_main<- celltype_encode_main$pruned.labels
HGSOC$celltype_fine<- celltype_encode_fine$pruned.labels
p2<-DimPlot(HGSOC, reduction = 'umap', group.by = 'celltype_fine', label = T, repel = T)
p3<-DimPlot(HGSOC, reduction = 'umap', group.by = 'celltype_main', label = T, repel = T)
x11()
p2
es.max = sctype_score(scRNAseqData = HGSOC[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(HGSOC@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(HGSOC@meta.data[HGSOC@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(HGSOC@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])
HGSOC@meta.data$ScType_cell = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  HGSOC@meta.data$ScType_cell[HGSOC@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
p5<-DimPlot(HGSOC, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'ScType_cell')     
saveRDS(HGSOC, file = "./HGSOC_annotated.rds")
##################################10X##############################
HGSOC_Data<- read_tsv('~/R/GSE146026_Izar_HGSOC_ascites_10x_log.tsv')
view(head(HGSOC_Data,20))
column
library(janitor)
HGSOC_count<-HGSOC_Data[-(1:7),]
izar_meta_10 <- as.data.frame(t(HGSOC_Data[1:7,]))
izar_meta_10 <- izar_meta_10 %>%
  row_to_names(row_number = 1)

izar_10 <- HGSOC_count %>%
  distinct(`10x_barcode`, .keep_all= TRUE) %>%
  remove_rownames %>%
  column_to_rownames(var="10x_barcode")
izar_10 <- as.matrix(izar_10)
class(izar_10) <- "numeric"
izar_10 <- as(izar_10, "dgCMatrix")
HGSOC_10x<- CreateSeuratObject(counts = izar_10, project = "IZAR", meta.data =izar_meta_10)

HGSOC_10x[["percent.mt"]] <-PercentageFeatureSet(HGSOC_10x, pattern = "^MT-")
VlnPlot(HGSOC_10x,features = c("nFeature_RNA",
                           "nCount_RNA",
                           "percent.mt"),ncol = 3)
plot1<- FeatureScatter(HGSOC_10x,feature1 ="nCount_RNA",
                       feature2 = "percent.mt")
plot2<- FeatureScatter(HGSOC_10x,feature1 ="nCount_RNA", 
                       feature2 = "nFeature_RNA")
view(HGSOC_10x@meta.data)
plot1 + plot2  
VlnPlot(HGSOC_10x, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = 'Patient')
ggplot(HGSOC_10x@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
  geom_point(size = 0.8) +
  scale_colour_viridis_c(option = 'F') + 
  lims(x = c(0, NA), y = c(0, NA)) +
  facet_wrap(~Patient, nrow = 5, scales = 'free') +
  theme_minimal()
HGSOC_10x<-subset(HGSOC_10x,subset = nFeature_RNA >500)
HGSOC_10x<-NormalizeData(HGSOC_10x, normalization.method = "LogNormalize",
                     scale.factor = 10000)
HGSOC_10x<-FindVariableFeatures(HGSOC_10x,selection.method = "vst",
                            nfeatures = 5000)

HGSOC_10x<-ScaleData(HGSOC_10x, features = rownames(HGSOC_10x))
HGSOC_10x<-RunPCA(HGSOC_10x, features = VariableFeatures(object = HGSOC_10x))
#print(HGSOC_10x[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(HGSOC_10x, dims = 1:2, reduction = "pca")
DimPlot(HGSOC_10x, reduction = "pca")
ElbowPlot(HGSOC_10x)
HGSOC_10x<-FindNeighbors(HGSOC_10x,dims = 1:10)
HGSOC_10x<-FindClusters(HGSOC_10x,resolution= 0.5)
HGSOC_10x<-RunUMAP(HGSOC_10x, dims = 1:10)
DimPlot(HGSOC_10x,reduction = "umap", label=T)
#DimPlot(HGSOC_10x,group.by = 'RNA_snn_res.0.3')
DimPlot(HGSOC_10x, reduction = 'umap', group.by = c('Patient', 'Time'))
sc_data<- HGSOC_10x@assays$RNA@data
ref_encode<-BlueprintEncodeData() 
celltype_encode_main<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.main)
celltype_encode_fine<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.fine)
HGSOC_10x$celltype_main<- celltype_encode_main$pruned.labels
HGSOC_10x$celltype_fine<- celltype_encode_fine$pruned.labels
p2<-DimPlot(HGSOC_10x, reduction = 'umap', group.by = 'celltype_fine', label = T, repel = T)
p3<-DimPlot(HGSOC_10x, reduction = 'umap', group.by = 'celltype_main', label = T, repel = T)
x11()

es.max = sctype_score(scRNAseqData = HGSOC_10x[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(HGSOC_10x@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(HGSOC_10x@meta.data[HGSOC_10x@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(HGSOC_10x@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])
HGSOC_10x@meta.data$ScType_cell = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  HGSOC_10x@meta.data$ScType_cell[HGSOC_10x@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
p5<-DimPlot(HGSOC_10x, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'ScType_cell')     
saveRDS(HGSOC_10x, file = "./HGSOC_10x_10X_annotated.rds")

my_cell<-HGSOC_10x@meta.data$celltype_main
#data<-subset(pbmc,subset = nFeature_RNA >200)
ns.data = prep_cell_mat(ge = HGSOC@assays$RNA,celltypes = HGSOC$seurat_clusters, logfc.threshold = 0.5)
table(ns.data$celltypes)
ns.data$markers
set.seed(0) # for reproducibility
mt = multiplet_types(ns.data$celltypes)
am = artificial_multiplets(cell.mat = ns.data$cell.mat, 
                           celltypes = ns.data$celltypes, 
                           multiplet_classes = mt)
rf = multiplet_rf(am)
#view(rf$test)
#> [1]  train-mlogloss:3.545750
mroc = mroc_format(rf$test, rf$pred)
mroc_plot(mroc)
#p5
#view(HGSOC_10x)
#library(ggraph)
pred = xgpred(rf, ns.data$cell.mat)
result = multiplet_pop_sig(pred = pred)
view(result)
plot_interactions(ns_10x$combined_result, 
                  legend.position = 'right', 
                  width_range = c(0.5,1))
ns_10x$combined_result
ns_10x = neighborseq(ns.data$cell.mat, ns.data$celltypes,iter =10)

ns.data = prep_cell_mat(ge = HGSOC@assays$RNA, logfc.threshold = 0.5, celltypes = HGSOC$seurat_clusters)
ns_HGSOC = neighborseq(ns.data$cell.mat, ns.data$celltypes, do.mroc =F)
plot_interactions(ns_HGSOC$result, 
                  legend.position = 'right', 
                  width_range = c(0.5,1))

######################################################OVC####
OVC_Data <- Read10X(data.dir = '/home/sodiq/R/2100-Ovariancancer_counts/export/OvC_counts')
saveRDS(OVC_Data, file = "./OVC_count.rds")
OVC_meta<- read_csv('~/R/2101-Ovariancancer_metadata.csv')
OVC_meta<- as.data.frame(OVC_meta)
rownames(OVC_meta)<- OVC_meta$Cell
OVC<- CreateSeuratObject(OVC_Data, object='ov',meta.data = OVC_meta, names.field = 1,
                         names.delim          = "_",
                         display.progress     = TRUE)
table(OVC$PatientNumber)
OVC[["percent.mt"]] <-PercentageFeatureSet(OVC, pattern = "^MT-")
VlnPlot(OVC,features = c("nFeature_RNA",
                               "nCount_RNA",
                               "percent.mt"),ncol = 3)
OVC<-subset(OVC,subset = nFeature_RNA >1000)
OVC<-NormalizeData(OVC, normalization.method = "LogNormalize",
                         scale.factor = 10000)
OVC<-FindVariableFeatures(OVC,selection.method = "vst",
                                nfeatures = 5000)

OVC<-ScaleData(OVC, features = rownames(OVC))
OVC<-RunPCA(OVC, features = VariableFeatures(object = OVC))
DimPlot(OVC, reduction = "pca")
ElbowPlot(OVC)
OVC<-FindNeighbors(OVC,dims = 1:20)
OVC<-FindClusters(OVC,resolution= 0.3)
OVC<-RunUMAP(OVC, dims = 1:20)
DimPlot(OVC,reduction = "umap", label=T)
sc_data<- OVC@assays$RNA@data
ref_encode<-BlueprintEncodeData() 
celltype_encode_main<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.main)
celltype_encode_fine<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.fine)
OVC$celltype_main<- celltype_encode_main$pruned.labels
OVC$celltype_fine<- celltype_encode_fine$pruned.labels
p2<-DimPlot(OVC, reduction = 'umap', group.by = 'celltype_fine', label = T, repel = T)
p3<-DimPlot(OVC, reduction = 'umap', group.by = 'celltype_main', label = T, repel = T)

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
## DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
es.max = sctype_score(scRNAseqData = OVC[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(OVC@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(OVC@meta.data[OVC@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(OVC@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])
OVC@meta.data$ScType_cell = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  OVC@meta.data$ScType_cell[OVC@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(OVC, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'ScType_cell')     
saveRDS(OVC, file = "./OVC_annotated.rds")

ns.data = prep_cell_mat(ge = OVC_annotated@assays$RNA, logfc.threshold = 0.5)
ns_10x = neighborseq(ns.data$cell.mat, ns.data$celltypes, celltypes = OVC_annotated$seurat_clusters)
plot_interactions(ns_10x$result, 
                  legend.position = 'right', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('combined tissue')
view(OVC_annotated@meta.data)
unique(as.data.frame(OVC_annotated$seurat_clusters, OVC_annotated$ScType_cell))

Regner_dirs <- list.dirs(path = '~/Downloads/GSE173682_RAW', recursive = F, full.names = F)
Regner_list<- as.list(Regner_dirs)
for (i in Regner_dirs){
  Regner_file<- list.files(paste0('~/Downloads/GSE173682_RAW/',i))
}
Regner_file
for (i in Regner_dirs){
  
}
GSM5276938<- Read10X('GSM5276938')
saveRDS(GSM5276938, 'matrix')
GSM5276939<-Read10X('GSM5276939')
saveRDS(GSM5276939, 'matrix1')
GSM5276940<- Read10X('GSM5276940')
saveRDS(GSM5276940, 'matrix2')
GSM5276941<- Read10X('GSM5276941')
saveRDS(GSM5276941, 'matrix3')
GSM5276942<-Read10X('GSM5276942')
saveRDS(GSM5276942, 'matrix4')
GSM5276943<-Read10X('GSM5276943')
saveRDS(GSM5276943, 'matrix5')
GSM5276938
GSM5276938[1:10,1:10]
colnames(GSM5276938[1:10])

CreateSeuratObject(GSM5276938)

##########omentum###################

omentum_dir<- list.files(path = '~/Downloads/GSE147082_RAW')
omentum_dir_list<- as.list(rep(NA, length(omentum_dir)))
omentum_matrix_list<- as.list(rep(NA, length(omentum_dir)))
for (x in 1:length(omentum_dir_list)){
  omentum_matrix_list[[x]]<- read_csv(omentum_dir[x])
}
omentum_combined<- Reduce(cbind, omentum_matrix_list)
omentum_obj_list<- as.list(rep(NA, length(omentum_matrix_list)))
for (x in 1:length(omentum_obj_list)){
  omentum_obj_list[[x]]<- CreateSeuratObject(omentum_matrix_list[[x]], 
                                             min.cells = 3, min.features = 300)
}
omentum_obj_1<- omentum_obj_list[[1]]
omentum_obj_2<- omentum_obj_list[[2]]@meta.data
omentum_obj_3<- omentum_obj_list[[3]]@meta.data
omentum_obj_4<- omentum_obj_list[[4]]@meta.data
omentum_obj_5<- omentum_obj_list[[5]]@meta.data
omentum_obj_6<- omentum_obj_list[[6]]@meta.data
Omentum_cancer<- merge(x=omentum_obj_1, y=omentum_obj_2)

omentum_obj_1[["percent.mt"]] <-PercentageFeatureSet(omentum_obj_1, pattern = "^MT-")
VlnPlot(omentum_obj_1,features = c("nFeature_RNA",
                         "nCount_RNA",
                         "percent.mt"),ncol = 3)
omentum_obj_1<-subset(omentum_obj_1,subset = nFeature_RNA >1500)
omentum_obj_1<-NormalizeData(omentum_obj_1, normalization.method = "LogNormalize",
                   scale.factor = 10000)
omentum_obj_1<-FindVariableFeatures(omentum_obj_1,selection.method = "vst",
                          nfeatures = 5000)

omentum_obj_1<-ScaleData(omentum_obj_1, features = rownames(omentum_obj_1))
omentum_obj_1<-RunPCA(omentum_obj_1, features = VariableFeatures(object = omentum_obj_1))
DimPlot(omentum_obj_1, reduction = "pca")
ElbowPlot(omentum_obj_1)
omentum_obj_1<-FindNeighbors(omentum_obj_1,dims = 1:20)
omentum_obj_1<-FindClusters(omentum_obj_1,resolution= 0.3)
omentum_obj_1<-RunUMAP(omentum_obj_1, dims = 1:20)
DimPlot(omentum_obj_1,reduction = "umap", label=T)
sc_data<- omentum_obj_1@assays$RNA@data
ref_encode<-BlueprintEncodeData() 
celltype_encode_main<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.main)
celltype_encode_fine<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.fine)
omentum_obj_1$celltype_main<- celltype_encode_main$pruned.labels
omentum_obj_1$celltype_fine<- celltype_encode_fine$pruned.labels
p2<-DimPlot(omentum_obj_1, reduction = 'umap', group.by = 'celltype_fine', label = T, repel = T)
p3<-DimPlot(omentum_obj_1, reduction = 'umap', group.by = 'celltype_main', label = T, repel = T)

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
## DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
es.max = sctype_score(scRNAseqData = omentum_obj_1[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(omentum_obj_1@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(omentum_obj_1@meta.data[omentum_obj_1@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(omentum_obj_1@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])
omentum_obj_1@meta.data$ScType_cell = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  omentum_obj_1@meta.data$ScType_cell[omentum_obj_1@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(omentum_obj_1, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'ScType_cell')     
saveRDS(omentum_obj_1, file = "./omentum_obj_1_annotated.rds")

pancreas.query <- pancreas.list[["indrop"]]
DefaultAssay(ref.combined) <- 'integrated'
pancreas.anchors <- FindTransferAnchors(reference = ref.combined, query = pancreas.query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = ref.combined$celltype,
                            dims = 1:30)
view(predictions)


######Marker refinement###############################################################
tubo_ovarian<- tubo_ovarian_annotated
VlnPlot(tubo_ovarian, features = c('STAR','FOXL2','DLK1', 'ARX', 'LUM'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('PDGFRA','COL1A1','COL1A2', 'BGN', 'DCN','POSTN','COL6A1'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('FABP4','ADIPOQ','ACSL1', 'PDGFRA', 'PDGFRB'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('PDGFRA','COL1A1','COL1A2', 'BGN', 'DCN','COL6A1'), log = T, group.by = 'ScType_cell')
VlnPlot(tubo_ovarian, features = c('CD3D','CD3E', 'TRAC', 'CD19', 'CD79A','IGKC'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('CD1C','CD207','CLEC9A', 'CD1A', 'HLA-DRA','CLEC10A','LILRA4'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('CD1C','ITGAX', 'CD1A', 'CD1B','CD207','ITGAM','CD209'), log = T, group.by ='seurat_clusters')
VlnPlot(tubo_ovarian, features = c('CD68','CD14','FCGR3A', 'LYZ', 'MARCO','CD86','HLA-DRB1'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('CD8A','PRF1','NCR1', 'KLRC1', 'KLRC3','CD8'), log = T, group.by = 'ScType_cell')
VlnPlot(tubo_ovarian, features = c('EPCAM','KRT18', 'KRT6A','TNNT2','BAMBI'), log = T, group.by = 'ScType_cell')
VlnPlot(tubo_ovarian, features = c('CLDN5','PECAM1','VWF'), log = T, group.by = 'ScType_cell')
VlnPlot(tubo_ovarian, features = c('IGKC','CD19', 'CD24','SDC1','JCHAIN'), log = T, group.by = 'ScType_cell')

library(scuttle)
BiocManager::install('SingleCellExperiment()')
tubo_ovarian_annotated@meta.data %>% ggplot(aes(x=celltype_main, fill= tissue))+
  geom_bar() + facet_wrap(facets = 'patient_id' ) + 
 theme(axis.text.x = element_text(angle = 90))

tubo_ovarian_annotated@meta.data %>% group_by(celltype_main) %>% tally() %>%
  filter(n >10)
plot_df <- tubo_ovarian_annotated@meta.data %>% 
  #dplyr::select(tissue, celltype_main) %>% 
  group_by(tissue, celltype_main) %>% 
  tally() %>% ungroup() %>% group_by(tissue) %>% 
  mutate(per = n/sum(n)) %>% select(-n)
ggplot2:: ggplot(plot_df, aes(x=tissue, fill= celltype_main)) + geom_bar()

my_Obj<- readRDS('Obj_Counts.rds')
2+2
head(my_Obj)
dim(my_Obj)
write.csv(my_Obj, file = "Regner.csv", append = FALSE, quote = TRUE, sep = " ")
obj<-read_csv('Regner.csv')
readRDS('Obj_Counts.rds')
MY_BC<-readRDS('BC_Counts.rds')
write.csv(as.matrix(MY_BC), file='MY_BC.csv')
dim(MY_BC)
library(tidyverse)
data_lda<- readRDS('Obj_Counts.rds')
#####################################################################
if (!require(ReactomeGSA))
BiocManager::install("ReactomeGSA", force = T)
#> Loading required package: ReactomeGSA
# install the ReactomeGSA.data package for the example data
if (!require(ReactomeGSA.data))
  BiocManager::install("ReactomeGSA.data")
library(ReactomeGSA.data)
data(jerby_b_cells)
library(ReactomeGSA)
gsva_result <- analyse_sc_clusters(jerby_b_cells, verbose = TRUE)
pathway_expression <- pathways(gsva_result)
# simplify the column names by removing the default dataset identifier
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

pathway_expression[1:3,]

########to identify most relevant pathways, assess the maximum difference in expression for every pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))
max_difference$diff <- max_difference$max - max_difference$min

# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

head(max_difference)
####plotting
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[23])
# Additional parameters are directly passed to gplots heatmap.2 function
dev.off()
plot_gsva_heatmap(gsva_result, max_pathways = 10, margins = c(6,20))
# limit to selected B cell related pathways
relevant_pathways <- c("R-HSA-983170", "R-HSA-388841", "R-HSA-2132295", "R-HSA-983705", "R-HSA-5690714")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,30), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left

#####PCA based on pathways
plot_gsva_pca(gsva_result)
#
###########scGSVA###############
library(devtools)
install_github("guokai8/scGSVA", force = T)
set.seed(123)   
library(scGSVA)   
hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
res<-scgsva(jerby_b_cells,hsko)
vlnPlot(res,features="Cell.cycle")
dotPlot(res,features="Cell.cycle")

install.packages("ggVennDiagram")
library(ggVennDiagram)
# Example data
set1 <- c(1, 2, 3, 4, 5, 6, 7, 8:20)
set2 <- c(4, 5, 6, 7, 8:25)
set3 <- c(5, 6, 8, 9:40)
set4<-c(15, 16:50)
set5<- c(14, 14:30)
# Generate Venn diagram
x = list(set1 = set1, set2 = set2, set3 = set3, myset= set4)
ggVennDiagram(x)
ggVennDiagram(x, label = 'percent')
ggVennDiagram(x, label = 'count')
ggVennDiagram(x, label = 'percent', label_alpha =0)
# Venn diagram with custom colors
ggVennDiagram(x) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
# Venn diagram with custom border
ggVennDiagram(x, color = "black", lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")


# Venn diagram without legend
ggVennDiagram(x, color = 1, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

###from ggvenn
venn_vec<- c(set1, set2, set3, set4)
ven_df<- data.frame(venn_vec)

ven_df$set1<- ven_df$venn_vec %in% set1
ven_df$set2<- ven_df$venn_vec %in% set2
ven_df$set3<- ven_df$venn_vec %in% set3
ven_df$set4<- ven_df$venn_vec %in% set4

ggplot(ven_df, aes(A=set1, B= set2, C= set3, D= set4)) + geom_venn(show_percentage = F)
BiocManager::install("DEP")
library("DEP")
BiocManager::install('ncdf4')
# Loading a package required for data handling
library("dplyr")


#########olbretch neignorseq######################
library(tidyverse)
library(Seurat)
library(Neighborseq)
library(igraph)
result_new <- readRDS('result_survival.rds')
result_PR<-result_new %>% filter(sample== 'PR')
result_CR<-result_new %>% filter(sample== 'CR')
result_PD<-result_new %>% filter(sample== 'PD')

#table(Tubo_tumor_annotated$sample_site)
plot_interactions(result_PR,
                      legend.position = 'left', 
                      width_range = c(0.5,1), min_color = 'red', max_color = 'green') + 
  ggtitle( 'Partial remission') +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

plot_interactions(result_CR,
                  legend.position = 'left', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green') + 
  ggtitle( 'Complete remission') +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

plot_interactions(result_PD,
                  legend.position = 'left', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green') + 
  ggtitle( 'Progressive disease') +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

#######################neighborseq functions############################
#' Prepare Neighbor-seq inputs: gene expression matrix, clusters, and cluster marker genes
#'
#' @param ge Gene by barcode single-cell gene expression matrix
#' @param celltypes Predetermined cell cluster annotations - optional
#' @param logfc.threshold Minimum average log fold change of a gene in a cell-type compared to all other cell-types to be tested for differential expression. Default is 1; decreaseing the threshold will find weaker signals but will increase computation time. See documentation for Seurat FindMarkers function for more details.
#' @param min.pct Minimum percent of cells from a cluster to express a gene for that gene to be tested for differential expression. Default is 0.25. See documentation for Seurat FindMarkers function for more details.
#' @param max.cells Maximum number of cells to randomly select from each cell cluster during differential gene expression testing. Default is 200.
#' @param topn Number of marker genes to select for each cell-type. The union of these - there may be overlap-  will be the genes returned in the output gene expression matrix. Top arker genes are selected based on their fold-change rank.
#' @param res Resolution parameter for Seurat clustering. Larger values will find a larger number of clusters. See Seurat FindClusters function for more details.
#' @param ... Arguments passed to other methods
#'
#' @return A list containing the cell matrix, cell-types, and marker gene output
#' @export
#'
#' @import Seurat
#' @import tidyverse
#' @import xgboost
#' @import multiROC
#' @import parallel
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import ggraph
#' @import metap
#' @import RColorBrewer
#' @import stringr
#' @import ggpubr
#' @import ggplot2
#' @importFrom tidyr "separate"
#' @importFrom dplyr "sample_frac"
#' @importFrom dplyr "group_by"
#' @importFrom dplyr "top_n"
#' @importFrom dplyr "anti_join"
#' @importFrom dplyr "bind_rows"
#' @importFrom dplyr "summarize"
#' @importFrom gtools "combinations"
#' @importFrom tibble "rownames_to_column"
#' @importFrom tibble "column_to_rownames"
#'
prep_cell_mat = function(ge, celltypes = NULL, logfc.threshold = 1, min.pct = 0.25, max.cells = 200, topn = 50, res = 0.8, ...){
  
  ge = ge %>%
    CreateSeuratObject() %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 5000)
  
  if(is.null(celltypes)){
    ge = ScaleData(ge) %>% RunPCA() %>% FindNeighbors() %>% FindClusters()
    while(nlevels(ge$seurat_clusters) > 12){
      res = 0.5*res; ge = FindClusters(ge, resolution = res)
    }
    celltypes = ge$seurat_clusters
  }
  
  Idents(ge) = celltypes
  markers = FindAllMarkers(ge, logfc.threshold = logfc.threshold,
                           only.pos = T,
                           min.pct = min.pct,
                           max.cells.per.ident = max.cells, ...)
  
  wt = ifelse('avg_log2FC' %in% colnames(markers), 'avg_log2FC', 'avg_logFC')
  m = markers %>% group_by(cluster) %>% top_n(n = topn, wt = wt)
  m = unique(m$gene)
  
  list(cell.mat = ge@assays$RNA@counts[m, ], celltypes = celltypes, markers = markers)
}

#' Enumerate all possible cell-type combinations of n cells
#'
#' @param celltypes A vector of cell-type annotations. There must be one annotation for each barcode.
#' @param n Multiplet degree, i.e. the maximum number of cells in each multiplet. Combinations will be included for integers from 2 to n.
#' @param exclude Optional cell-type(s) to exclude when enumerating multiplet-types.
#'
#' @return A dataframe of n columns with each row containing a multiplet type. Empty columsn are filled with NA.
#' @export
#'

multiplet_types = function(celltypes, n = 2, exclude = NULL){
  celltypes = celltypes %>% as.character() %>% unique()
  
  if(!is.null(exclude)){celltypes = setdiff(celltypes, exclude)}
  d = list()
  for(i in 2:n){
    d[[i]] = combinations(length(celltypes), i, celltypes, repeats.allowed = T) %>% data.frame()
  }
  rbindlist(d, fill = T) %>% data.frame()
}

#' Construct a dataset of artificial multiplets
#'
#' @param cell.mat Marker gene by barcode matrix
#' @param celltypes A vector of cell-type annotations
#' @param multiplet_classes A dataframe of multiplet types
#' @param n Number of artificial multiplets to construct for each multiplet type. Will also include n singlets from each cell-type.
#'
#' @return A dataframe where each row is a singlet or multiplet and the columns contain genes counts. The last column called 'type' contains the barcode type.
#' @export
#'
artificial_multiplets = function(cell.mat, celltypes, multiplet_classes, n = 100){
  celltypes = as.character(celltypes)
  
  i = which(celltypes %in% (c(multiplet_classes) %>% unlist() %>% as.character() %>% unique()))
  celltypes = celltypes[i]
  cell.mat = cell.mat[, i]
  
  multiplet_type = apply(multiplet_classes, 1, function(x) paste(x[!is.na(x)], collapse = '_'))
  idx = sapply(unique(celltypes) %>% sort(), function(x) which(celltypes == x), simplify = F)
  
  am = list()
  counter = 0
  for(i in 1:nrow(multiplet_classes)){
    ii = which(!is.na(multiplet_classes[i,]))
    for(j in 1:n){
      # print(paste(i,j))
      iii = sapply(ii, function(x) idx[[multiplet_classes[i,x]]] %>% sample(1))
      xx = data.frame(cell.mat[, iii] %>% rowSums() %>% t(), type = multiplet_type[i], check.names = F)
      counter = counter + 1
      am[[counter]] = xx
    }
  }
  am = rbindlist(am)
  
  if(any(table(celltypes) < n)){
    print('warning: artificial multiplets constructed by sampling with replacement'); replace = T
  } else {replace = F}
  
  xx = lapply(idx, function(x) cell.mat[, sample(x, n, replace = replace)] %>% t() %>% data.frame(check.names=F)) %>% rbindlist()
  xx = cbind(xx, type = rep(names(idx), each = n))
  
  colnames(am) = str_replace_all(colnames(am), '-', '.')
  colnames(xx) = str_replace_all(colnames(xx), '-', '.')
  am = rbind(am, xx)
  am$type = factor(am$type)
  am
}

#' Train and test a random forest classifier to predict barcode composition
#'
#' @param artificial_multiplets_mat Dataframe containing artificial multiplets. Each row is a singlet or multiplet and the columns are gene counts. The last column contains the barcode type.
#' @param f Fraction of multiplets to use for training the random forest
#' @param ... Arguments passed to other methods
#'
#' @return A list containing the trained random forest, barcode-type prediction probabilities for each barcode in the test dataset, the training dataset, and the test dataset
#' @export
#'
multiplet_rf = function(artificial_multiplets_mat, f = 0.8, ...){
  artificial_multiplets_mat$type = factor(artificial_multiplets_mat$type)
  artificial_multiplets_mat = artificial_multiplets_mat %>% rownames_to_column('n')
  train = artificial_multiplets_mat %>% group_by(type) %>% sample_frac(f, replace = F)
  test = anti_join(artificial_multiplets_mat, train, by = 'n')
  artificial_multiplets_mat = column_to_rownames(artificial_multiplets_mat, 'n')
  train = column_to_rownames(train, 'n')
  test = column_to_rownames(test, 'n')
  
  train.type = train[, 'type']
  train.typex = train[, 'type'] %>% as.numeric(); train.typex = train.typex - 1
  train = train[,-ncol(train)] %>% as.matrix()
  
  ncore = detectCores()
  
  rf <- xgboost(data = train, label = train.typex, objective = "multi:softprob",
                eval_metric = "mlogloss",
                num_class = nlevels(artificial_multiplets_mat$type),
                nthread = ncore, nround = 1, max_depth = 20,
                num_parallel_tree = 200, subsample = 0.632,
                colsample_bytree  = sqrt(ncol(train))/ncol(train),
                colsample_bynode  = sqrt(ncol(train))/ncol(train))
  
  test.type = test[, 'type']
  test = test[, -ncol(test)] %>% as.matrix()
  
  pred = predict(rf, test)
  pred <- matrix(pred,
                 nrow = nlevels(artificial_multiplets_mat$type),
                 ncol = length(pred)/nlevels(artificial_multiplets_mat$type)) %>%
    t() %>%
    data.frame(check.names = F)
  
  colnames(pred) = levels(artificial_multiplets_mat$type)
  
  list(rf = rf, pred = pred, train = data.frame(train, type = train.type, check.names = F), test = data.frame(test, type = test.type, check.names = F))
}

#' Calculate multivariate receiver operator sensitivities and specificities
#'
#' @param test Test dataset used to evalute the multiplet random forest. The last column must be labeled 'type' and contain the barcode type.
#' @param pred Prediction probabilities for each barcode.
#'
#' @return A dataframe containing the sensitivity, specificity, and AUC for each group. Also contains the micro and macro avareages of all groups. See multiROC package documentaiton for more details.
#' @export
#'
mroc_format = function(test, pred){
  test$type = factor(test$type)
  classes = levels(test$type)
  mat = matrix(0, nrow = nrow(pred), ncol = length(classes))
  colnames(mat) = classes
  for(i in 1:nrow(mat)){
    mat[i, test$type[i]] = 1
  }
  colnames(mat) = colnames(mat) %>% str_replace_all(' ', '_') %>% paste0('_true')
  pred.mat = pred
  colnames(pred.mat) = colnames(pred.mat) %>% str_replace_all(' ', '_') %>% paste0('_pred_RF')
  
  df = data.frame(mat, pred.mat, check.names = F)
  mroc = suppressWarnings(multi_roc(df, force_diag = T))
  
  # format data.frame for plotting
  n_method <- length(unique(mroc$Methods))
  n_group <- length(unique(mroc$Groups))
  mroc_df <- data.frame(Specificity= numeric(0), Sensitivity= numeric(0),
                        Group = character(0), AUC = numeric(0),
                        Method = character(0))
  for (i in 1:n_method) {
    for (j in 1:n_group) {
      temp_data_1 <- data.frame(Specificity=mroc$Specificity[[i]][j],
                                Sensitivity=mroc$Sensitivity[[i]][j],
                                Group=unique(mroc$Groups)[j],
                                AUC=mroc$AUC[[i]][j],
                                Method = unique(mroc$Methods)[i])
      colnames(temp_data_1) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
      mroc_df <- rbind(mroc_df, temp_data_1)
      
    }
    temp_data_2 <- data.frame(Specificity=mroc$Specificity[[i]][n_group+1],
                              Sensitivity=mroc$Sensitivity[[i]][n_group+1],
                              Group= "Macro",
                              AUC=mroc$AUC[[i]][n_group+1],
                              Method = unique(mroc$Methods)[i])
    temp_data_3 <- data.frame(Specificity=mroc$Specificity[[i]][n_group+2],
                              Sensitivity=mroc$Sensitivity[[i]][n_group+2],
                              Group= "Micro",
                              AUC=mroc$AUC[[i]][n_group+2],
                              Method = unique(mroc$Methods)[i])
    colnames(temp_data_2) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
    colnames(temp_data_3) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
    mroc_df <- rbind(mroc_df, temp_data_2)
    mroc_df <- rbind(mroc_df, temp_data_3)
  }
  
  mroc_df
}

#' Plot multivariate receiver operator curves
#'
#' @param mroc_df multiROC result dataframe obtained from mroc_format function
#' @param c_size linewidth for individual barcode type curves
#' @param a_size linewidth for average ROC curve
#'
#' @return A ggplot object graphing multiROC curves
#' @export
#'
mroc_plot = function(mroc_df, c_size = 0.3, a_size = 1){
  macro.auc = mroc_df$AUC[mroc_df$Group == 'Macro'][1] %>% round(3)
  ggplot(mroc_df, aes(x = 1-Specificity, y=Sensitivity)) +
    geom_line(data = subset(mroc_df, Group != 'Macro'), aes(color = Group), size = c_size, show.legend = F) +
    scale_color_manual(values = colorRampPalette(brewer.pal(9, 'Greys'))(length(unique(mroc_df$Group)))) +
    geom_line(data = subset(mroc_df, Group == 'Macro'), color = 'Firebrick', size = a_size, show.legend = T) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    theme_classic() +
    scale_x_continuous(limits = c(-0.05, 1.05), expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0,0)) +
    ggtitle(paste('Artificial droplet classification\nAvg AUC =', macro.auc)) +
    theme(axis.text = element_text(colour = 'black', size = 10),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          title = element_text(size=9),
          legend.position = 'none')
}


#' Format dataframe to create a clustered dot-plot
#'
#' @param df A dataframe in long format containing data for the dotplot
#' @param rows Name of dataframe column containing row labels for dotplot
#' @param columns Name of dataframe column containing column labels for dotplot
#' @param values Name of dataframe column containing values for dotplot
#'
#' @return A long format dataframe with factored and leveled row and column labels
#' @export
#'
cluster_dotplot = function(df, rows, columns, values){
  xx = df %>% pivot_wider(id_cols = rows,
                          names_from = columns,
                          values_from = values) %>%
    column_to_rownames(rows)
  xx[is.na(xx)] = 0
  row.order = hclust(dist(xx))$order
  column.order = hclust(dist(t(xx)))$order
  df = df %>% mutate(!!columns := factor(!!as.name(columns), levels = colnames(xx)[column.order]))
  df = df %>% mutate(!!rows := factor(!!as.name(rows), levels = rownames(xx)[row.order]))
  return(df)
}

#' Barcode-type prediction using a trained random forest
#'
#' @param rf Random forest trained on artificial multiplets. Output of multiplet_rf
#' @param cell.mat Gene by cell matrix
#'
#' @return Dataframe of barcode-type prediction probabilities
#' @export
#'
xgpred = function(rf, cell.mat){
  # format gene expresion matrix and gene names for the random forest
  pred.mat = t(cell.mat)
  colnames(pred.mat) = str_replace_all(colnames(pred.mat), '-', '.')
  # colnames(pred.mat)[str_which(colnames(pred.mat), '^[0-9]')] =
  #   paste0('X', colnames(pred.mat)[str_which(colnames(pred.mat), '^[0-9]')])
  
  pred = predict(rf$rf, pred.mat, type = 'prob')
  pred <- matrix(pred,
                 nrow = nlevels(rf$train$type),
                 ncol = length(pred)/nlevels(rf$train$type)) %>%
    t() %>%
    data.frame()
  
  colnames(pred) = levels(rf$train$type)
  pred
}

#' Compute multiplet type enrichment scores and statistical significance
#'
#' @param pred Multiplet type prediction probabilities for each barcode
#' @param sample Optional sample label annotation for each barcode. Statistical significance is assessed relative to cell-type proportions from an individual sample only.
#' @param n Number of simulations to run. In ach simulation multiplets are randomly created from the underlying cell-type distributins.
#' @param homotypic True or False; whether to include homotypic multiplets when calculating enrichment and statistical significance
#'
#' @return A dataframe containing the multiplet significance results. Each row contains data for a multiplet type in a sample. The sample, multiplet type, enrichment score, p-value, and adjusted p-value are reported
#' @export
#'
multiplet_pop_sig = function(pred, sample = NULL, n = 100, homotypic = T){
  
  if(is.null(sample)){sample = rep('Sample.1', nrow(pred))}
  
  out = list()
  for(i in unique(sample)){
    
    pred.sample = pred[sample == i, ] %>% apply(1, function(x) names(x)[which.max(x)])
    celltypes = pred.sample %>% strsplit('_') %>% unlist()
    pred.sample = pred.sample[str_which(pred.sample, '_')]
    
    if(length(pred.sample) == 0){next}
    
    df = table(pred.sample) %>% data.frame() %>% separate('pred.sample', into=c('Cell_1','Cell_2'))
    colnames(df)[3] = 'Counts'
    
    if(homotypic == F){
      j = apply(df, 1, function(x) length(unique(x[1:(length(x)-1)])))
      df = df[j>1, ]
      rownames(df) = NULL
      pred.sample = pred.sample[j>1]
    }
    
    y = sapply(1:n, function(x){
      sample(celltypes, max(4, 2*length(pred.sample)), replace = F) %>%
        matrix(ncol = 2) %>%
        apply(1, sort) %>%
        t() %>%
        apply(1, paste0, collapse = '_') %>%
        table() %>%
        as.matrix() %>%
        t() %>%
        data.frame(check.names = F)
    }) %>% bind_rows(); y[is.na(y)] = 0
    
    
    types = paste0(df$Cell_1, '_', df$Cell_2)
    total_edges = c()
    for(j in unique(celltypes)){
      total_edges[j] = sum(df$Counts[str_which(types, j)])
    }
    
    e2 = c()
    for(j in 1:nrow(df)){
      e2[j] = df$Counts[j]/(total_edges[df$Cell_1[j]] * total_edges[df$Cell_2[j]])
    }
    
    names(e2) = types; e2=e2*10^4
    
    df$EnrichmentScore = e2
    
    # enrichment score for randomized data
    cn = colnames(y) %>% data.frame() %>% separate(col='.', into = c('cell1', 'cell2'))
    ct = unique(c(cn$cell1, cn$cell2))
    yedges = matrix(0, nrow = nrow(y), ncol = length(ct))
    colnames(yedges) = ct
    for(j in 1:nrow(y)){
      for(k in ct) {
        yedges[j,k] = y[j, str_which(colnames(y), k)] %>% sum()
      }
    }
    
    for(j in 1:nrow(y)){
      for(k in 1:ncol(y)){
        y[j,k] = y[j,k]/(yedges[j, cn$cell1[k]] * yedges[j, cn$cell2[k]])
      }
    }
    y = y*10^4
    
    
    sdiff = setdiff(names(e2), colnames(y))
    if(length(sdiff) > 0){y = cbind(y, matrix(0, nrow = nrow(y), ncol = length(sdiff),
                                              dimnames = list(rownames=NULL, colnames=sdiff)))}
    
    pval = sapply(names(e2), function(x) wilcox.test(y[, x], mu = e2[x], alternative = 'less')$p.value)
    padj = p.adjust(pval)
    
    
    df = data.frame(sample = i, df, pval = pval, padj); rownames(df) = NULL
    out[[i]] = df
  }
  out = rbindlist(out) %>% data.frame()
  
  out
}

#' Function to run the entire Neighbor-seq pipeline
#'
#' @param cell.mat Gene by cell matrix. Performance is improved if only marker genes are kept.
#' @param celltypes Vector containing cell-type annotations for each barcode
#' @param sample Optional vector containing sample annotations for each barcode
#' @param iter Optional number of iterations to train the random forest
#' @param exclude Optional cell-type to exclude from the artificial multiplet dataset. Can be uesful when known multiplets are present.
#' @param multiplet.degree Maximum number of cells in a multiplet
#' @param n.am Number of artificial multiplets to create for each multiplet type
#' @param f Fraction of artificial multiplets data to use for training the random forest
#' @param nsim Number of simulations to run to calculate multiplet enrichment significance
#' @param do.mroc True or False; whether to calculate multivariate receiver operator curves
#' @param homotypic True or False; whether to include homotypic multiplets when assessing enrichment and significance
#'
#' @return A list containing: 1. A combined result dataframe with multiplet enrichment scores and Fisher combined p-values if iter>1, 2. a list or dataframe with multiplet enrichment scores and p-values from each iteration, 3. a list or dataframe containing barcode type prediction probabilities for all barcodes in the dataset, 4. a list or dataframe containing the multiROC analysis
#' @export
#'
neighborseq = function(cell.mat, celltypes, sample = NULL, iter = 1, exclude = NULL, multiplet.degree = 2, n.am = 100, f = 0.8, nsim = 100, do.mroc = T, homotypic = T){
  start = Sys.time()
  mt = multiplet_types(celltypes, n = multiplet.degree, exclude = exclude)
  
  print('Creating artificial multiplets')
  am = artificial_multiplets(cell.mat, celltypes, mt, n = n.am)
  
  rf = list()
  mroc = list()
  pred = list()
  result = list()
  for(i in 1:iter){
    print(paste('Training random forest', i, '/', iter))
    rf[[i]] = multiplet_rf(am, f = f)
    
    if(do.mroc == T){
      print('Computing multi-ROC')
      mroc[[i]] = mroc_format(rf[[i]]$test, rf[[i]]$pred)
    }
    else{mroc = NULL}
    
    pred[[i]] = xgpred(rf[[i]], cell.mat)
    
    print('Computing population significance')
    result[[i]] = multiplet_pop_sig(pred[[i]], sample = sample, n = nsim, homotypic = homotypic)
    
  }
  
  if(iter > 1){
    combined_result = rbindlist(result) %>%
      group_by(sample, Cell_1, Cell_2) %>%
      summarize(MeanCounts = mean(Counts),
                CombinedP = suppressWarnings(sumlog(padj)$p))
  } else {
    result = result[[1]]
    pred = pred[[1]]
    rf = rf[[1]]
    combined_result = NULL
    if(!is.null(mroc)){mroc = mroc[[1]]}
  }
  
  print(round(Sys.time() - start), 1)
  
  list(combined_result = combined_result, result = result, pred = pred, rf = rf, mroc = mroc)
}

#####plot function########################################################
#' Plot the cell-cell interaction network
#'
#' @param df Dataframe from containing the multiplet enrichment and significance results
#' @param minCounts_sample Minimum number of multiplet type counts in a sample to include the plot
#' @param minCounts_total Minimum number of multiplet type counts in all samples to include the plot
#' @param maxP Maximum multiplet type p-value to include in the plot
#' @param padding Line padding around cell-type labels
#' @param combined True or False; whether the dataframe contains an ensemble Neighbor-seq result or a single run
#' @param layout Network graph layout; see ggraph documentation for more information
#' @param color A named vector of color lables for each cell-type
#' @param color_breaks Edge color breaks. See ggplot documentation for more details
#' @param width_breaks Edge width breaks. See ggplot documentation for more details
#' @param width_range Edge width range. See ggplot documentation for more details
#' @param min_color Edge color for smallest edge / multiplet type with fewest counts
#' @param max_color Edge color for largest edge / multiplet type with highest counts
#' @param legend.position Legend position; either 'top', 'bottom', 'left', 'right', 'none'
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object containing the cell-cell interaction network
#' @export
#'
plot_interactions = function(df, minCounts_sample = 10, minCounts_total = 10, maxP = 0.05, padding = 1, combined = F,
                             layout = 'stress', color = NULL, color_breaks = NULL, width_breaks = NULL,
                             width_range = c(0.1,1), min_color = 'grey50', max_color = 'black',
                             legend.position = 'none', ...){
  
  require(ggraph)
  
  df$type = paste0(df$Cell_1, '_', df$Cell_2)
  
  n = padding
  if(combined == T){
    df = subset(df, MeanCounts > minCounts_sample & CombinedP < maxP)
    df = df %>%
      group_by(type, Cell_1, Cell_2) %>%
      summarize(MeanCounts = sum(MeanCounts, na.rm=T),
                CombinedP = min(CombinedP, na.rm=T),
                .groups = 'keep')
    edge_width = -1*log10(df$CombinedP)
    edge_color = df$MeanCounts
  } else {
    df = subset(df, Counts > minCounts_sample & padj < maxP)
    df = df %>%
      group_by(type, Cell_1, Cell_2) %>%
      summarize(Counts = sum(Counts, na.rm=T),
                padj = min(padj),
                .groups = 'keep')
    edge_width = -1*log10(df$padj)
    edge_color = df$Counts
  }
  edge_width = edge_width[df$Cell_1 != df$Cell_2]
  edge_color = edge_color[df$Cell_1 != df$Cell_2]
  df = df[df$Cell_1 != df$Cell_2, ]
  
  graph = graph_from_data_frame(df[,2:4])
  
  if(!is.null(color)){
    color = color[V(graph) %>% as.factor() %>% names()]
  } else {
    color = rep('black', length(V(graph)))
  }
  
  if(is.null(color_breaks)){color_breaks = waiver()}
  if(is.null(width_breaks)){width_breaks = waiver()}
  
  ggraph(graph, layout = layout) +
    geom_edge_link(aes(start_cap = label_rect(node1.name, padding = ggplot2::margin(n,n,n,n,'mm')),
                       end_cap = label_rect(node2.name, padding = ggplot2::margin(n,n,n,n,'mm')),
                       width = edge_width,
                       color = edge_color)) +
    scale_edge_width(range = width_range, name = '-log10(P)', breaks = width_breaks) +
    scale_edge_color_gradient(low = min_color, high = max_color, name = 'Counts', breaks = color_breaks) +
    geom_node_text(aes(label = name), color = color, size = 7, repel = T, point.padding = NA, box.padding = 0, force = 0.001)+
    theme(panel.background = element_rect(fill = 'white'),
          legend.position = legend.position,
          legend.text = element_text(size = 14),
          legend.key.size = unit(1.5, "cm"),
          legend.title = element_text(size=15),
          legend.key=element_blank())
}



