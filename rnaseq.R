### PART 1- Pre-processing and Integration----
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(usethis) 
usethis::edit_r_environ()
library(cowplot)
library(RColorBrewer)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
install.packages("metap")
library(metap)
library(data.table)
library(magrittr)
library(Seurat)


###Creating seurat object for the first round data
Invivo.HDM_1 <- Read10X("round1")

Invivo <- CreateSeuratObject(counts = Invivo.HDM_1$`Gene Expression`, project = "Round1", min.features = -1)
Invivo@assays$RNA@counts[1:10,1:10]

##Adding HTO as an independent assay
Invivo [["HTO"]] <- CreateAssayObject(counts = Invivo.HDM_1$`Antibody Capture`, project = "Round1")

##Normalise data HTO using centered log-ratio (CLR)
Invivo <- NormalizeData(Invivo, assay = "HTO", normalization.method = "CLR")
Invivo = HTODemux(Invivo, assay = "HTO", positive.quantile = 0.99)

View(Invivo)

##Visualise demultiplexing result
VlnPlot(Invivo, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Global classification results
table(Invivo$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(Invivo) <- "HTO_maxID"
RidgePlot(Invivo, assay = "HTO", features = rownames(Invivo[["HTO"]])[1:10], ncol = 5)

##Visualize HTO signals in a heatmap
HTOHeatmap(Invivo, assay = "HTO")


###Compare number of UMIs for singlets, doublets and negative cells

Idents(Invivo) <- "HTO_classification.global"
Idents(Invivo)
VlnPlot(Invivo, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
table(Invivo$HTO_classification.global)

###Generate a two dimensional tSNE embedding for HTOs

# First, we will remove negative cells from the object
Invivo.subset <- subset(Invivo, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(Invivo.subset) <- "HTO"
Invivo.subset <- ScaleData(Invivo.subset, features = rownames(Invivo.subset),
                           verbose = FALSE)
Invivo.subset <- RunPCA(Invivo.subset, features = rownames(Invivo.subset), approx = FALSE)
Invivo.subset <- RunUMAP(Invivo.subset, dims = 1:10, perplexity = 100)
DimPlot(Invivo.subset, group.by = "HTO_classification.global")
DimPlot(Invivo.subset, group.by = "hash.ID")



####Changing to RNA
DefaultAssay(Invivo) <- "RNA"

# Subset singlet
Idents(Invivo) <- "HTO_classification.global"
Invivo.singlet <- subset(Invivo, idents = "Singlet")
##perfom standard work flow (Using the new object)
Invivo.singlet[["percent.mt"]] <- PercentageFeatureSet(Invivo.singlet, pattern = "^mt-")
VlnPlot(Invivo.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Invivo.singlet <- subset(Invivo.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
Invivo.singlet <- NormalizeData(Invivo.singlet)
Invivo.singlet <- FindVariableFeatures(Invivo.singlet, selection.method = "vst", nfeatures = 2000)
Invivo.singlet <- ScaleData(Invivo.singlet)
Invivo.singlet <- RunPCA(Invivo.singlet, features = VariableFeatures(object = Invivo.singlet))
ElbowPlot(Invivo.singlet)
Invivo.singlet <- RunUMAP(Invivo.singlet, dims = 1:20)
Invivo.singlet <- FindNeighbors(Invivo.singlet, dims = 1:20)
Invivo.singlet <- FindClusters(Invivo.singlet, resolution = 1)
DimPlot(Invivo.singlet, label = TRUE)


#####Removing Doublets for Experiment 1
##pK Identification (no ground-truth)
sweep.res = paramSweep_v3(Invivo.singlet, PCs = 1:20, sct = FALSE)
sweep.stats = summarizeSweep(sweep.res, GT = FALSE)
bcmvn = find.pK(sweep.stats)

ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + 
  geom_point() + 
  geom_line()
####Select the pK that corresponds to the max bcmvn to optimize the doubet detection
pK = bcmvn %>%
  filter(BCmetric ==max(BCmetric)) %>%
  select(pK)
pK = as.numeric(as.character(pK[[1]]))

##Homotypic Doublet Proportion Estimate
annotation = Invivo.singlet@meta.data$seurat_clusters
homotypic.prop = modelHomotypic(annotation)
nExp_poi = round(0.075*nrow(Invivo.singlet@meta.data)) ##Assuming 7.5% doublets formed
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))

###Run doublet finder
Invivo.singlet = doubletFinder_v3(Invivo.singlet,
                                  PCs = 1:20,
                                  pN = 0.25,
                                  pK = pK,
                                  nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = FALSE)

names(Invivo.singlet@meta.data)

###Visualize doublets
DimPlot(Invivo.singlet, reduction = 'umap', group.by = "DF.classifications_0.25_0.26_380")

###To get the number of singlets and doublets
table(Invivo.singlet@meta.data$DF.classifications_0.25_0.26_380)

####Subsetting the singlets for downstream analyses
Invivo.singlets_new <- subset(Invivo.singlet, subset = DF.classifications_0.25_0.26_380 == "Singlet")

Invivo.singlets_new <- NormalizeData(object = Invivo.singlets_new)
Invivo.singlets_new <- FindVariableFeatures(object = Invivo.singlets_new, selection.method = "vst", nfeatures = 2000)
Invivo.singlets_new <- ScaleData(object = Invivo.singlets_new)
Invivo.singlets_new <- RunPCA(object = Invivo.singlets_new)
ElbowPlot(Invivo.singlets_new)
Invivo.singlets_new <- RunUMAP(Invivo.singlets_new, dims = 1:20)
Invivo.singlets_new <- FindNeighbors(Invivo.singlets_new, dims = 1:20)
Invivo.singlets_new <- FindClusters(Invivo.singlets_new, resolution = 1)
DimPlot(Invivo.singlets_new, label = TRUE)



##Creating seurat object for the second round data
Invivo.HDM_2 <- Read10X("round2")

Invivo2 <- CreateSeuratObject(counts = Invivo.HDM_2$`Gene Expression`, project = "Round2", min.features = -1)
Invivo2@assays$RNA@counts[1:10,1:10]
##Adding HTO as an independent assay
Invivo2 [["HTO"]] <- CreateAssayObject(counts = Invivo.HDM_2$`Antibody Capture`, project = "Round2")

##Normalise data HTO using centered log-ratio (CLR)
Invivo2 <- NormalizeData(Invivo2, assay = "HTO", normalization.method = "CLR")
Invivo2 = HTODemux(Invivo2, assay = "HTO", positive.quantile = 0.99)

View(Invivo2)

##Visualise demultiplexing result
VlnPlot(Invivo2, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Global classification results
table(Invivo2$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(Invivo2) <- "HTO_maxID"
RidgePlot(Invivo2, assay = "HTO", features = rownames(Invivo2[["HTO"]])[1:10], ncol = 5)

##Visualize HTO signals in a heatmap
HTOHeatmap(Invivo2, assay = "HTO")



###Compare number of UMIs for singlets, doublets and negative cells

Idents(Invivo2) <- "HTO_classification.global"
Idents(Invivo2)
VlnPlot(Invivo, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
table(Invivo2$HTO_classification.global)

###Generate a two dimensional tSNE embedding for HTOs

# First, we will remove negative cells from the object
Invivo.subset2 <- subset(Invivo2, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(Invivo.subset2) <- "HTO"
Invivo.subset2 <- ScaleData(Invivo.subset2, features = rownames(Invivo.subset2),
                            verbose = FALSE)
Invivo.subset2 <- RunPCA(Invivo.subset2, features = rownames(Invivo.subset2), approx = FALSE)
Invivo.subset2 <- RunUMAP(Invivo.subset2, dims = 1:10, perplexity = 100)
DimPlot(Invivo.subset2, group.by = "HTO_classification.global")
DimPlot(Invivo.subset2, group.by = "hash.ID")

DefaultAssay(Invivo2)

####Changing to RNA
DefaultAssay(Invivo2) <- "RNA"

# Subset singlet
Idents(Invivo2) <- "HTO_classification.global"
Invivo.singlet2 <- subset(Invivo2, idents = "Singlet")
##perfom standard work flow (Using the new object)
Invivo.singlet2[["percent.mt"]] <- PercentageFeatureSet(Invivo.singlet2, pattern = "^mt-")
VlnPlot(Invivo.singlet2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

Invivo.singlet2 <- subset(Invivo.singlet2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
Invivo.singlet2 <- NormalizeData(Invivo.singlet2)
Invivo.singlet2 <- FindVariableFeatures(Invivo.singlet2, selection.method = "vst", nfeatures = 2000)
Invivo.singlet2 <- ScaleData(Invivo.singlet2)
Invivo.singlet2 <- RunPCA(Invivo.singlet2, features = VariableFeatures(object = Invivo.singlet2))
ElbowPlot(Invivo.singlet2)
Invivo.singlet2 <- RunUMAP(Invivo.singlet2, dims = 1:20)
Invivo.singlet2 <- FindNeighbors(Invivo.singlet2, dims = 1:20)
Invivo.singlet2 <- FindClusters(Invivo.singlet2, resolution = 1)
DimPlot(Invivo.singlet2, label = TRUE)


#####Removing Doublets for Experiment 2
##pK Identification (no ground-truth)
sweep.res = paramSweep_v3(Invivo.singlet2, PCs = 1:20, sct = FALSE)
sweep.stats = summarizeSweep(sweep.res, GT = FALSE)
bcmvn = find.pK(sweep.stats)

ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + 
  geom_point() + 
  geom_line()
####Select the pK that corresponds to the max bcmvn to optimize the doubet detection
pK = bcmvn %>%
  filter(BCmetric ==max(BCmetric)) %>%
  select(pK)
pK = as.numeric(as.character(pK[[1]]))

##Homotypic Doublet Proportion Estimate
annotation = Invivo.singlet2@meta.data$seurat_clusters
homotypic.prop = modelHomotypic(annotation)
nExp_poi = round(0.075*nrow(Invivo.singlet2@meta.data)) ##Assuming 7.5% doublets formed
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))

###Run doublet finder
Invivo.singlet2 = doubletFinder_v3(Invivo.singlet2,
                                   PCs = 1:20,
                                   pN = 0.25,
                                   pK = pK,
                                   nExp = nExp_poi.adj,
                                   reuse.pANN = FALSE, sct = FALSE)

names(Invivo.singlet2@meta.data)

###Visualize doublets
DimPlot(Invivo.singlet2, reduction = 'umap', group.by = "DF.classifications_0.25_0.28_319")

###To get the number of singlets and doublets
table(Invivo.singlet2@meta.data$DF.classifications_0.25_0.28_319)

####Subsetting the singlets for downstream analyses
Invivo.singlets2_new <- subset(Invivo.singlet2, subset = DF.classifications_0.25_0.28_319 == "Singlet")

Invivo.singlets2_new <- NormalizeData(object = Invivo.singlets2_new)
Invivo.singlets2_new <- FindVariableFeatures(object = Invivo.singlets2_new, selection.method = "vst", nfeatures = 2000)
Invivo.singlets2_new <- ScaleData(object = Invivo.singlets2_new)
Invivo.singlets2_new <- RunPCA(object = Invivo.singlets2_new)
ElbowPlot(Invivo.singlets2_new)
Invivo.singlets2_new <- RunUMAP(Invivo.singlets2_new, dims = 1:20)
Invivo.singlets2_new <- FindNeighbors(Invivo.singlets2_new, dims = 1:20)
Invivo.singlets2_new <- FindClusters(Invivo.singlets2_new, resolution = 1)
DimPlot(Invivo.singlets2_new, label = TRUE)





####Merging the seurat objects

merged_seurat = merge(x = Invivo.singlets_new, y = Invivo.singlets2_new, add.cell.ids = c("Round1", "Round2"))
P1 = DimPlot(merged_seurat, reduction = "umap", label = TRUE)


######Perform downstream processing for the merged data
merged_seurat <- NormalizeData(object = merged_seurat)
merged_seurat <- FindVariableFeatures(object = merged_seurat, selection.method = "vst", nfeatures = 2000)
merged_seurat <- ScaleData(object = merged_seurat)
merged_seurat <- RunPCA(object = merged_seurat)
ElbowPlot(merged_seurat)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20)
merged_seurat <- FindClusters(merged_seurat, resolution = 1)

###Creating sample column
merged_seurat$sample = rownames(merged_seurat@meta.data)
View(merged_seurat@meta.data)

##Split samples
merged_seurat@meta.data = separate(merged_seurat@meta.data, col = 'sample', into = c('Round', 'Barcode'), sep = '_')
View(merged_seurat@meta.data)

HDM.list = SplitObject(merged_seurat, split.by = "Round")

Before_Integration = DimPlot(merged_seurat, label = FALSE, group.by = 'Round')
Before_Integration

###Normalize and identify variable features for each dataset independently
for(i in 1:length(HDM.list)){
  HDM.list[[i]] = NormalizeData(object = HDM.list[[i]])
  HDM.list[[i]] = FindVariableFeatures(object = HDM.list[[i]])
}

###Performing integration
###Find integration anchors (CCA), i.e. selecting features that are repeatedly variable across datasets for integration

features = SelectIntegrationFeatures(object.list = HDM.list)

###find integration features

anchors = FindIntegrationAnchors(object.list = HDM.list,
                                 anchor.features = features)

combined.dataset = IntegrateData(anchorset = anchors)

###Scaling the new data
combined.dataset@assays$integrated
DefaultAssay(combined.dataset) = "integrated"
combined.dataset = ScaleData(object = combined.dataset)
combined.dataset = RunPCA(object = combined.dataset)
ElbowPlot(combined.dataset)
combined.dataset <- FindNeighbors(combined.dataset, dims = 1:15)
combined.dataset <- FindClusters(combined.dataset, resolution = 0.5)
combined.dataset = RunUMAP(object = combined.dataset, dims = 1:15)

DimPlot(combined.dataset, label = TRUE)
DimPlot(combined.dataset, label = TRUE) + NoLegend()

###Comparison before and after data integration

Before_Integration = DimPlot(merged_seurat, label = FALSE, group.by = 'Round')
Before_Integration

After_Integration = DimPlot(combined.dataset, label = FALSE, group.by = 'Round')
After_Integration

Before_Integration + After_Integration

####Switching from the default assay (Integrated) to RNA
DefaultAssay(combined.dataset)
DefaultAssay(combined.dataset) = "RNA"
combined.dataset = ScaleData(object = combined.dataset)


####Saving RDS file
saveRDS(combined.dataset, file = "~/Desktop/DATA/scRNA seq/scRNAseq_Integration/scRNAseq_Data_NEW.rds")

##Save SessionInfo
sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
?save.image





###PART 2- Downstream analysis----

library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(usethis) 
usethis::edit_r_environ()
library(cowplot)
library(RColorBrewer)
library(metap)
library(data.table)
library(magrittr)
library(Seurat)
library("gridExtra")
library("scCustomize")
library(ggpubr)

combined.dataset = readRDS("~/Desktop/DATA/scRNA seq/scRNAseq_Integration/scRNAseq_Data_NEW.rds")
View(combined.dataset@meta.data)

View(combined.dataset@meta.data)
table(combined.dataset$Round)
table(combined.dataset$HTO_classification.global)

DefaultAssay(combined.dataset)
DefaultAssay(combined.dataset) = "integrated"


####Reoder levels
levels(combined.dataset)
my_levels = c( "0" , "1","2" , "3" , "4" , "5" , "6" , "7" , "8" , "9", "10", "11", "12", "13")
combined.dataset@meta.data$seurat_clusters = factor(x=combined.dataset@meta.data$seurat_clusters, levels = my_levels)
##Umap
DimPlot(combined.dataset, reduction = "umap", label = TRUE, group.by = 'seurat_clusters')
DimPlot(combined.dataset, reduction = "umap", label = FALSE, group.by = 'HTO_classification')

####Getting the cells for each cluster and Sample
Cell.nums = table(combined.dataset@meta.data$HTO_classification,
                  combined.dataset@meta.data$seurat_clusters)
Cell.nums

write.csv(x=Cell.nums, file = "Number of cells per all cluster.csv", quote = FALSE)


DefaultAssay(combined.dataset) = "RNA"

###Dotplots for for lineage specific genes
##clusters
##Receptors
genes_to_plot1 = c("Csf1r","Cd14",'Pecam1',"Icam1","Vcam1",'Tcf7','Pdgfra',"Trbc2","Cd4","Cd44", "Sell","Pdcd1","Ctla4","Icos","Havcr2","Itgae","Cd69","Il2ra","Il2rg", "Il7r","Il15ra","Ccr2","Ccr3","Ccr4","Ccr6","Ccr7", "Ccr8","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Klrg1", "Il10ra","Il1rl1","Il17rb","Crlf2","Cd200r1","Lrrc32","Entpd1",
                   "Nt5e", "Lag3","Nrp1","Tnfrsf4","Tnfrsf9","Tnfrsf18", "Vdr", "Il21r", "Ptgdr2", "Klrb1", "Klrd1", "Fcgr3", "Ncam1", "Slamf1",'Slamf6',"Plin2", "Ffar3", "Cd27", "Tgfbr1", "Tgfbr2", "Tgfbr3", "Mki67")

DP1 = DotPlot(combined.dataset, features = genes_to_plot1, group.by = "seurat_clusters", cols = c("blue", "red")) + RotatedAxis() + theme(legend.position="bottom") +
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("Lineage-specific Receptors") + theme(plot.title = element_text(hjust = 0.5))

DP1

###Transcription factors
genes_to_plot2 = c("Tbx21", "Gata3", "Foxp3","Rorc","Pparg","Bhlhe40","Ahr","Ikzf2","Spi1","Irf4", "Eomes", "Zeb2","Stat1","Stat2","Stat3", "Stat4","Stat5a", "Stat5b", "Stat6", "Nfil3", "Junb", "Batf", "Klf2", "Prdm1", "Bcl6",
                   "Nfatc1", "Nfatc2", "Nfatc3", "Nfatc4", "Nfat5", "Foxo1", "Smad2", "Smad3", "Smad4","Zbtb16", "Hif1a","Nfkb1", "Nfkb2", "Maf", "Ikzf3")


DP2 = DotPlot(combined.dataset, features = genes_to_plot2, group.by = "seurat_clusters", cols = c("blue", "red")) + RotatedAxis() + theme(legend.position="bottom") +
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("Lineage-specific TFs") + theme(plot.title = element_text(hjust = 0.5))


Receptors_and_TFs_of_Th_lineages = DP1/DP2
Receptors_and_TFs_of_Th_lineages
Receptors_and_TFs_of_Th_lineages = ggsave("Dot plot for lineage-specific receptors and TFs_clusters.pdf",width = 15, height = 10, units = "in", dpi = 300)


####Using scCustomized

DP1_1 = DotPlot_scCustom(seurat_object = combined.dataset, features = genes_to_plot1, x_lab_rotate = TRUE)
DP1_2 = DotPlot_scCustom(seurat_object = combined.dataset, features = genes_to_plot2, x_lab_rotate = TRUE)

Receptors_and_TFs_of_Th_lineages = DP1_1/DP1_2
Receptors_and_TFs_of_Th_lineages
Receptors_and_TFs_of_Th_lineages = ggsave("Dot plot for lineage-specific receptors and TFs_clusters.pdf",width = 15, height = 10, units = "in", dpi = 300)

###Subsetting the data to exclude clusters 9,10,11,12,13
DefaultAssay(combined.dataset) = "integrated"

Subdata= subset(x=combined.dataset, idents = c(0,1,2,3,4,5,6,7,8))

####Reoder levels
levels(Subdata)
my_levels = c( "0" , "1", "2" , "3" , "4" , "5" , "6" , "7" , "8")
Subdata@meta.data$seurat_clusters = factor(x=Subdata@meta.data$seurat_clusters, levels = my_levels)

Umap_Sub = DimPlot(Subdata, reduction = "umap",group.by = 'seurat_clusters' ,label = TRUE) + NoLegend()
Umap_Sub
Umap_Sub = ggsave("Umap for subsetted data.pdf",width = 8, height =5, units = "in", dpi = 300)

###UMAP HTO
Umap_Sub_HTO = DimPlot(Subdata, reduction = "umap",group.by = 'HTO_classification' ,label = FALSE) + NoLegend()
Umap_Sub_HTO

View(Subdata@meta.data)

DefaultAssay(Subdata) = "RNA"

##Sccale Subsetted data
Subdata <- ScaleData(Subdata)
###Dotplots for for selected cytokines and trancription factors
genes_to_plot3 = c("Il2","Il4", "Il6","Il5","Il9","Il10","Il12a","Il13","Il17a","Il17f","Il21","Il22","Il24","Areg","Ifng","Tgfb1", "Csf2","Tnfsf11",
                   "Gzmb","Tbx21", "Gata3", "Pparg","Bhlhe40","Bcl6","Prdm1","Foxp3","Rorc","Ikzf2")
DP3 = DotPlot(Subdata, features = genes_to_plot3, group.by = "seurat_clusters", cols = c("blue", "red")) + RotatedAxis() + theme(legend.position="bottom") +
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("") + theme(plot.title = element_text(hjust = 0.5))


DP3

DP3 = ggsave("Dot plot for selected Cytokines and Tfs_Clusters.pdf",width = 8, height = 4, units = "in", dpi = 300)


###scCustomised
DP3_1 = DotPlot_scCustom(seurat_object = Subdata, features = genes_to_plot3, x_lab_rotate = TRUE)
DP3_1
DP3_1 = ggsave("Dot plot for selected Cytokines and Tfs_Clusters.pdf",width = 10, height = 4, units = "in", dpi = 300)



####Finding Marker genes
Combined.dataset_markers = FindAllMarkers(Subdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.dataset_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(x=Combined.dataset_markers, file = "Marker genes for subsetted clusters clusters_NEW.csv", quote = FALSE)

###A heatmap of the top 10 marker genes
top10 = Combined.dataset_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

Heatmap = DoHeatmap(Subdata, features = top10$gene, size = 8, angle =0)+ scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(axis.text.y = element_text(size = 15))+
  theme(axis.text = element_text(face="bold"))
Heatmap
Heatmap = ggsave("Heatmap for top 10 marker genes_Subsetted clusters.jpg",width = 18, height = 25, units = "in", dpi = 300)



####Finding markers distinguishing cluster 2 from all
cluster2.markers <- FindMarkers(Subdata, ident.1 = "2", min.pct = 0.25)
head(cluster2.markers, n = 5)
write.csv(x=cluster2.markers, file = "Marker genes for clusters 2 vs all.csv", quote = FALSE)

####Finding markers distinguishing cluster 5 from all
cluster5.markers <- FindMarkers(Subdata, ident.1 = "5", min.pct = 0.25)
head(cluster5.markers, n = 5)
write.csv(x=cluster5.markers, file = "Marker genes for clusters 5 vs all.csv", quote = FALSE)
####Finding markers distinguishing cluster 2 from 5
cluster2_from5.markers <- FindMarkers(Subdata, ident.1 = "2" , ident.2 =  "5", min.pct = 0.25)
head(cluster2_from5.markers, n = 5)
write.csv(x=cluster2_from5.markers, file = "Marker genes distinct to clusters 2 vs 5.csv", quote = FALSE)


#####Dot and feature plots for pathogenic Th2 marker genes
##Feature plots
Featureplot_for_pth2_markers = FeaturePlot(Subdata, features = c("Il1rl1","Gata3", "Bhlhe40" ,"Pparg", "Il4", "Il5", "Il13", "Tnfsf11", "Areg"), ncol = 5)
Featureplot_for_pth2_markers
Featureplot_for_pth2_markers = ggsave("Featureplot for pth2_markers_3.pdf",width = 15, height =6, units = "in", dpi = 300)
#Dotplots
genes_to_plot7 = c("Il1rl1","Il17rb","Gata3", "Bhlhe40" ,"Pparg", "Il4", "Il5", "Il13", "Tnfsf11", "Areg")
DP7 = DotPlot(Data, features = genes_to_plot7, group.by = "seurat_clusters", cols = c("blue", "red")) + RotatedAxis() + 
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("Pathogenic Th2 marker genes") + theme(plot.title = element_text(hjust = 0.5))
DP7

genes_to_plot7_1 = c("Cd44", "Sell","Il1rl1", "Il17rb", "Ptgdr2","Cd200r1","Pdcd1", "Il4", "Il5","Il6", "Il9", "Il13","Areg","Tnfsf8","Tnfsf11","Csf2", "Gata3","Pparg", "Bhlhe40")
DP7_1 = DotPlot(Data, features = genes_to_plot7_1, group.by = "HTO_classification", cols = c("blue", "red")) + RotatedAxis() + 
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("Pathogenic Th2 marker genes") + theme(plot.title = element_text(hjust = 0.5))
DP7_1
Pathogenic_Th2_marker_genes = DP7/DP7_1
Pathogenic_Th2_marker_genes
Pathogenic_Th2_marker_genes = ggsave("Dot plot Pathogenic Th2 marker genes.pdf",width = 10, height = 10, units = "in", dpi = 300)





#### Violin plot for for  Selected pTh2 genes Violin with boxplot

pTh2_Vln1 = VlnPlot(obj = Subdata, features = 'Il1rl1', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln2 = VlnPlot(obj = Subdata, features = 'Gata3', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0)  + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln3 = VlnPlot(obj = Subdata, features = 'Bhlhe40', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln4 = VlnPlot(obj = Subdata, features = 'Pparg', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0)  + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln5 = VlnPlot(obj = Subdata, features = 'Il4', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0)  + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln6 = VlnPlot(obj = Subdata, features = 'Il5', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln7 = VlnPlot(obj = Subdata, features = 'Il13', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln8 = VlnPlot(obj = Subdata, features = 'Tnfsf11', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2_Vln9 = VlnPlot(obj = Subdata, features = 'Areg', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())





## Apply grid.arrange function
pTh2_Plots = grid.arrange(pTh2_Vln1, pTh2_Vln2, pTh2_Vln3, pTh2_Vln4, pTh2_Vln5, pTh2_Vln6, pTh2_Vln7, pTh2_Vln8, pTh2_Vln9, ncol = 1)
pTh2_Plots
pTh2_Plots = ggsave("Violin plot for T helper subsets.pdf",width = 6, height =12, units = "in", dpi = 300)




###Dot and Feature Plots for Identification of pathogenic Th2 cells_seurat clusters
genes_to_plot4 = c("Cd44", "Il1rl1","Cxcr6", "Pdcd1","Cd200r1","Crlf2","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Il17rb", 'Il17re',
                   "Ifnar1","Ifnar2","Cd69","Klrg1","Cd27", "S1pr1", "Il9r","Ffar3","Cxcr3",'Cxcr4',"Cxcr5",'Slamf1','Slamf6',"Il23r")
DP4 = DotPlot(Subdata, features = genes_to_plot4, group.by = "seurat_clusters", cols = c("blue", "red")) + RotatedAxis() + 
  # set transparency
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA)
  ) + xlab(NULL) + ylab(NULL) + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="bottom")
DP4
DP4 = ggsave("Dot plot Pathogenic Th2 Identification.pdf",width =10, height =5, units = "in", dpi = 300)


###scCustomised
DP4_1 = DotPlot_scCustom(seurat_object = Subdata, features = genes_to_plot4, x_lab_rotate = TRUE)
DP4_1
DP4_1 = ggsave("Dot plot Pathogenic Th2 Identification.pdf",width = 10, height = 4, units = "in", dpi = 300)




#### Violin plot for for Identification of pth2 cells Violin with boxplot

pTh2.1_Vln1 = VlnPlot(obj = Subdata, features = 'Cd44', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2.1_Vln2 = VlnPlot(obj = Subdata, features = 'Il1rl1', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0)  + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2.1_Vln3 = VlnPlot(obj = Subdata, features = 'Cd27', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2.1_Vln4 = VlnPlot(obj = Subdata, features = 'Klrg1', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0)  + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2.1_Vln5 = VlnPlot(obj = Subdata, features = 'Pdcd1', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0)  + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
pTh2.1_Vln6 = VlnPlot(obj = Subdata, features = 'Cd69', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())

## Apply grid.arrange function
Combined_Plots = grid.arrange(pTh2.1_Vln1, pTh2.1_Vln2, pTh2.1_Vln3, pTh2.1_Vln4, pTh2.1_Vln5, pTh2.1_Vln6, ncol = 3)
Combined_Plots
Combined_Plots = ggsave("Violin plot for T helper subsets.pdf",width = 6, height =12, units = "in", dpi = 300)


####Violin plot for Trm-genes
Trm1 = VlnPlot(obj = Subdata, features = 'Cd69', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Trm2 = VlnPlot(obj = Subdata, features = 'Cxcr6', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Trm3 = VlnPlot(obj = Subdata, features = 'Itgae', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Trm4 = VlnPlot(obj = Subdata, features = 'Itga1', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Trm5 = VlnPlot(obj = Subdata, features = 'Itgal', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Trm6 = VlnPlot(obj = Subdata, features = 'S1pr1', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Trm7 = VlnPlot(obj = Subdata, features = 'Klf2', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Trm8 = VlnPlot(obj = Subdata, features = 'Ccr7', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())

## Apply grid.arrange function
Combined_Plots_Trm = grid.arrange(Trm1, Trm2, Trm3, Trm4, Trm5, Trm6, Trm7, Trm8, ncol = 2)

Combined_Plots_Trm = ggsave("Violin plot for T helper subsets.pdf",width = 6, height =12, units = "in", dpi = 300)


####Assigning cell type identity to clusters (Subsetted data)
levels(Subdata)
Idents(object = Subdata)
new.cluster.ids = c( "Naive 1 (0)", "Naive 2 (1)", "Dying cells (11)", "pTh2 (2)", "Th1 (3)", "Th17 (4)", "Cd69 high pTh2 (5)", "Th2/Treg (6)", "tTreg (7)", "Naive 3 (8)", "Proliferating cells (9)")
names(new.cluster.ids) <- levels(Subdata)
Subdata<- RenameIdents(Subdata, new.cluster.ids)

# Stash cell identity classes
Subdata[["Cell_type"]] <- Idents(object = Subdata)
Subdata<- StashIdent(object = Subdata, save.name = "Cell_type")




Dim1 = DimPlot(Subdata, reduction = "umap", label = FALSE, group.by = 'Cell_type') + NoLegend()
Dim1
Dim1 = ggsave("Feature Plot_Clusters_With names 3.pdf",width = 6, height =5, units = "in", dpi = 300)

Dim2 = DimPlot(Subdata, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters") + NoLegend()
Dim2
Dimplots = Dim2|Dim1
Dimplots

Dimplots = ggsave("Clusters_With names 2 .pdf",width = 8, height =5, units = "in", dpi = 300)


###Making stacked vln plots
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(Subdata, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
features <-c("Cd69", "Cxcr6", "Itgae", "Itga1",'Itgal', "S1pr1", "Klf2", "Ccr7")
features1<- c("Il1rl1","Gata3", "Bhlhe40" ,"Pparg", "Il4", "Il5", "Il13", "Tnfsf11", "Areg")
features <- c("Mki67","Sell","Cd44","Cd69","Pdcd1", "Cxcr6" ,'Cxcr3',"Klf2", "S1pr1","Ctla4","Icos" ,"Cdh1", "Itgae", "Itga1", 'Itgal', 'Ccr7')
features <- c("Il4","Il5","Il10","Il13", "Areg" ,"Tnfsf11", "Tgfb1")

Vln1 = StackedVlnPlot(obj = Subdata, features = features, group.by = "seurat_clusters")
Vln1
Vln1 = ggsave("Violin plot for TRM genes.pdf",width = 6, height =14, units = "in", dpi = 300)

Vln2 = StackedVlnPlot(obj = Subdata, features = features, group.by = "Cell_type")
Vln2
Vln2= ggsave("Violin plot for pathogenic Th2 genes.pdf",width = 6, height =14, units = "in", dpi = 300)

Vln3= StackedVlnPlot(obj = Subdata, features = features1, group.by = "seurat_clusters")
Vln3
Vln3= ggsave("Violin plot for pathogenic Th2 genes.pdf",width = 6, height =14, units = "in", dpi = 300)

Vln4= StackedVlnPlot(obj = Subdata, features = features, group.by = "HTO_classification")
Vln4
Vln4= ggsave("Violin plot for Pathogenic Th2 genes_HTO_classification.pdf",width = 6, height =12, units = "in", dpi = 300)

Vln5= StackedVlnPlot(obj = Subdata, features = features, group.by = "seurat_clusters")
Vln5
Vln5= ggsave("Violin plot for marker genes distinguishing lung T helper cells genes.pdf",width = 6, height =12, units = "in", dpi = 300)





##Saving the metadata

view(Subdata@meta.data)
write.csv(combined.dataset@meta.data,"metadata.csv")


#### Violin plot for T helper subset genes
Th_Vln1 = VlnPlot(obj = Subdata, features = 'Klf2', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) 
Th_Vln2 = VlnPlot(obj = Subdata, features = 'Ccr7', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln3 = VlnPlot(obj = Subdata, features = 'Sell', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln4 = VlnPlot(obj = Subdata, features = 'Cd44', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln5 = VlnPlot(obj = Subdata, features = 'Ifng', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln6 = VlnPlot(obj = Subdata, features = 'Tbx21', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln7 = VlnPlot(obj = Subdata, features = 'Ccl5', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln8 = VlnPlot(obj = Subdata, features = 'Cxcr3', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln9 = VlnPlot(obj = Subdata, features = 'Il4', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln10 = VlnPlot(obj = Subdata, features = 'Il5', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln11 = VlnPlot(obj = Subdata, features = 'Il13', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln12 = VlnPlot(obj = Subdata, features = 'Gata3', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln13 = VlnPlot(obj = Subdata, features = 'Foxp3', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln14 = VlnPlot(obj = Subdata, features = 'Ikzf2', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln15 = VlnPlot(obj = Subdata, features = 'Klrg1', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln16 = VlnPlot(obj = Subdata, features = 'Il2ra', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln17 = VlnPlot(obj = Subdata, features = 'Il17a', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln18 = VlnPlot(obj = Subdata, features = 'Il17re', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln19 = VlnPlot(obj = Subdata, features = 'Ccr6', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln20 = VlnPlot(obj = Subdata, features = 'Rorc', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())


## Apply grid.arrange function
Combined_Plots = grid.arrange(Th_Vln1, Th_Vln2, Th_Vln3, Th_Vln4, Th_Vln5, Th_Vln6, Th_Vln7, Th_Vln8, Th_Vln9,Th_Vln10,
                              Th_Vln11, Th_Vln12, Th_Vln13, Th_Vln14, Th_Vln15, Th_Vln16, Th_Vln17, Th_Vln18, Th_Vln19,Th_Vln20,ncol = 4)

Combined_Plots = ggsave("Violin plot for T helper subsets.pdf",width = 6, height =12, units = "in", dpi = 300)


#### Violin plot for T helper subset genes without dots
Th_Vln1 = VlnPlot(obj = Subdata, features = 'Klf2', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln2 = VlnPlot(obj = Subdata, features = 'Ccr7', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln3 = VlnPlot(obj = Subdata, features = 'Sell', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln4 = VlnPlot(obj = Subdata, features = 'Cd44', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln5 = VlnPlot(obj = Subdata, features = 'Ifng', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln6 = VlnPlot(obj = Subdata, features = 'Tbx21', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln7 = VlnPlot(obj = Subdata, features = 'Ccl5', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln8 = VlnPlot(obj = Subdata, features = 'Cxcr3', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln9 = VlnPlot(obj = Subdata, features = 'Il4', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln10 = VlnPlot(obj = Subdata, features = 'Il5', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln11 = VlnPlot(obj = Subdata, features = 'Il13', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln12 = VlnPlot(obj = Subdata, features = 'Il1rl1', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln13 = VlnPlot(obj = Subdata, features = 'Foxp3', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln14 = VlnPlot(obj = Subdata, features = 'Ikzf2', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln15 = VlnPlot(obj = Subdata, features = 'Klrg1', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln16 = VlnPlot(obj = Subdata, features = 'Ctla4', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln17 = VlnPlot(obj = Subdata, features = 'Il17a', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln18 = VlnPlot(obj = Subdata, features = 'Il17re', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln19 = VlnPlot(obj = Subdata, features = 'Ccr6', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())
Th_Vln20 = VlnPlot(obj = Subdata, features = 'Rorc', group.by = "seurat_clusters", pt.size = 0) + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line())


## Apply grid.arrange function
Combined_Plots = grid.arrange(Th_Vln1, Th_Vln2, Th_Vln3, Th_Vln4, Th_Vln5, Th_Vln6, Th_Vln7, Th_Vln8, Th_Vln9,Th_Vln10,
                              Th_Vln11, Th_Vln12, Th_Vln13, Th_Vln14, Th_Vln15, Th_Vln16, Th_Vln17, Th_Vln18, Th_Vln19,Th_Vln20,ncol = 4)

Combined_Plots = ggsave("Violin plot for T helper subsets.pdf",width = 6, height =12, units = "in", dpi = 300)





Combined_Plots = ggsave("Violin plot for AP-1 genes.pdf",width = 6, height =12, units = "in", dpi = 300)
####Subsetting cluster 2 to distinguish WT and KO
Subdata_cluster_2= subset(x=Subdata, idents = "2")


##Scale Subsetted data
Subdata_cluster_2 <- ScaleData(Subdata_cluster_2)
View(Subdata_cluster_2@meta.data)

##add genotype to column in metadata
metadata = Subdata_cluster_2@meta.data
metadata = metadata%>%mutate(genotype=ifelse(grepl("WT", HTO_classification), "WT", "HDAC1cKO"))
Subdata_cluster_2@meta.data = metadata

###Set Idents to genotypes
Idents(Subdata_cluster_2) = "genotype"
Idents(Subdata_cluster_2)
levels(Subdata_cluster_2) <- c("WT", "HDAC1cKO")
levels(Subdata_cluster_2)
##Finding the DE markers between two genotypes
Markergenesfor_Cluster2 = FindMarkers(Subdata_cluster_2, ident.1 = "HDAC1cKO", ident.2 = "WT", min.pct = 0.25)

head(Markergenesfor_Cluster2)
write.csv(x=Markergenesfor_Cluster2, file = "Marker genes distinct to clusters 2_HDAC1-cKO compared WT.csv", quote = FALSE)

###Vln for selected genes without
Vln1 = VlnPlot(obj = Subdata_cluster_2, features = 'Il4', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln2 = VlnPlot(obj = Subdata_cluster_2, features = 'Il5', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln3 = VlnPlot(obj = Subdata_cluster_2, features = 'Il13', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln4 = VlnPlot(obj = Subdata_cluster_2, features = 'Tnfsf11', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln5 = VlnPlot(obj = Subdata_cluster_2, features = 'Tgfb1', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln6 = VlnPlot(obj = Subdata_cluster_2, features = 'Areg', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln7 = VlnPlot(obj = Subdata_cluster_2, features = 'Furin', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln8 = VlnPlot(obj = Subdata_cluster_2, features = 'Gata3', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln9 = VlnPlot(obj = Subdata_cluster_2, features = 'Pparg', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln10 = VlnPlot(obj = Subdata_cluster_2, features = 'Bhlhe40', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")




## Apply grid.arrange function
Combined_Plots = grid.arrange(Vln1, Vln2, Vln3, Vln4, Vln5, Vln6, Vln7, Vln8, Vln9,Vln10, ncol = 5)



###Vln for selected genes with p values
Vln1 = VlnPlot(obj = Subdata_cluster_2, features = 'Il4', group.by = "genotype", cols = c("red", "blue"), y.max = 7) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln2 = VlnPlot(obj = Subdata_cluster_2, features = 'Il5', group.by = "genotype", cols = c("red", "blue"), y.max = 7) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln3 = VlnPlot(obj = Subdata_cluster_2, features = 'Il13', group.by = "genotype", cols = c("red", "blue"), y.max = 7) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln4 =  VlnPlot(obj = Subdata_cluster_2, features = 'Tnfsf11', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln5 = VlnPlot(obj = Subdata_cluster_2, features = 'Tgfb1', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln6 = VlnPlot(obj = Subdata_cluster_2, features = 'Areg', group.by = "genotype", cols = c("red", "blue"), y.max = 8) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln7 = VlnPlot(obj = Subdata_cluster_2, features = 'Furin', group.by = "genotype", cols = c("red", "blue"), y.max = 7) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln8 = VlnPlot(obj = Subdata_cluster_2, features = 'Gata3', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln9 = VlnPlot(obj = Subdata_cluster_2, features = 'Pparg', group.by = "genotype", cols = c("red", "blue"), y.max = 5) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln10 = VlnPlot(obj = Subdata_cluster_2, features = 'Bhlhe40', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))




## Apply grid.arrange function
Combined_Plots = grid.arrange(Vln1, Vln2, Vln3, Vln4, Vln5, Vln6, Vln7, Vln8, Vln9,Vln10, ncol = 5)


####Subsetting cluster 5 to distinguish WT and KO
Subdata_cluster_5= subset(x=Subdata, idents = "5")
Subdata_cluster_5 <- ScaleData(Subdata_cluster_5)
View(Subdata_cluster_5@meta.data)

##add genotype to column in metadata
metadata2 = Subdata_cluster_5@meta.data
metadata2 = metadata2%>%mutate(genotype=ifelse(grepl("WT", HTO_classification), "WT", "HDAC1cKO"))
Subdata_cluster_5@meta.data = metadata2

###Set Idents to genotypes
Idents(Subdata_cluster_5) = "genotype"
Idents(Subdata_cluster_5)
levels(Subdata_cluster_5) <- c("WT", "HDAC1cKO")
levels(Subdata_cluster_5)
##Finding the DE markers between two genotypes
Markergenesfor_Cluster5 = FindMarkers(Subdata_cluster_5, ident.1 = "HDAC1cKO", ident.2 = "WT", min.pct = 0.25)
head(Markergenesfor_Cluster5)
write.csv(x=Markergenesfor_Cluster5, file = "Marker genes distinct to clusters 5_HDAC1-cKO compared to WT.csv", quote = FALSE)



###Vln for selected genes without p values
Vln2_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Il5', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln3_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Il13', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln4_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Tnfsf11', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln5_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Tgfb1', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln6_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Areg', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln7_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Furin', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln8_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Gata3', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln9_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Pparg', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln10_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Bhlhe40', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")




## Apply grid.arrange function
Combined_Plots = grid.arrange(Vln1_1, Vln2_1, Vln3_1, Vln4_1, Vln5_1, Vln6_1, Vln7_1, Vln8_1, Vln9_1,Vln10_1, ncol = 5)

Combined_Plots = ggsave("Violin plot for selected cluster 2 genes_WT and KO.pdf",width = 6, height =12, units = "in", dpi = 300)


###Vln for selected genes with p values
Vln1 = VlnPlot(obj = Subdata_cluster_5, features = 'Il4', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln2 = VlnPlot(obj = Subdata_cluster_5, features = 'Il5', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln3 = VlnPlot(obj = Subdata_cluster_5, features = 'Il13', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln4 =  VlnPlot(obj = Subdata_cluster_5, features = 'Tnfsf11', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln5 = VlnPlot(obj = Subdata_cluster_5, features = 'Tgfb1', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln6 = VlnPlot(obj = Subdata_cluster_5, features = 'Areg', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln7 = VlnPlot(obj = Subdata_cluster_5, features = 'Furin', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln8 = VlnPlot(obj = Subdata_cluster_5, features = 'Gata3', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln9 = VlnPlot(obj = Subdata_cluster_5, features = 'Pparg', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))

Vln10 = VlnPlot(obj = Subdata_cluster_5, features = 'Bhlhe40', group.by = "genotype", cols = c("red", "blue"), y.max = 6) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none") + 
  geom_boxplot(width=0.1, fill='white', outlier.colour = NULL, outlier.size = 0) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) + 
  stat_compare_means(comparisons = list(c("WT", "HDAC1cKO", method = "wilcox.test")))




## Apply grid.arrange function
Combined_Plots = grid.arrange(Vln1, Vln2, Vln3, Vln4, Vln5, Vln6, Vln7, Vln8, Vln9,Vln10, ncol = 5)


#####AP-1 protein

Vln1_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Jun', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO")) + theme(legend.position="none")
Vln2_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Junb', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln3_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Jund', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln4_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Fos', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln5_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Fosb', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln6_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Fosl2', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln7_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Zfp36l2', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")
Vln8_1_1 = VlnPlot(obj = Subdata_cluster_5, features = 'Rinl', group.by = "genotype", cols = c("red", "blue")) + scale_x_discrete(limits=c("WT", "HDAC1cKO"))+ theme(legend.position="none")

## Apply grid.arrange function
Combined_Plots = grid.arrange(Vln1_1_1, Vln2_1_1, Vln3_1_1, Vln4_1_1, Vln5_1_1, Vln6_1_1, Vln7_1_1, Vln8_1_1, ncol = 4)

Combined_Plots = ggsave("Violin plot for selected cluster 2 genes_WT and KO.pdf",width = 6, height =12, units = "in", dpi = 300)






###Heatmap for HALLMARK_IL-2/STAT5 and TNFA-NFKB pathways leading-edge genes

genes_heatmap = c('Ahnak', 'Areg',	'Atf3',	'Batf',	'Bcl2a1d',	'Bcl2l1',	'Bhlhe40',	'Capg',	'Ccnd2',	'Ccr4',	'Cd44',	'Cdkn1a',	'Cish',	'Ctla4',	'Dennd5a',	'Dhrs3',	'Dusp1',	'Dusp4',	'Dusp5',	'Ehd1',	'Fgl2',	
                  'Fosb',	'Fosl2',	'Furin',	'Gadd45b',	'Gpr65',	'Gpx4',	'Icos',	'Id2',	'Ier5',	'Ifitm3',	'Il10',	'Il10ra',	'Il13',	'Il1rl1',	'Il2ra',	'Irf4',	'Itgav',	'Junb',	'Kdm6b',	'Litaf',	'Map2k3',	
                  'Mapkapk2',	'Myo1e',	'Ndrg1',	'Nfil3',	'Nfkb1',	'Nfkbia',	'Nfkbiz',	'Nr4a1',	'Nr4a2',	'Nr4a3',	'Nrp1',	'Nt5e',	'Odc1',	'Phlda1',	'Pim1',	'Plin2',	'Rel',	'Rgs16',	'Rnf19b',	'Rora',	
                  'Serpinb6a',	'Socs2',	'Tnfaip3',	'Tnfrsf18',	'Tnfrsf1b',	'Tnfrsf4',	'Tnfrsf9',	'Tnfsf11',	'Traf1',	'Zc3h12a')

DoHeatmap(Subdata, features = genes_heatmap, size = 10) + scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(axis.text.y = element_text(size = 15))+
  theme(axis.text = element_text(face="bold")) 

ggsave("Heatmap of IL-2 STAT5 and TNFa-Nfkb pathways genes_Seurat.jpg",width = 15, height =20, units = "in", dpi = 300)





###Feature plots for AP-1, MAPK, NF-KB

genes_for_Featureplots1= c('Fos','Fosb',  'Jun')

genes_for_Featureplots1_1= 'Fos'
genes_for_Featureplots1_2= 'Fosb'
genes_for_Featureplots1_3= 'Jun'


# Create Plots using scCustomize for Fos, Fosb and Jun

FP1 = FeaturePlot_scCustom(seurat_object = Subdata, features = genes_for_Featureplots1_1, order = F)
FP2 = FeaturePlot_scCustom(seurat_object = Subdata, features = genes_for_Featureplots1_2, order = F)
FP3 = FeaturePlot_scCustom(seurat_object = Subdata, features = genes_for_Featureplots1_3, order = F)


FP4 = DimPlot(Subdata, reduction = "umap",group.by = 'seurat_clusters' ,label = TRUE) + NoLegend()

Combined_Plots = grid.arrange(FP4, FP1, FP2, FP3, ncol = 2, nrow = 2) 



# Create Plots using scCustomize for TNF, MAPK, AP1, NFKB

genes_for_Featureplots2= c('Tnfrsf4', 'Tnfrsf9', 'Tnfrsf18', 'Traf1', 'Map2k3',	'Mapkapk2','Fosl2','Batf', 'Nfkb1','Rel') 

FeaturePlot_scCustom(seurat_object = Subdata, features = genes_for_Featureplots2, order = F) +
  patchwork::plot_layout(ncol = 2)


####Showing notch4, Rbpj, FoxP3, Gata3

Vln_notch4 = VlnPlot(obj = Subdata, features = 'Notch4', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) 

Vln_Rbpj = VlnPlot(obj = Subdata, features = 'Rbpj', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) 

Vln_FoxP3 = VlnPlot(obj = Subdata, features = 'Foxp3', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) 

Vln_Gata3 = VlnPlot(obj = Subdata, features = 'Gata3', group.by = "seurat_clusters") + theme(legend.position="none") + theme(axis.text.x=element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_line()) 


## Apply grid.arrange function
Combined_Plots_Rbpj_Notch4 = grid.arrange(Vln_notch4, Vln_Rbpj,Vln_FoxP3,Vln_Gata3, ncol = 2)

###Saving SessionInfo

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")




###PART 3-GSEA analysis of lung pTh2 vs Day 15 Th2 from Tibbitt et al----

library(fgsea)
library(msigdbr)
library(data.table)
library(ggplot2)
library(readxl)
library(stringr)
library(dplyr)


####GSEA of lung pTh2 vs Day 15 Th2 (Tibbitt et al. 2019)
hallmark <- fgsea::gmtPathways("~/Desktop/DATA/scRNA seq/scRNAseq_Integration/New Analysis/scRNAseq Analyses/GSEA (Comparison with Day 15 Th2 from Tibbitt et al)/Day_15_Th2_Cells_1.gmt")

head(hallmark)


peTh2 = read_xlsx("peTh2_genelist.xlsx")
head(peTh2)



###Organize and rank genes by log2FC

peTh2 <- peTh2 %>% arrange(desc(avg_log2FC))
genelist1 <- peTh2$avg_log2FC
names(genelist1) = peTh2$gene
fold_changes_peTh2 = genelist1



fgseaRes_peTh2<- fgsea(pathways = hallmark, 
                       stats    = fold_changes_peTh2,
                       minSize  = 1,
                       maxSize  = Inf)
View(fgseaRes_peTh2)

plotEnrichment(hallmark[["Day_15_Th2_genelist"]],
               fold_changes_peTh2) + labs(title="peTh2 vs Day 15 Th2")

Enrichmentplot = ggsave("Comparison of peTh2 to Day 15 Th2 from Tibbitt et al.pdf",width = 5, height =3, units = "in", dpi = 300)


### collapse the leading edge gene list to vector
fgseaRes_peTh2$leadingEdge <- vapply(fgseaRes_peTh2$leadingEdge, paste, collapse = ",", character(1L))

stats1 <- data.frame(pathway=fgseaRes_peTh2$pathway) %>%
  data.frame(pvalue=fgseaRes_peTh2$pval) %>%
  data.frame(adjpvalue=fgseaRes_peTh2$padj) %>%
  data.frame(Size=fgseaRes_peTh2$size) %>%
  data.frame(NES=fgseaRes_peTh2$NES) %>%
  data.frame(genes=fgseaRes_peTh2$leadingEdge)

write.csv(stats1, file = "~/Desktop/DATA/scRNA seq/scRNAseq_Integration/New Analysis/scRNAseq Analyses/GSEA (Comparison with Day 15 Th2 from Tibbitt et al)/peTh2 vs Day 15 Th2 from Tibbitt et al 2019.csv", quote = TRUE, col.names = TRUE, row.names = FALSE) #quote=T to avoid splitting by ','



#####GSEA comparison of In vivo Th2 Trm with Day 15 Th2 (Tibbitt et al. 2019)
Th2_Trm = read_xlsx("Th2 Trm_genelist.xlsx")
head(Th2_Trm)



###Organize and rank genes by log2FC

Th2_Trm <- Th2_Trm %>% arrange(desc(avg_log2FC))
genelist1 <- Th2_Trm$avg_log2FC
names(genelist1) = Th2_Trm$gene
fold_changes_Th2_Trm = genelist1



fgseaRes_Th2_Trm<- fgsea(pathways = hallmark, 
                         stats    = fold_changes_Th2_Trm,
                         minSize  = 1,
                         maxSize  = Inf)
View(fgseaRes_Th2_Trm)

plotEnrichment(hallmark[["Day_15_Th2_genelist"]],
               fold_changes_Th2_Trm) + labs(title="Th2 Trm vs Day 15 Th2")

Enrichmentplot = ggsave("Comparison of Th2 Trm to Day 15 Th2 from Tibbitt et al.pdf",width = 5, height =3, units = "in", dpi = 300)


###collapse the leading edge gene list to vector
fgseaRes_Th2_Trm$leadingEdge <- vapply(fgseaRes_Th2_Trm$leadingEdge, paste, collapse = ",", character(1L))

stats1 <- data.frame(pathway=fgseaRes_Th2_Trm$pathway) %>%
  data.frame(pvalue=fgseaRes_Th2_Trm$pval) %>%
  data.frame(adjpvalue=fgseaRes_Th2_Trm$padj) %>%
  data.frame(Size=fgseaRes_Th2_Trm$size) %>%
  data.frame(NES=fgseaRes_Th2_Trm$NES) %>%
  data.frame(genes=fgseaRes_Th2_Trm$leadingEdge)

write.csv(stats1, file = "~/Desktop/DATA/scRNA seq/scRNAseq_Integration/New Analysis/scRNAseq Analyses/GSEA (Comparison with Day 15 Th2 from Tibbitt et al)/Th2 Trm vs Day 15 Th2 from Tibbitt et al 2019.csv", quote = TRUE, col.names = TRUE, row.names = FALSE) #quote=T to avoid splitting by ','







###PART 4- Volcano plots for comparison of different scRNA-seq clusters----

library(EnhancedVolcano)
library(ggplot2)
library(readxl)







####Cluster 2 vs all
Clu2vsall = read.csv("Marker genes for clusters 2 vs all.csv")
head(Clu2vsall)
####


EnhancedVolcano(Clu2vsall,
                lab = Clu2vsall$X,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = c('Klf2', 'Klf3','Satb1', 'Lef1','Txk','S1pr1','Gimap6', 'Gimap3', 'Gramd3', 'Il4',
                              'Il5', 'Areg',  'Gadd45b','Zeb2', 'Pparg','Map2k3', 'Nfat5','Fos', 'Fosb', 'Fosl2',
                              'Tnfsf11','Mt1', 'Il10', 'Il1rl1', 'Il13', 'Bhlhe40', 'Gata3', 'Nfkbia','Ccr7'),
                pCutoff = 0.05,
                FCcutoff =0.58,
                pointSize = 2,
                labSize = 14,
                legendLabels=c("NS", "avg_log2FC", "p_val_adj", "p_val_adj & avg_log2FC"),
                title = "Cluster all/cluster 2", subtitle = '',
                xlab = bquote("avg_log"[2]*"FC"),
                ylab = bquote("-log"[10]*"(padj)"),
                legendPosition = "",
                legendLabSize = 14,
                boxedLabels = FALSE,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1,
                borderColour = 'black',
                legendIconSize = 4,
                axisLabSize = 40,
                cutoffLineType = 'solid',
                cutoffLineWidth = 1)
####Saving the data
ggsave("Clu2vsall_Volcanoplot.pdf", units="in", width=15, height=15, dpi=300)

####Cluster 5 vs all
Clu5vsall = read.csv("Marker genes for clusters 5 vs all.csv")
head(Clu5vsall)
####


EnhancedVolcano(Clu5vsall,
                lab = Clu5vsall$X,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = c('Ifngr1','Klf2','Satb1','S1pr1','Gimap6', 'Gimap3', 'Gramd3','Ccr7', 'Zeb1','Odc1',
                              'Il5', 'Fos', 'Fosb', 'Fosl2','Gadd45b','Il13', 'Cd69','Lpar6', 'Pprag', 'Maf', 'Il4',
                              'Tnfsf11', 'Il1rl1', 'Il13', 'Bhlhe40', 'Gata3',  'Areg', 'Il13', 'Mns1',
                              'Furin', 'Jun', 'Klf6', 'Rbpj', 'Dusp1', 'Egr1', 'Zeb2', 'Klf6', 'Map2k3', 'Nfat5'),
                pCutoff = 0.05, 
                FCcutoff =0.58,
                pointSize = 2,
                labSize = 14,
                legendLabels=c("NS", "avg_log2FC", "p_val_adj", "p_val_adj & avg_log2FC"),
                title = "Cluster all/cluster 5", subtitle = '',
                xlab = bquote("avg_log"[2]*"FC"),
                ylab = bquote("-log"[10]*"(padj)"),
                legendPosition = "",
                legendLabSize = 14,
                boxedLabels = FALSE,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1,
                borderColour = 'black',
                legendIconSize = 4,
                axisLabSize = 40,
                cutoffLineType = 'solid',
                cutoffLineWidth = 1)
####Saving the data
ggsave("Clu5vsall_volcanoplot.pdf", units="in", width=15, height=15, dpi=300)

####Cluster 2 vs cluster 5
Clu2vs5 = read.csv("Marker genes distinct to clusters 2 vs 5.csv")
head(Clu2vs5)
####


EnhancedVolcano(Clu2vs5,
                lab = Clu2vs5$X,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = c('Il5', 'Fos', 'Fosb', 'Il13', 'Cd69','Lpar6', 'Il4',
                              'Il1rl1', 'Il13',  'Gata3', 'Tnfrsf4','Tnfrsf9' ,'Areg', 'Il13',
                              'Tnfrsf18', 'Furin', 'Jun', 'Egr1', 'Pdcd1', 'Calca', 'Il10',
                              'Zfp36l2', 'Slc38a2', 'Nfkb1', 'Nfkb2', 'Map2k3', 'Map2k6', 'Zeb2', 'Ccr7'),
                pCutoff = 0.05,
                FCcutoff =0.58,
                pointSize = 2,
                labSize = 14,
                legendLabels=c("NS", "avg_log2FC", "p_val_adj", "p_val_adj & avg_log2FC"),
                title = "Cluster 5/cluster 2", subtitle = '',
                xlab = bquote("avg_log"[2]*"FC"),
                ylab = bquote("-log"[10]*"(padj)"),
                legendPosition = "",
                legendLabSize = 14,
                boxedLabels = FALSE,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1.5,
                borderColour = 'black',
                legendIconSize = 4,
                axisLabSize = 40,
                cutoffLineType = 'solid',
                cutoffLineWidth = 1)
####Saving the data
ggsave("Clu2vs5_Volcanoplot.pdf", units="in", width=15, height=15, dpi=300)

##




###PART 5- Gene set enrichment analysis for pathways----

library(fgsea)
library(msigdbr)
library(data.table)
library(ggplot2)
library(readxl)
library(stringr)
library(dplyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
library(pathview)
library(VennDiagram)
library(clusterProfiler)
library(enrichplot)

hallmark <- fgsea::gmtPathways("~/Desktop/DATA/scRNA seq/scRNAseq_Integration/New Analysis/mh.all.v0.3.symbols.gmt") #50 gene sets


head(hallmark)

#####GSEA ANALYSIS_CLUSTER 2 vs all clusters
Cluster2_GSEA = read.csv ("Marker genes for clusters 2 vs all.csv")
head(Cluster2_GSEA)



### Organise and rank genes  by log2FC

Cluster2_GSEA <- Cluster2_GSEA %>% arrange(desc(avg_log2FC))
genelist <- Cluster2_GSEA$avg_log2FC
names(genelist) = Cluster2_GSEA$X
fold_changes_cluster2 = genelist



###Running GSEA to determine genesets that are upregulated or downregualted in cluster 2 vs all clusters


fgseaRes_Cluster2 <- fgsea(pathways = hallmark, 
                           stats    = fold_changes_cluster2,
                           minSize  = 15,
                           maxSize  = 500)
View(fgseaRes_Cluster2)

head(fgseaRes_Cluster2[order(pval), ])

#### up/down regulated pathways are generated using the log10(pvalue).
filtered_pathway_cluster2 <- subset(fgseaRes_Cluster2, pval < 0.05)

filt_up_cluster2 <- subset(filtered_pathway_cluster2, NES > 0.0)
filt_up_cluster2$log10 <- -log10(filt_up_cluster2$pval)
filt_up_cluster2$pathway <- str_replace_all(filt_up_cluster2$pathway, c("HALLMARK"="", "_"=" "))



Cluster_up_Pathways = ggplot(data=filt_up_cluster2, aes(x=reorder(pathway,log10), y=log10)) + 
  geom_bar(stat="identity", fill="red", width = 0.7) +
  coord_flip() +
  labs(x="", y="pvalue (-log10)") +
  theme_minimal() + ggtitle("Enriched pathways in Cluster 2 vs all")
Cluster_up_Pathways = ggsave("Top pathways for cluster 2.pdf",width = 8, height =4, units = "in", dpi = 300)



### collapse the leading edge gene list to vector
fgseaRes_Cluster2$leadingEdge <- vapply(fgseaRes_Cluster2$leadingEdge, paste, collapse = ",", character(1L))

stats1 <- data.frame(pathway=fgseaRes_Cluster2$pathway) %>%
  data.frame(pvalue=fgseaRes_Cluster2$pval) %>%
  data.frame(adjpvalue=fgseaRes_Cluster2$padj) %>%
  data.frame(NES=fgseaRes_Cluster2$NES) %>%
  data.frame(genes=fgseaRes_Cluster2$leadingEdge)

write.csv(stats1, file = "~/Desktop/DATA/scRNA seq/scRNAseq_Integration/New Analysis/Cluster 2 Pathways.csv", quote = TRUE, col.names = TRUE, row.names = FALSE) #quote=T to avoid splitting by ','




###Making enrichment plot for a pathway:

##IL2 STAT5
plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]],
               fold_changes_cluster2) + labs(title="IL2 STAT5 SIGNALING")

Enrichmentplot = ggsave("IL2 STAT5 pathway for cluster 2.pdf",width = 8, height =4, units = "in", dpi = 300)

###TNFA
plotEnrichment(hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               fold_changes_cluster2) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")

Enrichmentplot = ggsave("TNFA pathway pathway for cluster 2.pdf",width = 8, height =4, units = "in", dpi = 300)



#####GSEA ANALYSIS_CLUSTER 5 vs all clusters
Cluster5_GSEA = read.csv ("Marker genes for clusters 5 vs all.csv")
head(Cluster5_GSEA)



###Organise and rank genes  by log2FC

Cluster5_GSEA <- Cluster5_GSEA %>% arrange(desc(avg_log2FC))
genelist <- Cluster5_GSEA$avg_log2FC
names(genelist) = Cluster5_GSEA$X
fold_changes_cluster5 = genelist



###Running GSEA to determine genesets that are upregulated or downregualted

fgseaRes_Cluster5 <- fgsea(pathways = hallmark, 
                           stats    = fold_changes_cluster5,
                           minSize  = 15,
                           maxSize  = 500)

head(fgseaRes_Cluster5[order(pval), ])
View(fgseaRes_Cluster5)

#### up/down regulated pathways are generated using the log10(pvalue).
filtered_pathway_cluster5 <- subset(fgseaRes_Cluster5, pval < 0.05)

filt_up_cluster5 <- subset(filtered_pathway_cluster5, NES > 0.0)
filt_up_cluster5$log10 <- -log10(filt_up_cluster5$pval)
filt_up_cluster5$pathway <- str_replace_all(filt_up_cluster5$pathway, c("HALLMARK"="", "_"=" "))



Cluster5_up_Pathways = ggplot(data=filt_up_cluster5, aes(x=reorder(pathway,log10), y=log10)) + 
  geom_bar(stat="identity", fill="red", width = 0.7) +
  coord_flip() +
  labs(x="", y="pvalue (-log10)") +
  theme_minimal() + ggtitle("Enriched pathways in Cluster 5 vs all")

Cluster5_up_Pathways

Cluster5_up_Pathways = ggsave("Top pathways for cluster 5.pdf",width = 8, height =4, units = "in", dpi = 300)



### collapse the leading edge gene list to vector
fgseaRes_Cluster5$leadingEdge <- vapply(fgseaRes_Cluster5$leadingEdge, paste, collapse = ",", character(1L))

stats2 <- data.frame(pathway=fgseaRes_Cluster5$pathway) %>%
  data.frame(pvalue=fgseaRes_Cluster5$pval) %>%
  data.frame(adjpvalue=fgseaRes_Cluster5$padj) %>%
  data.frame(NES=fgseaRes_Cluster5$NES) %>%
  data.frame(genes=fgseaRes_Cluster5$leadingEdge)



write.csv(stats2, file = "~/Desktop/DATA/scRNA seq/scRNAseq_Integration/New Analysis/Cluster 5 Pathways.csv", quote = TRUE, col.names = TRUE, row.names = FALSE) #quote=T to avoid splitting by ','




#####GSEA ANALYSIS_CLUSTER 2 vs 5
Cluster2vs5_GSEA = read.csv ("Marker genes distinct to clusters 2 vs 5.csv")
head(Cluster2vs5_GSEA)



###Organise and rank genes  by log2FC

Cluster2vs5_GSEA <- Cluster2vs5_GSEA %>% arrange(desc(avg_log2FC))
genelist <- Cluster2vs5_GSEA$avg_log2FC
names(genelist) = Cluster2vs5_GSEA$X
fold_changes_cluster2vs5 = genelist



###Running GSEA to determine genesets that are upregulated or downregualted


fgseaRes_Cluster2vs5 <- fgsea(pathways = hallmark, 
                              stats    = fold_changes_cluster2vs5,
                              minSize  = 15,
                              maxSize  = 500)

head(fgseaRes_Cluster2vs5[order(pval), ])

#### up/down regulated pathways are generated using the log10(pvalue).
filtered_pathway_cluster2vs5 <- subset(fgseaRes_Cluster2vs5, pval < 0.05)

filt_up_cluster2vs5 <- subset(filtered_pathway_cluster2vs5, NES > 0.0)
filt_up_cluster2vs5$log10 <- -log10(filtered_pathway_cluster2vs5$pval)
filt_up_cluster2vs5$pathway <- str_replace_all(filtered_pathway_cluster2vs5$pathway, c("HALLMARK"="", "_"=" "))



Cluster2vs5_up_Pathways = ggplot(data=filt_up_cluster2vs5, aes(x=reorder(pathway,log10), y=log10)) + 
  geom_bar(stat="identity", fill="red", width = 0.7) +
  coord_flip() +
  labs(x="", y="pvalue (-log10)") +
  theme_minimal() + ggtitle("Enriched pathways in Cluster 2 vs 5")

Cluster2vs5_up_Pathways

Cluster2vs5_up_Pathways = ggsave("Top pathways for cluster 2 vs 5.pdf",width = 8, height =4, units = "in", dpi = 300)



### collapse the leading edge gene list to vector
fgseaRes_Cluster2vs5$leadingEdge <- vapply(fgseaRes_Cluster2vs5$leadingEdge, paste, collapse = ",", character(1L))

stats3 <- data.frame(pathway=fgseaRes_Cluster2vs5$pathway) %>%
  data.frame(pvalue=fgseaRes_Cluster2vs5$pval) %>%
  data.frame(adjpvalue=fgseaRes_Cluster2vs5$padj) %>%
  data.frame(NES=fgseaRes_Cluster2vs5$NES) %>%
  data.frame(genes=fgseaRes_Cluster2vs5$leadingEdge)


write.csv(stats3, file = "~/Desktop/DATA/scRNA seq/scRNAseq_Integration/Cluster 2 vs 5 Pathways.csv", quote = TRUE, col.names = TRUE, row.names = FALSE) #quote=T to avoid splitting by ','



