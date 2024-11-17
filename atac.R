
if(FALSE){
  BiocManager::install("chromVAR")
  devtools::install_github("GreenleafLab/chromVARmotifs")
}


library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(Signac)
library(Seurat)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(hdf5r)
library(plotly)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(umap)
library(devtools)
library(reticulate)
library(cowplot)
library(sceasy)
library(sqldf)
library(reshape2)
library(motifmatchr)
library(chromVARmotifs)


list_bsgenome <- list()
list_bsgenome[["mouse"]] <- BSgenome.Mmusculus.UCSC.mm10
list_annotation <- list()
list_annotation[["mouse"]] <- EnsDb.Mmusculus.v79
list_blacklist <- list()
list_blacklist[["mouse"]] <- blacklist_mm10

###### Genes of interest
gene_il4 <-  "ENSMUSG00000000869"
gene_il13 <- "ENSMUSG00000020383"
gene_il5 <-  "ENSMUSG00000036117"





################################################################################
###################### Genes in Th2 CRISPR screen ##############################
################################################################################

if(FALSE){
  crisprdata <- read.csv("th2crispr_mageck_output.csv")
  #cr2_Il4 cr2_Il13 cr2_Irf4 cr2_Xbp1 cr2_Gata3
  
  colnames(crisprdata)[1] <- "ensid"
  
  map_ensid_symbol <- read.csv(file.path('map_uniprot_mousesymbol.txt'), header = T, sep = '\t')[,c(1,3)]
  colnames(map_ensid_symbol) <- c('ensid','symbol')
  map_ensid_symbol <- unique(map_ensid_symbol)
  
  crisprdata <- merge(map_ensid_symbol, crisprdata)
  
  cr_il4 <- crisprdata$symbol[crisprdata$cr2_Il4<1000]
  cr_il13 <- crisprdata$symbol[crisprdata$cr2_Il13<1000]
  cr_gata3 <- crisprdata$symbol[crisprdata$cr2_Gata3<1000]
}




################################################################################
###################### load RNA ################################################
################################################################################


rna_adata <- readRDS("/husky/otherdataset/nicole_pth2_atac/rna/scRNAseq_Data_NEW.rds")


######## Figure out which genes are expressed
explevel_per_gene <- data.frame(
  symbol=rownames(rna_adata),
  mean_exp=rowMeans(rna_adata@assays$RNA$counts)
)
plot(sort(log10(1+explevel_per_gene$mean_exp)))


explevel_per_gene[explevel_per_gene$symbol %in% c("Gata3","Stat1","Stat2","Stat3","Stat5a","Irf4"),]
#0.05 is a sensible cutoff; 0.02

sum(explevel_per_gene$mean_exp>0.05)
expressed_genes <- explevel_per_gene[explevel_per_gene$mean_exp>0.01,]


rna_adata$is_ko <- str_detect(rna_adata$HTO_classification,"HDAC")
rna_adata$is_il13 <- str_detect(rna_adata$HTO_classification,"IL13plus")
rna_adata$is_naive <- str_detect(rna_adata$HTO_classification,"naive")
rna_adata$is_pathogenic <- rna_adata$seurat_clusters=="2"      #peTh2

if(FALSE){
  DimPlot(rna_adata, group.by = "HTO_classification")
  DimPlot(rna_adata, group.by = "is_ko")
  DimPlot(rna_adata, group.by = "is_il13")
  DimPlot(rna_adata, group.by = "is_naive")
  DimPlot(rna_adata) #cluster 2 is pTh2
  rna_adata$clusters <- Idents(rna_adata)
  
  
  
  list_favo <- c(
    "Il4",
    "Il5",
    "Il13",
    "Bhlhe40", "Zeb2", "Pparg",
    "Fosl2",
    "Junb",
    "Mbnl1",
    "Jarid2"
  )
  
  FeaturePlot(rna_adata, c(
    "Mbnl1"
  )) +scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 2.5
  )
  
  
  
  FeaturePlot(rna_adata, c(
    "Il4","Il13","Il5"
  )) +scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 2.5
  )
  
  VlnPlot(rna_adata, "Mbnl1", group.by = "HTO_classification")
  
  
}



################################################################################
###################### Volcano plot ############################################
################################################################################



list_cyto <- c("Il4", "Il13", "Il5")
list_tf <- read.csv("animaltfdb4.csv",sep="\t")[,"Symbol"]

unique(rna_adata$HTO_classification)


de_ko <- FindMarkers(rna_adata[,!rna_adata$is_naive], group.by = "is_pathogenic", ident.1 = "TRUE", ident.2 = "FALSE")
de_ko$gene <- rownames(de_ko)
ggplot(de_ko, aes(avg_log2FC, -log10(p_val), label=gene)) + 
  geom_point(color="gray") + 
  geom_text(data = de_ko[de_ko$p_val<1e-50 & de_ko$gene %in% c(list_tf),], color="red") +   ########## show this plot!!
  geom_text(data = de_ko[de_ko$gene %in% c(list_cyto),], color="blue") +   
  xlab("Log2(FC), pTh2 vs other non-naive Th2") + 
  ylab("-Log10(p)")
ggsave("new/new_volcano.svg", width = 8, height = 6)



rna_adata$is_5 <- rna_adata$seurat_clusters=="5"
rna_adata$is_25 <- rna_adata$seurat_clusters %in% c("2","5")

de_ko <- FindMarkers(rna_adata[,!rna_adata$is_naive], group.by = "is_5", ident.1 = "TRUE", ident.2 = "FALSE")
de_ko$gene <- rownames(de_ko)
ggplot(de_ko, aes(avg_log2FC, -log10(p_val), label=gene)) + 
  geom_point(color="gray") + 
  geom_text(data = de_ko[de_ko$p_val<1e-50 & de_ko$gene %in% c(list_tf),], color="red") +   
  geom_text(data = de_ko[de_ko$gene %in% c(list_cyto),], color="blue") +   
  xlab("Log2(FC), cl5 vs other non-naive Th2") + 
  ylab("-Log10(p)")
ggsave("new/volcano5.svg", width = 8, height = 6)


de_ko <- FindMarkers(rna_adata[,!rna_adata$is_naive], group.by = "is_25", ident.1 = "TRUE", ident.2 = "FALSE")
de_ko$gene <- rownames(de_ko)
ggplot(de_ko, aes(avg_log2FC, -log10(p_val), label=gene)) + 
  geom_point(color="gray") + 
  geom_text(data = de_ko[de_ko$p_val<1e-50 & de_ko$gene %in% c(list_tf),], color="red") +  
  geom_text(data = de_ko[de_ko$gene %in% c(list_cyto),], color="blue") +   
  xlab("Log2(FC), cl2+5 vs other non-naive Th2") + 
  ylab("-Log10(p)")
ggsave("new/volcano25.svg", width = 8, height = 6)





################################################################################
###################### Get gene-jaspar mapping #################################
################################################################################

### Load mapping motif-uniprot
map_jaspar_uniprot <- read.csv(file.path('jaspar_uniprot.csv'), header = T, sep = ',')[c(3,5,6)]
colnames(map_jaspar_uniprot) <- c("jaspar_baseid","jasparname","uniprot")
map_jaspar_uniprot <- map_jaspar_uniprot[map_jaspar_uniprot$uniprot!="",]

### Get uniprot to mouse mapping
map_uniprot_symbol <- read.csv(file.path('map_uniprot_mousesymbol.txt'), header = T, sep = '\t')[,c(2,3)]
colnames(map_uniprot_symbol) <- c('uniprot','symbol')
map_uniprot_symbol <- map_uniprot_symbol[map_uniprot_symbol$uniprot!="",]

### Get uniprot to human mapping
map_uniprot_humansymbol <- read.csv(file.path('/corgi/websites/tcellnet/shiny/data/uniprot_symbol_ensid.csv'), header = T, sep = '\t')[,c(2,3)]
colnames(map_uniprot_humansymbol) <- c('uniprot','human_symbol')
map_uniprot_humansymbol <- map_uniprot_humansymbol[map_uniprot_humansymbol$uniprot!="",]

### Use human-mouse orthology
map_human_mouse <- read.csv(file.path('map_mouse_human_symbol.csv'), header = T, sep = '\t')#[,c(2,3)]
colnames(map_human_mouse) <- c("human_symbol","symbol")
map_human_mouse <- map_human_mouse[map_human_mouse$human_symbol!="" & map_human_mouse$symbol!="",]

map_uniprot_symbol2 <- merge(map_uniprot_humansymbol, map_human_mouse)[,c("uniprot","symbol")]

## Only keep genes which are expressed
map_uniprot_symbol_all <- unique(rbind(map_uniprot_symbol, map_uniprot_symbol2))
map_uniprot_symbol_all <- merge(explevel_per_gene,map_uniprot_symbol_all)
#map_uniprot_symbol_all <- map_uniprot_symbol_all[map_uniprot_symbol$symbol %in% expressed_genes,]
map_jaspar_genesym <- unique(merge(map_uniprot_symbol_all, map_jaspar_uniprot))

unique(sort(map_jaspar_genesym$symbol)) #601 genes if no filtering



################################################################################
###################### load ATAC ###############################################
################################################################################


##########################
## Helper function to load from atacseq fragments
loadAtacFromFragments <- function(fragpath, organism){
  
  organism<- "mouse"
  fragpath <- "/husky/otherdataset/nicole_pth2_atac/for_signac/chr11.atac_fragments.tsv.gz"
  
  
  #### Index if needed
  fragpath_index <- paste(fragpath,".tbi",sep="")
  if(!file.exists(fragpath_index)){
    print("Indexing fragment file")
    system(paste("tabix -p vcf ",fragpath))
    #"tabix -p vcf fragments.tsv.gz"
    #"fragments.tsv.gz.tbi"
  }
  
  #### Count fragments
  counts <- CountFragments(fragments = fragpath)
  colnames(counts)[1] <- "cb"
  counts <- counts[order(counts$reads_count, decreasing = TRUE),]
  
  #knee plot analysis
  plot(sort(log10(counts$reads_count), decreasing = TRUE))
  #keep_cells <- counts$cb[log10(counts$reads_count)>4]
  keep_cells <- counts$cb[log10(counts$reads_count)>0]
  
  
  #### Create a dummy assay
  stupidmat <- matrix(0, nrow = 2, ncol=length(keep_cells))
  rownames(stupidmat) <- c("chr1:1-100", "chr1:200-300")
  colnames(stupidmat) <- keep_cells
  
  chrom_assay <- CreateChromatinAssay(
    counts = as.sparse(stupidmat),
    sep = c(":", "-"),
    fragments = fragpath,
    min.cells = 0,
    min.features = 0
  )
  
  adata <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC"
  )
  
  # call peaks using MACS2
  peaks <- CallPeaks(adata)
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
  
  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(adata),
    features = peaks,
    cells = colnames(adata)
  )
  the_annotation <- GetGRangesFromEnsDb(ensdb = list_annotation[[organism]])
  
  #seqnames here are 3... MT; no chr in front. can we trust alignment?
  seqlevels(the_annotation) <- paste0("chr",seqlevels(the_annotation)) #add chr to each chromosome name
  
  # create a new assay using the MACS2 peak set and add it to the Seurat object. Overwrite previous dummy set
  adata[["ATAC"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = fragpath,
    annotation =  the_annotation
  )
  
  
  adata
}



## Load the ATAC fragments
adata <- loadAtacFromFragments("/husky/otherdataset/nicole_pth2_atac/for_signac/chr11.atac_fragments.tsv.gz")

## Pull out metadata from file names
adata$genotype <- str_split_fixed(colnames(adata),"_",5)[,1]
adata$ct <- str_split_fixed(colnames(adata),"_",5)[,2]
adata$time <- str_split_fixed(colnames(adata),"_",5)[,3]
adata$replicate <- paste0("#",str_sub(str_split_fixed(colnames(adata),"R00",2)[,2],1,1))

##Subset to 48h time point
#adata48 <- adata[,adata$time=="48h"]


################################################################################
###################### dimensional reduction  ##################################  
################################################################################


## Dimensional reduction
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 3, verbose = T )
adata <- RunSVD(adata, verbose = T, reduction.name = 'SVD_ATAC')  #, n=ncol(adata)) #irlba.work=500, n=)

DepthCor(adata, reduction = 'SVD_ATAC')#,n = 20)
use_dim <- 2:5 # Could consider removing one more axis

adata <- RunUMAP(object = adata, reduction = 'SVD_ATAC', dims = use_dim, reduction.name = 'UMAP_ATAC', n.neighbors = ncol(adata))
adata <- FindNeighbors(object = adata, reduction = 'SVD_ATAC', dims = use_dim)

### Plot metadata on umap
adata$genotype_ct_time <- paste(adata$genotype, adata$ct, adata$time)

adata48 <- adata[,adata$time=="48h"]


p1 <- DimPlot(object = adata48, reduction = 'UMAP_ATAC', group.by = "genotype", pt.size = 2) 
p2 <- DimPlot(object = adata48, reduction = 'UMAP_ATAC', group.by = "ct", pt.size = 2)
#p3 <- DimPlot(object = adata48, reduction = 'UMAP_ATAC', group.by = "time", pt.size = 2)
p3 <- DimPlot(object = adata48, reduction = 'UMAP_ATAC', group.by = "genotype_ct_time", pt.size = 2)
ptot <- egg::ggarrange(p3,p2,p1, nrow=1)
ptot

ggsave(plot = ptot, "new/new_umap.svg", width = 10, height = 3)





p1 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "genotype", pt.size = 2) 
p2 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "ct", pt.size = 2)
p3 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "genotype_ct_time", pt.size = 2)
p4 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "time", pt.size = 2)
ptot <- egg::ggarrange(p3,p2,p1,p4, nrow=1)
ptot

#ggsave(plot = ptot, "new/new_umap.svg", width = 10, height = 3)

################################################################################
###################### inferred gene expression ################################
################################################################################

if(FALSE){
  gene.activities <- GeneActivity(adata)
  
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  adata[['RNA']] <- CreateAssayObject(counts = gene.activities)
  adata <- NormalizeData(
    object = adata,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(adata$nCount_RNA)
  )
  
  ## Plot our genes of interest on UMAP
  DefaultAssay(adata) <- "RNA"
  FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c(
    "Il5","Il13","Il4"
  ), ncol=3, pt.size = 5)
  
  
  ####### Produce gene expression heatmap 
  get_long_rna_matrix <- function(adata){
    
    #Get expression, normalize per gene
    tonorm <- as.matrix(adata@assays$RNA$data)
    tonorm <- t(scale(t(tonorm), center = FALSE))
    
    
    #Add metadata
    the_meta <- adata@meta.data
    the_meta$sample <- rownames(the_meta)
    long_rna <- melt(tonorm)
    colnames(long_rna) <- c("symbol","sample","value")
    long_rna <- merge(long_rna, the_meta)
    long_rna
  }
  
  long_rna <- get_long_rna_matrix(adata)
  
  
  
  
  #### Show our genes of interest in heatmap
  ggplot(long_rna[long_rna$symbol %in% c(
    "Il4","Il13","Il5"
  ),], aes(symbol, paste(time,ct,genotype,replicate), fill=value)) + 
    geom_tile() + 
    ylab("") + xlab("") + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 1
    )
  ggsave("heatmap_inferred_exp.svg", width = 4, height = 7)
  
}


################################################################################
###################### Load motifs #############################################
################################################################################


## Subset to only use standard chromosomes
DefaultAssay(adata) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- which(as.character(seqnames(granges(adata))) %in% main.chroms)
adata[["ATAC"]] <- subset(adata[["ATAC"]], features = rownames(adata[["ATAC"]])[keep.peaks])

## Get and add motifs from JASPAR
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
adata <- AddMotifs(adata, genome = BSgenome.Mmusculus.UCSC.mm10,pfm = pfm,verbose = T)


## Load motif data
motif2names <- data.frame(
  row.names=names(pfm),
  id=names(pfm),
  name=sapply(1:length(pfm), function(x) name(pfm[[x]]))
)
names2motifs <- motif2names
rownames(names2motifs) <- names2motifs$name




################################################################################
###################### ChromVAR analysis #######################################
################################################################################


library(BiocParallel)
register(MulticoreParam(2)) #avoid swap trashing by setting it somewhat low

DefaultAssay(adata) <- "ATAC"

## Compute motif activity; this is deviations in counts for each motif across peaks
## as compared to total count. it does not care /which/ peaks the motif is in, just the total sum
##https://www.nature.com/articles/nmeth.4401
adata <- RunChromVAR(
  object = adata,
  genome = BSgenome.Mmusculus.UCSC.mm10, assay = 'ATAC'
)




#############################
## Function to check correlation of motif activity with a variable
corr_chromvar <- function(adata, with_var, keep=NULL){
  DefaultAssay(adata) <- "chromvar"
  if(!is.null(keep)){
    adata <- adata[,keep]
    with_var <- with_var[keep]
  }
  thecor <- as.data.frame(Rfast::correls(
    with_var,
    t(as.matrix(adata@assays$chromvar@data))
  ))
  thecor$id <- rownames(thecor)
  thecor <- merge(names2motifs,thecor)
  thecor <- thecor[,c("name","id","correlation","p-value")]  
  colnames(thecor) <- c("jasparname","id","correlation","pvalue")
  thecor <- merge(map_jaspar_genesym, thecor)  #only expressed genes in map_jaspar_genesym
  unique(thecor)
}

topbottom_cor_chromvar <- function(thecor){
  thecor <- thecor[!is.na(thecor$correlation),]
  thecor <- thecor[thecor$mean_exp>0.01,]
  thecor <- unique(thecor[,c("id","jasparname","symbol","correlation","pvalue","mean_exp")])
  
  rbind(
    head(thecor[order(thecor$correlation),],n=20),
    tail(thecor[order(thecor$correlation),],n=20)
  )
}



######################### 
## Function to produce heatmap from chromvar
make_chromvar_heatmap <- function(adata){
  the_meta <- adata@meta.data
  the_meta$sample <- rownames(the_meta)
  long_activity <- melt(adata@assays$chromvar$data)
  colnames(long_activity) <- c("id","sample","value")
  long_activity <- merge(long_activity, the_meta)
  long_activity <- merge(long_activity, names2motifs)
  long_activity
}


############################## Look at different axes

#adata_orig <- adata
#adata <- adata[,adata$time=="48h"]
#adata <- adata_orig

adata48 <- adata[,adata$time=="48h"]
thecor <- corr_chromvar(adata48, adata48$ct=="pTh2")

w <- thecor
w <- unique(w[,colnames(w) %in% c("correlation","pvalue","jasparname")])
if(FALSE){
  ggplot(w, aes(correlation, -log10(pvalue), label=jasparname)) + 
    geom_point() + 
    ggrepel::geom_text_repel(max.overlaps = 50) +
    xlab("corr, pTh2 vs Th2")
}

w <- merge(w,map_jaspar_genesym)
w$is_de <- w$symbol %in% de_ko$gene[de_ko$p_val<1e-30]

w2 <- merge(w, data.frame(pval_de=de_ko$p_val, symbol=de_ko$gene))
w2$is_de <- w2$pval_de<1e-30

ggplot(w2, aes(correlation, -log10(pval_de), label=jasparname)) + 
  geom_point() + 
  ggrepel::geom_text_repel(data = w2[w2$is_de,], max.overlaps = 50) +
  xlab("Motif activity correlation, pTh2 vs Th2") + 
  ylab("-log10(p), Differentially expressed")
ggsave("new/chromvar vs volcano.svg",width = 10,height = 6)



################################################################################
###################### Peaks in the Th2 locus ##################################
################################################################################

region_peak_il5_1 <- "chr11-53713681-53713936"  #up in all pTh2 at 48h
region_peak_il5_2 <- "chr11-53693226-53693769"  #up in all pTh2 at 48h, for KO only
sum(rownames(adata)==region_peak_il5_2)


### Subset by coordinates. Picked by manual inspection in subset of Th2 locus
DefaultAssay(adata) <- "ATAC"
region_peak_il4_13_start <- 53589596  # 53592335
region_peak_il4_13_end   <- 53663727
isin_bigdip_region <- str_split_fixed(rownames(adata),"-",3)[,1]=="chr11" &
  as.integer(str_split_fixed(rownames(adata),"-",3)[,3])>=region_peak_il4_13_start &
  as.integer(str_split_fixed(rownames(adata),"-",3)[,2])<=region_peak_il4_13_end
sum(isin_bigdip_region)




### Subset by coordinates. Picked by manual inspection in the UCSC genome browser
DefaultAssay(adata) <- "ATAC"
isin_th2_locus <- str_split_fixed(rownames(adata),"-",3)[,1]=="chr11" &
  as.integer(str_split_fixed(rownames(adata),"-",3)[,3])>53505121 &
  as.integer(str_split_fixed(rownames(adata),"-",3)[,2])<53903660

#Number of peaks selected in region
sum(isin_th2_locus)


#isin_th2_locus <- isin_bigdip_region

#############################
# Function to extract count matrix for peaks in Th2 locus; long format for ggplot
get_atac_peaks_for_region <- function(adata, which_peaks){
  #Get values in long format
  DefaultAssay(adata) <- "ATAC"
  atac_count <- melt(as.matrix(adata@assays$ATAC$data[which_peaks, ]))
  colnames(atac_count) <- c("pos","sample","value")
  atac_count$start <- as.integer(str_split_fixed(atac_count$pos,"-",3)[,2])
  atac_count$pos <- factor(atac_count$pos,naturalsort::naturalsort(unique(atac_count$pos)))
  
  #Add metadata
  the_meta <- adata@meta.data
  the_meta$sample <- rownames(the_meta)
  atac_count <- merge(atac_count, the_meta)
  atac_count
}


th2_atac_count <- get_atac_peaks_for_region(adata, isin_th2_locus)

#### Plot the region as a heatmap, 
ggplot(th2_atac_count, aes(pos, paste(time,ct,genotype,replicate), fill=value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("") + 
  xlab("Coordinate on chr11") #+ coord_flip()
ggsave("new/new_heatmap_of_peaks.svg", width = 8, height = 8)  ### updated 666



if(FALSE){
  #### Plot the region in absolute coordinates, using dots
  ggplot(th2_atac_count, aes(start, paste(time,ct,genotype,replicate), color=value)) + 
    geom_point() + 
    ylab("") + 
    xlab("Coordinate on chr11")
  
  
  #Show the real coordinates and where our genes of interest are
  plot_map_bin_to_realcoord <- function(th2_atac_count){
    the_annotation <- adata@assays$ATAC@annotation
    p1 <- ggplot(th2_atac_count, aes(pos, start)) + geom_point() 
    p1 <- p1 + geom_hline(yintercept = min(start(the_annotation[the_annotation$gene_name=="Il4"])),color="red") + 
      geom_hline(yintercept = max(end(the_annotation[the_annotation$gene_name=="Il4"])),color="red")
    p1 <- p1 + geom_hline(yintercept = min(start(the_annotation[the_annotation$gene_name=="Il13"])),color="blue") + 
      geom_hline(yintercept = max(end(the_annotation[the_annotation$gene_name=="Il13"])),color="blue")
    p1 <- p1 + geom_hline(yintercept = min(start(the_annotation[the_annotation$gene_name=="Il5"])),color="green") + 
      geom_hline(yintercept = max(end(the_annotation[the_annotation$gene_name=="Il5"])),color="green")
    p1 + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      ylab("Start coordinate") + 
      xlab("Peak") 
  }
  
  plot_map_bin_to_realcoord(th2_atac_count)
  
  
}



################################################################################
###################### ChromVAR analysis, Th2 locus only #######################
################################################################################



DefaultAssay(adata) <- 'ATAC'


adata$is_wt_th2 <- adata$genotype_ct_time=="WT Th2 48h"

da_peaks <- FindMarkers(
  object = adata[,adata$time=="48h"],
  ident.1 = "TRUE",
  ident.2 = "FALSE",
  group.by = "is_wt_th2",
  test.use = 't',
  min.pct = 0.1
)  

str_th2_locus <- "chr11-53559303-53751444"

da_peaks$start <- start(ranges(StringToGRanges(rownames(da_peaks))))
da_peaks$end <- end(ranges(StringToGRanges(rownames(da_peaks))))
da_peaks <- da_peaks[da_peaks$start>53550000 & da_peaks$end<53800000,]
da_peaks <- da_peaks[order(da_peaks$start),]
da_peaks <- da_peaks[da_peaks$p_val<0.01,]
da_peaks


CoveragePlot(
  object = adata[,adata$time=="48h"],
  region = str_th2_locus,
  region.highlight = StringToGRanges(rownames(da_peaks)),
  extend.upstream = 0,
  extend.downstream = 0,
  group.by="genotype_ct_time"
) 
ggsave("new/covplot.svg", width = 15, height = 6)




