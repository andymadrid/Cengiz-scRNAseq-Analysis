# The following script reads in the preprocessed samples and performs downstream analyses
# This includes finding variable genes, normalizing and scaling samples, clustering nuclei, and more...

# load packages
library(Seurat)
library(dplyr)
library(SingleR)
library(remotes)
library(Matrix)
library(ggplot2)
library(scRNAseq)
library(scuttle)
library(BiocParallel)
library(tibble)
library(DESeq2)

# Load in the preprocessed data
fh <- as.list(dir(pattern="rdata"))
lapply(fh,load,.GlobalEnv)

# Add some sample information into the seurat objects
seurat_M_Sham_1$Sex <- "Male"
seurat_M_Sham_2$Sex <- "Male"
seurat_M_Sham_3$Sex <- "Male"
seurat_M_HiC_1$Sex <- "Male"
seurat_M_HiC_2$Sex <- "Male"
seurat_M_HiC_3$Sex <- "Male"
seurat_M_HiC_4$Sex <- "Male"
seurat_F_Sham_1$Sex <- "Female"
seurat_F_Sham_2$Sex <- "Female"
seurat_F_Sham_3$Sex <- "Female"
seurat_F_HiC_1$Sex <- "Female"
seurat_F_HiC_2$Sex <- "Female"
seurat_F_HiC_3$Sex <- "Female"
seurat_F_HiC_4$Sex <- "Female"
seurat_F_HiT_1$Sex <- "Female"

seurat_M_Sham_1$Treatment <- "Sham"
seurat_M_Sham_2$Treatment <- "Sham"
seurat_M_Sham_3$Treatment <- "Sham"
seurat_M_HiC_1$Treatment <- "HiC"
seurat_M_HiC_2$Treatment <- "HiC"
seurat_M_HiC_3$Treatment <- "HiC"
seurat_M_HiC_4$Treatment <- "HiC"
seurat_F_Sham_1$Treatment <- "Sham"
seurat_F_Sham_2$Treatment <- "Sham"
seurat_F_Sham_3$Treatment <- "Sham"
seurat_F_HiC_1$Treatment <- "HiC"
seurat_F_HiC_2$Treatment <- "HiC"
seurat_F_HiC_3$Treatment <- "HiC"
seurat_F_HiC_4$Treatment <- "HiC"
seurat_F_HiT_1$Treatment <- "HiT"

seurat_M_Sham_1$Group <- "M_Sham"
seurat_M_Sham_2$Group <- "M_Sham"
seurat_M_Sham_3$Group <- "M_Sham"
seurat_M_HiC_1$Group <- "M_HiC"
seurat_M_HiC_2$Group <- "M_HiC"
seurat_M_HiC_3$Group <- "M_HiC"
seurat_M_HiC_4$Group <- "M_HiC"
seurat_F_Sham_1$Group <- "F_Sham"
seurat_F_Sham_2$Group <- "F_Sham"
seurat_F_Sham_3$Group <- "F_Sham"
seurat_F_HiC_1$Group <- "F_HiC"
seurat_F_HiC_2$Group <- "F_HiC"
seurat_F_HiC_3$Group <- "F_HiC"
seurat_F_HiC_4$Group <- "F_HiC"
seurat_F_HiT_1$Group <- "F_HiT"

seurat_M_Sham_1$Sample <- "M_Sham1"
seurat_M_Sham_2$Sample <- "M_Sham2"
seurat_M_Sham_3$Sample <- "M_Sham3"
seurat_M_HiC_1$Sample <- "M_HiC1"
seurat_M_HiC_2$Sample <- "M_HiC2"
seurat_M_HiC_3$Sample <- "M_HiC3"
seurat_M_HiC_4$Sample <- "M_HiC4"
seurat_F_Sham_1$Sample <- "F_Sham1"
seurat_F_Sham_2$Sample <- "F_Sham2"
seurat_F_Sham_3$Sample <- "F_Sham3"
seurat_F_HiC_1$Sample <- "F_HiC1"
seurat_F_HiC_2$Sample <- "F_HiC2"
seurat_F_HiC_3$Sample <- "F_HiC3"
seurat_F_HiC_4$Sample <- "F_HiC4"
seurat_F_HiT_1$Sample <- "F_HiT1"

# Change directory for plots and such
setwd("/media/Data/RNAseq/scRNA/Pelin_15_Mouse_Hip_July2022/15_Samples_Plots")

# Merge the objects

seuratMerged <- merge(x=seurat_M_Sham_1,y=list(seurat_M_Sham_2,seurat_M_Sham_3,seurat_F_Sham_1,seurat_F_Sham_2,seurat_F_Sham_3,seurat_M_HiC_1,seurat_M_HiC_2,seurat_M_HiC_3,seurat_M_HiC_4,seurat_F_HiC_1,seurat_F_HiC_2,seurat_F_HiC_3,seurat_F_HiC_4,seurat_F_HiT_1),add.cell.ids=c("M_Sham_1","M_Sham_2","M_Sham_3","F_Sham_1","F_Sham_2","F_Sham_3","M_HiC_1","M_HiC_2","M_HiC_3","M_HiC_4","F_HiC_1","F_HiC_2","F_HiC_3","F_HiC_4","F_HiT_1"))

# Remove objects to clean up
rm(seurat_F_Sham_1,seurat_F_Sham_2,seurat_F_Sham_3,seurat_F_HiC_1,seurat_F_HiC_2,seurat_F_HiC_3,seurat_F_HiC_4,seurat_F_HiT_1,seurat_M_Sham_1,seurat_M_Sham_2,seurat_M_Sham_3,seurat_M_HiC_1,seurat_M_HiC_2,seurat_M_HiC_3,seurat_M_HiC_4)

# Split the data
seuratMerged.list <- SplitObject(seuratMerged, split.by = "Sample")

# Normalize split data
seuratMerged.list <- lapply(X = seuratMerged.list, FUN = function(x) {
x <- NormalizeData(x, verbose = FALSE)
x <- FindVariableFeatures(x, verbose = FALSE)
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = seuratMerged.list)
#seuratMerged.list <- lapply(X = seuratMerged.list, FUN = function(x) {
#x <- ScaleData(x, features = features, verbose = FALSE)
#x <- RunPCA(x, features = features, verbose = FALSE)
#})

# Identify integration anchors and integrate data
anchors <- FindIntegrationAnchors(object.list = seuratMerged.list, dims = 1:50)
seuratMerged.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

# Scale and cluster data
DefaultAssay(seuratMerged.integrated) <- "integrated"
seuratMerged.integrated <- ScaleData(seuratMerged.integrated, verbose = FALSE)
seuratMerged.integrated <- RunPCA(seuratMerged.integrated, verbose = FALSE)
seuratMerged.integrated <- RunUMAP(seuratMerged.integrated, dims = 1:50)
seuratMerged.integrated <- FindNeighbors(seuratMerged.integrated, dims = 1:50)
seuratMerged.integrated <- FindClusters(seuratMerged.integrated, resolution = 0.6)
save(seuratMerged.integrated,file="seuratMerged.integrated.rdata")

# UMAP plot for data exploration
pdf("UMAP.pdf")
plot <- DimPlot(seuratMerged.integrated, reduction = "umap") + NoLegend()
LabelClusters(plot=plot,id="ident")
dev.off()

pdf("UMAP2.pdf")
plot <- DimPlot(seuratMerged.integrated, reduction = "umap")
LabelClusters(plot=plot,id="ident")
dev.off()

pdf("UMAP_PredictedCellType.pdf")
DimPlot(seuratMerged.integrated,dims=1:2,group.by="predIdent")
dev.off()

pdf("UMAP_Sex.pdf")
DimPlot(seuratMerged.integrated,dims=1:2,group.by="Sex")
dev.off()

pdf("UMAP_Treatment.pdf")
DimPlot(seuratMerged.integrated,dims=1:2,group.by="Treatment")
dev.off()

pdf("UMAP_Group.pdf")
DimPlot(seuratMerged.integrated,dims=1:2,group.by="Group")
dev.off()

pdf("UMAP_Group_Split.pdf",width=13)
DimPlot(seuratMerged.integrated,dims=1:2,split.by="Group")
dev.off()

goi <- c("Esr1","Esr2","Gper1","Ntrk2","Mtor","Bdnf","Mapk3","Mapk1","Akt1","Akt3","Gadd45b","Tet1","Tet2","Dnmt3a","Dnmt3b","Src","Fyn","Tnf","Casp3","Casp9","Parp1","Oxtr","Gsdme","Atp8b1","Itga3","Cd9","Fry","Nanos1","Crim1","Ntn4","Apln","Apela","Cyp19a1","Pgr","Trp53","Cav1","Creb1","Pik3ca","Il4","Il1b","Il10","Tgfb1")
goi1 <- goi[1:9]
goi2 <- goi[10:18]
goi3 <- goi[19:27]
goi4 <- goi[28:36]
goi5 <- goi[37:42]

# Make violin and features plots for the genes of interest
pdf("VlnPlot_GenesOfInterest1.pdf")
VlnPlot(seuratMerged,features=goi1) + theme(axis.text.x = element_text(angle = 90))
dev.off()
pdf("VlnPlot_GenesOfInterest2.pdf")
VlnPlot(seuratMerged,features=goi2) + theme(axis.text.x = element_text(angle = 90))
dev.off()
pdf("VlnPlot_GenesOfInterest3.pdf")
VlnPlot(seuratMerged,features=goi3) + theme(axis.text.x = element_text(angle = 90))
dev.off()
pdf("VlnPlot_GenesOfInterest4.pdf")
VlnPlot(seuratMerged,features=goi4) + theme(axis.text.x = element_text(angle = 90))
dev.off()
pdf("VlnPlot_GenesOfInterest5.pdf")
VlnPlot(seuratMerged,features=goi5) + theme(axis.text.x = element_text(angle = 90))
dev.off()
pdf("FeaturePlot_GenesOfInterest1.pdf")
FeaturePlot(seuratMerged,feature=goi1)
dev.off()
pdf("FeaturePlot_GenesOfInterest2.pdf")
FeaturePlot(seuratMerged,feature=goi2)
dev.off()
pdf("FeaturePlot_GenesOfInterest3.pdf")
FeaturePlot(seuratMerged,feature=goi3)
dev.off()
pdf("FeaturePlot_GenesOfInterest4.pdf")
FeaturePlot(seuratMerged,feature=goi4)
dev.off()
pdf("FeaturePlot_GenesOfInterest5.pdf")
FeaturePlot(seuratMerged,feature=goi5)
dev.off()

# Hierarchical Clustering of Clusters
library(ggtree)
seuratMerged.integrated <- BuildClusterTree(seuratMerged.integrated)
myPhyTree <- Tool(object=seuratMerged.integrated, slot = "BuildClusterTree")
pdf("/media/Data/tree.pdf")
ggtree(myPhyTree)+geom_tiplab()+theme_tree()+xlim(NA,400)
#PlotClusterTree(seuratMerged,direction="rightward")
dev.off()

# Dot plot using ordered clusters based on hierarchical clustering
seuratMerged.integrated$active.clust <- factor(seuratMerged.integrated$seurat_clusters,levels=rev(c(32,0,8,19,7,20,2,14,31,43,28,34,27,23,10,39,25,3,33,22,44,5,36,35,15,16,6,29,38,30,42,4,40,1,12,9,17,24,11,37,18,26,21,41,13)))

#seuratMerged.integrated$active.clust <- factor(Idents(seuratMerged.integrated),levels=rev(c("DG GC (Prox1)","CA1 Dorsal","CA1 Dorsal/Ventral","CA1 Ventral","CA2-4","EC L6 IT","EC L2","InN (Ssy Npy)","InN (Vip Nr2f2)","InN (Chol)" ,"InN (Lamp5 Lhx6)","Astro","Protoplasmic Astro","OPC","Oligo","Microglia","Macro","Immune (T cell)","Endo","Epend","VLMC")))

load("/media/Data/RNAseq/scRNA/Pelin_15_Mouse_Hip_July2022/andreCellMarkers.rdata")
goi <- andreCellMarkers

gaba.glut.genes <- c("Glul","Grin1","Grin2b","Slc17a7","Slc17a6","Gabbr1","Gabbr2","Gad1","Gad2","Slc6a1","Slc32a1")

#DefaultAssay(seuratMerged.integrated) <- "integrated"
DefaultAssay(seuratMerged.integrated) <- "RNA"
pdf("/media/Data/goi.pdf",width=20,height=12)
DotPlot(seuratMerged.integrated,features=goi,cols=c("lightgrey","red"),group.by="active.clust") + theme(axis.text.x = element_text(angle = 90))
dev.off()

# Identify cluster markers for annotation
DefaultAssay(seuratMerged.integrated) <- "RNA"
allMarkers <- FindAllMarkers(seuratMerged.integrated)

# Identify Conversed Markers for Cluster Annotation
DefaultAssay(seuratMerged.integrated) <- "RNA"
conservedMarkers <-c()
consMarker <- function(x) {
clusters <- levels(x$seurat_clusters)
for (i in 0:(length(clusters)-1)) {
cMarkers <- FindConservedMarkers(seuratMerged.integrated, ident.1 = i, grouping.var = "Group", verbose = FALSE)
cMarkers$cluster <- i
cMarkers <- rownames_to_column(cMarkers,"Gene")
conservedMarkers <- rbind(conservedMarkers,cMarkers)
prog <- paste0("Done with cluster ",i)
prog <- noquote(prog)}
return(conservedMarkers)}
conservedMarkers <- consMarker(seuratMerged.integrated)
save(conservedMarkers,file="conservedMarkers.rdata")

# Identify cluster markers
DefaultAssay(seuratMerged.integrated) <- "RNA"
clusterMarkers <- c()
clusters <- levels(seuratMerged.integrated$seurat_clusters)
for (i in 0:(length(clusters)-1)) {
prog <- paste0("Analyzing ",i,"\n")
prog <- noquote(prog)
cat(prog)
markers <- FindMarkers(seuratMerged.integrated,ident.1=i,only.pos=TRUE)
markers$cluster <- i
markers <- rownames_to_column(markers,"Gene")
clusterMarkers <- rbind(clusterMarkers,markers)
}
save(clusterMarkers,file="clusterMarkers.rdata")

