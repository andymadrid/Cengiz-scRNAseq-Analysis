# This script performs differential gene expression analysis of the processed samples
# There are several different tpyes of differential analyses to be done, such as sex-specific changes, and inter-group changes

# This is a function to get differential expression for each cluster between the two groups specified
# x = the seuratMerged object
# y = one group of interest (e.g., “M_Sham”)
# z = another group of interest to compare to (e.g., “F_Sham”)

DefaultAssay(seuratMerged.integrated) <- "RNA"
clusDiff <- c()
getDiff <- function(x,y,z) {
  clusters <- levels(x$seurat_clusters)
  for (i in 0:(length(clusters)-1)) {
    isClus <- which(x$seurat_clusters==i & seuratMerged.integrated$Group==y)
    if (length(isClus) > 3) {
      isClus2 <- which(x$seurat_clusters==i & seuratMerged.integrated$Group==z)
      if (length(isClus2) > 3) {
        markers <- FindMarkers(seuratMerged.integrated,subset.ident=i,group.by="Group",ident.1=y,ident.2=z,only.pos=FALSE,logfc.threshold = 0.3)
        markers$cluster <- i
        markers <- rownames_to_column(markers,"Gene")
        #clusDiff <- append(clusDiff,markers)
        clusDiff <- rbind(clusDiff,markers)
        }}
    prog <- paste0("Done with cluster ",i)
    prog <- noquote(prog)
  }
  return(clusDiff)
}

res_M_F_Sham <- getDiff(seuratMerged.integrated,"M_Sham","F_Sham")
res_M_F_HiC <- getDiff(seuratMerged.integrated,"M_Sham","F_HiC")
res_M_Sham_HiC <- getDiff(seuratMerged.integrated,"M_Sham","M_HiC")
res_F_Sham_HiC <- getDiff(seuratMerged.integrated,"F_Sham","F_HiC")
res_F_Sham_HiT <- getDiff(seuratMerged.integrated,"F_Sham","F_HiT")
res_F_HiC_HiT <- getDiff(seuratMerged.integrated,"F_HiC","F_HiT")

res_M_F_Sham <- res_M_F_Sham[which(res_M_F_Sham$p_val_adj < 0.05),]
res_M_F_HiC <- res_M_F_HiC[which(res_M_F_HiC$p_val_adj < 0.05),]
res_M_Sham_HiC <- res_M_Sham_HiC[which(res_M_Sham_HiC$p_val_adj < 0.05),]
res_F_Sham_HiC <- res_F_Sham_HiC[which(res_F_Sham_HiC$p_val_adj < 0.05),]
res_F_Sham_HiT <- res_F_Sham_HiT[which(res_F_Sham_HiT$p_val_adj < 0.05),]
res_F_HiC_HiT <- res_F_HiC_HiT[which(res_F_HiC_HiT$p_val_adj < 0.05),]
save(res_M_F_Sham,res_M_F_HiC,res_M_Sham_HiC,res_F_Sham_HiC,res_F_Sham_HiT,res_F_HiC_HiT,file="differentialResults.RNAassay.rdata")
write.table(res_M_F_Sham,file="de.M.F.Sham.txt",quote=F,row.names=F,sep='\t')
write.table(res_M_F_HiC,file="de.M.F.HiC.txt",quote=F,row.names=F,sep='\t')
write.table(res_M_Sham_HiC,file="de.M.Sham.HiC.txt",quote=F,row.names=F,sep='\t')
write.table(res_F_Sham_HiC,file="de.F.Sham.HiC.txt",quote=F,row.names=F,sep='\t')
write.table(res_F_Sham_HiT,file="de.F.Sham.HiT.txt",quote=F,row.names=F,sep='\t')
write.table(res_F_HiC_HiT,file="de.F.HiC.HiT.txt",quote=F,row.names=F,sep='\t')


##############################
# Pseudo-bulk RNA-seq analysis
##############################

# Sex-Group interaction Analysis
DefaultAssay(seuratMerged.integrated) <- "RNA"
seuratMerged.integrated.diet <- seuratMerged.integrated[,which(seuratMerged.integrated$Group!="F_HiT")]
pseudoDiff <- c()
pseudoBulkDE <- function(x) {
  clusters <- levels(x$seurat_clusters)
  for (i in 0:(length(clusters)-1)) {
    prog <- paste0(i,"\n")
    prog <- noquote(prog)
    cat(prog)
    seuratSub <- x[,which(x$seurat_clusters==i)]
    aggCounts <- as.data.frame(AggregateExpression(seuratSub,assays="RNA",group.by="Sample",slot="counts"))
    # remove pseudogenes from analysis
    aggCounts <- aggCounts[!grepl("Gm",rownames(aggCounts)),]
    if (length(levels(factor(seuratSub$Sample)))==14) {
      colData <- as.data.frame(rep(c("F","M"),c(7,7)))
      colnames(colData) <- c("Sex")
      colData$Group <- rep(c("HiC","Sham","HiC","Sham"),c(4,3,4,3))
      colData$Group <- factor(colData$Group)
      colData$Sex <- factor(colData$Sex)
      colData$condition <- factor(paste0(colData$Group,"_",colData$Sex))
      dds <- DESeqDataSetFromMatrix(countData=aggCounts,colData=colData,design= ~ Group + Sex + Group:Sex)
      keep <- rowSums(counts(dds)) > 1
      dds <- dds[keep,]
      keep <- rowSums(counts(dds) >= 10) >= 8
      dds <- dds[keep,]
      dds <- DESeq(dds)
      resultsNames(dds)
      res <- results(dds,name="GroupSham.SexM")
      res <- res[order(res$pvalue),]
      #res$log2FoldChange <- -1*res$log2FoldChange
      res <- as.data.frame(res)
      res$cluster <- i
      res <- rownames_to_column(res,"Gene")
      res <- res[which(res$pvalue<0.05),]
      pseudoDiff <- rbind(pseudoDiff,res)
      }}
  return(pseudoDiff)
}
pseudoDiffRes <- pseudoBulkDE(seuratMerged.integrated.diet)
save(pseudoDiffRes,file="pseudoDiffRes.clsuters.rdata")


# Look at only female samples for DE across treatment groups
DefaultAssay(seuratMerged.integrated) <- "RNA"
seuratMerged.integrated.diet <- seuratMerged.integrated[,which(seuratMerged.integrated$Sex!="Male")]
pseudoFemaleDiff <- c()
pseudoFemaleBulkDE <- function(x) {
  clusters <- levels(x$seurat_clusters)
  for (i in 0:(length(clusters)-1)) {
    prog <- paste0(i,"\n")
    prog <- noquote(prog)
    cat(prog)
    seuratSub <- x[,which(x$seurat_clusters==i)]
    aggCounts <- as.data.frame(AggregateExpression(seuratSub,assays="RNA",group.by="Sample",slot="counts"))
    # remove pseudogenes from analysis
    if (ncol(aggCounts)==8) {
      aggCounts <- aggCounts[!grepl("Gm",rownames(aggCounts)),]
      aggCounts <- aggCounts[,c(6:8,1:4,5)]
      colData <- as.data.frame(factor(rep(c("Sham","HiC","HiT"),c(3,4,1))))
      colnames(colData) <- "condition"
      colData$Group <- as.numeric(rep(c(1,2,3),c(3,4,1)))
      dds <- DESeqDataSetFromMatrix(countData=aggCounts,colData=colData,design= ~ Group)
      keep <- rowSums(counts(dds)) > 1
      dds <- dds[keep,]
      keep <- rowSums(counts(dds) >= 10) >= 3
      dds <- dds[keep,]
      dds <- DESeq(dds)
      res <- results(dds,name="Group")
      res <- res[order(res$pvalue),]
      res <- as.data.frame(res)
      res$cluster <- i
      res <- rownames_to_column(res,"Gene")
      pseudoFemaleDiff <- rbind(pseudoFemaleDiff,res)
      }}
  return(pseudoFemaleDiff)
}
pseudoDiffRes_F_Sham_HiC_HiT <- pseudoFemaleBulkDE(seuratMerged.integrated.diet)
pseudoDiffRes_F_Sham_HiC_HiT <- pseudoDiffRes_F_Sham_HiC_HiT[which(pseudoDiffRes_F_Sham_HiC_HiT$pvalue<0.05),]

# Another (better?) way for three way group rescue analysis
# Look at only female samples for DE across treatment groups
DefaultAssay(seuratMerged.integrated) <- "RNA"
seuratMerged.integrated.diet <- seuratMerged.integrated[,which(seuratMerged.integrated$Sex!="Male")]
pseudoFemaleDiff <- c()
pseudoFemaleBulkDE <- function(x) {
  clusters <- levels(x$seurat_clusters)
  for (i in 0:(length(clusters)-1)) {
    prog <- paste0(i,"\n")
    prog <- noquote(prog)
    cat(prog)
    seuratSub <- x[,which(x$seurat_clusters==i)]
    aggCounts <- as.data.frame(AggregateExpression(seuratSub,assays="RNA",group.by="Sample",slot="counts"))
    # remove pseudogenes from analysis
    if (ncol(aggCounts)==8) {
      aggCounts <- aggCounts[!grepl("Gm",rownames(aggCounts)),]
      aggCounts <- aggCounts[,c(6:8,1:4,5)]
      colData <- as.data.frame(factor(rep(c("Sham","HiC","Sham"),c(3,4,1))))
      colnames(colData) <- "condition"
      dds <- DESeqDataSetFromMatrix(countData=aggCounts,colData=colData,design= ~ condition)
      keep <- rowSums(counts(dds)) > 1
      dds <- dds[keep,]
      keep <- rowSums(counts(dds) >= 10) >= 3
      dds <- dds[keep,]
      dds <- DESeq(dds)
      res <- results(dds,name="condition_Sham_vs_HiC")
      res <- res[order(res$pvalue),]
      res <- as.data.frame(res)
      res$cluster <- i
      res <- rownames_to_column(res,"Gene")
      pseudoFemaleDiff <- rbind(pseudoFemaleDiff,res)
      }}
  return(pseudoFemaleDiff)
}
pseudoDiffRes_F_Sham_HiC_HiT <- pseudoFemaleBulkDE(seuratMerged.integrated.diet)
pseudoDiffRes_F_Sham_HiC_HiT <- pseudoDiffRes_F_Sham_HiC_HiT[which(pseudoDiffRes_F_Sham_HiC_HiT$pvalue<0.05),]
