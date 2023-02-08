# The following script preprocesses each of the different samples
# Steps include doing a preliminary cell-type identification using previous data from another publication
# Then reading removing poor quality nuclei and nuclei with too much mitochondrial/ribosomal RNA

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

# create a reference from Zeisel brain data to predict cell types
ref <- ZeiselBrainData()
ref <- ref[,!is.na(ref$level2class)]
ref <- logNormCounts(ref)

# make a function to read in data,predict cell types using the Zeisel brain data reference, and filter cells
# in this case, x object is a data.frame with three columns: 
# 1) the directory where data is stored
# 2) the name of the Seurat object that will be save with predicted cell labels
# 3) the name of the file where the Seurat object will be saved as an .rdata file
# the print out is kind of clunky….but what’s a girl to do? Too lazy to optimize it…

prepareSeuratData <- function(x) {
  for (i in 1:nrow(x)) {
    prog <- i
    prog <- noquote(prog)
    cat(prog)
    prog <- "\nBeginning preprocessing of "
    prog <- noquote(prog)
    cat(prog)
    prog <- noquote(x[i,2])
    cat(prog)
    prog <- " now...\n"
    prog <- noquote(prog)
    cat(prog)
    prog <- "\tReading in the data now...\n"
    prog <- noquote(prog)
    cat(prog)
    data <- Read10X(data.dir=as.character(x[i,1]))
    outputFH <- x[i,2]
    outputFH <- CreateSeuratObject(counts=data,project=as.character(outputFH),min.cells=3,min.features=200)
    prog <- "\tPredicting cell types now...\n"
    prog <- noquote(prog)
    cat(prog)
    #predictedCells <- SingleR(test=data,ref=ref,assay.type.test=1,labels=ref$level1class,de.method="wilcox")
    predictedCells <- SingleR(test=data,ref=ref,assay.type.test=1,labels=ref$level2class,de.method="wilcox")
    outputFH$predIdent <- predictedCells$labels
    outputFH[["percent.mt"]] <- PercentageFeatureSet(outputFH,pattern="^mt-")
    outputFH[["percent.rib"]] <- PercentageFeatureSet(outputFH,pattern="^Rp[sl]")
    outputFH[["percent.cyto"]] <- outputFH[["percent.mt"]] + outputFH[["percent.rib"]]
    prog <- "\tFiltering cells and removing mito and ribo genes now...\n"
    prog <- noquote(prog)
    cat(prog)
    outputFH <- subset(outputFH, subset = nFeature_RNA > 500 & percent.cyto < 5)
    cts <- GetAssayData(outputFH, assay = "RNA")
    cts <- cts[-(grep(pattern="^mt-",x=rownames(cts))),]
    cts <- cts[-(grep(pattern="^Rp[sl]",x=rownames(cts))),]
    outputFH <- subset(outputFH,features=rownames(cts))
    eval(call("<-",as.name(x[i,2]),outputFH))
    prog <- "\tSaving data now...\n"
    prog <- noquote(prog)
    cat(prog)
    save(list=x[i,2],file=as.character(x[i,3]))
    prog <- "\tCleaning up because I’m a dirty little thing…\n"
    cat(prog)
    rm(outputFH,data,predictedCells,cts)
    prog <- "\tDone! Let’s have some fun!\n\n"
    prog <- noquote(prog)
    cat(prog)
  }
}

# Read in the file containing the information for data preprocessing
 metasheet <- read.csv("seuratInput.csv",header=T)

# Run the function
prepareSeuratData(metasheet)
rm(list=ls())

