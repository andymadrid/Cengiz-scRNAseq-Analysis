# This script gets the proportions of cells in each of the different clusters, amongst the different treatment groups

# Get cell numbers for cell-type proportions analysis
cellCounts <- table(Idents(seuratMerged.integrated),seuratMerged.integrated$Sample)
totalCells <- colSums(cellCounts)
for (i in 1:length(totalCells)) {
  cellProp <- c()
  for (j in 1:nrow(cellCounts)) {
    prop <- as.data.frame(cellCounts[j,i]/totalCells[i])
    cellProp <- rbind(cellProp,prop)
  }
  if (i == 1) {
    totalCellProps <- cellProp
  }
  else {
    totalCellProps <- cbind(totalCellProps,cellProp)
}}
colnames(totalCellProps) <- colnames(cellCounts)
rownames(totalCellProps) <- 0:(nrow(totalCellProps)-1)
totalCellProps <- totalCellProps*100
propMeans <- c()
for (i in 1:nrow(totalCellProps)) {
  f_hic <- as.data.frame(mean(as.numeric(totalCellProps[i,1:4])))
  f_hit <- as.data.frame(mean(as.numeric(totalCellProps[i,5])))
  f_sham <- as.data.frame(mean(as.numeric(totalCellProps[i,6:8])))
  m_hic <- as.data.frame(mean(as.numeric(totalCellProps[i,9:12])))
  m_sham <- as.data.frame(mean(as.numeric(totalCellProps[i,13:15])))
  colnames(f_sham) <- "F_Sham"
  colnames(m_sham) <- "M_Sham"
  colnames(f_hic) <- "F_HiC"
  colnames(f_hit) <- "F_HiT"
  colnames(m_hic) <- "M_HiC"
  f_hic_sem <- as.data.frame(mean(as.numeric(sd(totalCellProps[i,1:4]))/sqrt(4)))
  f_hit_sem <- as.data.frame(mean(as.numeric(sd(totalCellProps[i,5]))))
  f_sham_sem <- as.data.frame(mean(as.numeric(sd(totalCellProps[i,6:8]))/sqrt(3)))
  m_hic_sem <- as.data.frame(mean(as.numeric(sd(totalCellProps[i,9:12]))/sqrt(4)))
  m_sham_sem <- as.data.frame(mean(as.numeric(sd(totalCellProps[i,13:15]))/sqrt(3)))
  colnames(f_sham_sem) <- "F_Sham_SEM"
  colnames(m_sham_sem) <- "M_Sham_SEM"
  colnames(f_hic_sem) <- "F_HiC_SEM"
  colnames(f_hit_sem) <- "F_HiT_SEM"
  colnames(m_hic_sem) <- "M_HiC_SEM"
  clustSEM <- cbind(f_sham_sem,f_hic_sem,f_hit_sem,m_sham_sem,m_hic_sem)
  clusMean <- cbind(f_sham,f_hic,f_hit,m_sham,m_hic,clustSEM)
  propMeans <- rbind(propMeans,clusMean)
}
rownames(propMeans) <- 0:41
# make bar plots for each cluster
plot_list = list()
for (i in 1:nrow(propMeans)) {
  data <- t(as.data.frame(propMeans[i,]))
  data <- as.data.frame(data)
  data2 <- data
  data <- as.data.frame(data[1:5,])
  data$SEM <- data2[6:10,1] 
  data$Group <- c("F_Sham","F_HiC","F_HiT","M_Sham","M_HiC")
  colnames(data) <- c("Value","SEM","Group")
  data$Group <- factor(data$Group,levels=data$Group)
  tcProps <- as.data.frame(t(as.data.frame(totalCellProps[i,c(6:8,1:4,5,13:15,9:12)])))
  tcProps$Group <- rep(c("F_Sham","F_HiC","F_HiT","M_Sham","M_HiC"),c(3,4,1,3,4))
  colnames(tcProps) <- c("Points","Group")
  ii <- i-1
  mainLab <- paste0("Cluster ",ii)
  p <- ggplot(data=data,aes(x=Group,y=Value,fill=Group)) +  geom_bar(stat="identity", width=0.5,fill=c("firebrick3","chocolate2","goldenrod3","forestgreen","mediumblue"))   + theme_classic() + theme(legend.position="none") + xlab("") + ylab("Percent of Cells in Cluster") + theme(text=element_text(size=20)) + geom_point(data=tcProps,aes(x=Group,y=Points)) + geom_errorbar( aes(x=Group, ymin=Value-SEM, ymax=Value+SEM), width=0.4, colour="black", alpha=0.9, size=1.3) + ggtitle(mainLab) + theme(axis.text.x = element_text(angle = 45,hjust=1))
  plot_list[[i]] <- p
}
pdf("clusterPercentages.pdf",onefile=TRUE)
for (i in 1:length(plot_list)) {
 print(plot_list[[i]])
}
dev.off()
