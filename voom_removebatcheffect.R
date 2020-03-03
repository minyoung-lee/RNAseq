setwd('E:/miRNA_NCI60/work_matlab')
library(edgeR)
rm(list=ls())


data <- as.matrix(read.table("readcount_cell_763miR.txt"))
group <- as.matrix(read.table("cancer.txt", sep="\t"))
batch <- as.matrix(read.table("batch_cell.txt", sep="\t"))
a=unique(batch)
bid=rep(1,dim(batch)[1])
for (i in 1:length(a)){
 bid[batch==a[i]]=i
}


libSizes <- as.vector(colSums(data))


d <- DGEList(counts=data,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
VST <- voom(d, design=NULL, plot=TRUE)
voom_matrix <- VST$E

voom_rb <- removeBatchEffect(voom_matrix, batch=batch)

tiff("plotMDS_voom_cell.tiff")
plotMDS(voom_matrix, labels=group, col=bid)
dev.off()

tiff("plotMDS_voom_removeBatchEffect_cell.tiff")
plotMDS(voom_rb, labels=group, col=bid)
dev.off()

fname="voom.norm.factors_cell.txt"
write.table(d$samples, file=fname, row.names = TRUE, col.names = TRUE, sep = "\t" )

fname="voom.logcpm_cell.txt"
write.table(voom_matrix, file=fname, row.names = TRUE, col.names = TRUE, sep = "\t" )

fname="voom.logcpm_removebatch_cell.txt"
write.table(voom_rb, file=fname, row.names = TRUE, col.names = TRUE, sep = "\t" )
