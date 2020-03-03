


d$samples$group <- as.factor(group)
d$samples$batch <- as.factor(batch)

design <- model.matrix(~0+group+batch)
colnames(design) <- gsub("group", "", colnames(design))
design


#x <- c("CNS-Breast","Colon-Breast","Leukemia-Breast","Melanoma-Breast")
#contr.matrix <- makeContrasts(contrasts=x,levels=design)


design.pairs <- function(levels) { 
n <- length(levels) 
design <- matrix(0,n,choose(n,2)) 
rownames(design) <- levels 
colnames(design) <- 1:choose(n,2) 
k <- 0 
for (i in 1:(n-1)) 
for (j in (i+1):n) { 
k <- k+1 
design[i,k] <- 1 
design[j,k] <- -1 
colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="") 
} 
design 
}


contr.matrix <- design.pairs(unique(group))
tmp <- matrix(0, nrow=2, ncol=dim(contr.matrix)[2])
rownames(tmp) <- unique(batch)[1:2]
contr.matrix <- rbind(contr.matrix, tmp)


v <- voom(d, design, plot=TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

topTable(efit, number=10, adjust.method="BH", sort.by="B", p.value=1, lfc=0)
topTableF(efit, number=10, adjust.method="BH", sort.by="F", p.value=1, lfc=0)

write.fit(efit, results=NULL, file="results.txt", digits=4, adjust="BH", method="separate",F.adjust="BH", sep="\t")
