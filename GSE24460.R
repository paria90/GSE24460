library(Biobase)
library(limma)
library(GEOquery)
library(ggplot2)
library(pheatmap)
library(plyr)
library(reshape2)
library(gplots)

#1:Get data from GEO
series <- "GSE24460"
gset <- getGEO(series,GSEMatrix = TRUE,AnnotGPL = TRUE,destdir = "data/")
gset <- gset[[1]]
gr <- c(rep("MCF.P",2),rep("MCF",2))

#2:Expression matrix
ex <- exprs(gset)

#3:Quality control
#3-1:
pdf("results/boxplpot.pdf",width = 4)           #boxplot(ex)
boxplot(ex)
dev.off()

#3-2:
pdf("results/corheatmap.pdf")                   #pheatmap(cor(ex),labels_row = gr,labels_col = gr)
pheatmap(cor(ex),labels_row = gr,labels_col = gr)
dev.off()

#3-3:
ex.scale <- t(scale(t(ex),scale = FALSE))
pc <- prcomp(ex.scale)
pdf("results/pc-scale.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()
 
#3-4:
pcr <- data.frame(pc$rotation[,1:3],Group=gr)
pdf("results/pca-samples.pdf")
ggplot(pcr,aes(PC1,PC2,color=Group)) +geom_point(size=3) +theme_bw()
dev.off()

#4:
gr <- as.factor(gr)
gset$description <- gr
design <- model.matrix(~description+0,gset)
colnames(design) <- levels(gr)
fit <- lmFit(gset,design)
cont.matrix <- makeContrasts(MCF.P-MCF,levels = design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2,0.01)
tT <- topTable(fit2,adjust.method = "fdr",number = Inf)
tT <- subset(tT,select=c("Gene.symbol","Gene.ID","logFC","adj.P.Val"))
write.table(tT,"results/MCF.P-MCF.txt",row.names = FALSE,sep = "\t",quote = FALSE)

#5:
MCF.up <- subset(tT,logFC > 1 & adj.P.Val < 0.05)
MCF.up.genes <- unique(as.character(strsplit2(MCF.up$Gene.symbol,"///")))
write.table(MCF.up.genes,file = "results/MCF-up.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

MCF.down <- subset(tT,logFC < -1 & adj.P.Val < 0.05)
MCF.down.genes <- unique(as.character(strsplit2(MCF.down$Gene.symbol,"///")))
write.table(MCF.down.genes,file = "results/MCF-down.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
