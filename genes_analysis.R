setwd("~/Desktop/University/term5/Bioinformatic/Bioinformatics_AML_Gene")
library(Biobase)
library(GEOquery)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(raster)
library(tidyverse)
library(palmerpenguins)
library(Rtsne)


Sys.setenv("VROOM_CONNECTION_SIZE" = 13107210)

series = "GSE48558"
platform = "GPL6244"

# load the data
gset = GEOquery::getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "Data/")
if (length(gset) > 1) idx = grep(platform, attr(gset, "names")) else idx <- 1
gset = gset[[idx]]

## group by source name
gr = c('AML_Patient', 'AML_Patient', 'AML_Patient', 'AML_Patient', 'AML_Patient', 'AML_Patient', 'AML_Patient', 
       'AML_Patient', 'AML_Patient', 'AML_Patient', 'AML_Patient', 'AML_Patient', 'AML_Patient', 'B_ALL_Cell_Line', 'B_ALL_Cell_Line', 
       'B_ALL_Cell_Line', 'B_ALL_Cell_Line', 'T_ALL_Cell_Line', 'T_ALL_Cell_Line', 'B_ALL_Cell_Line', 'T_ALL_Cell_Line', 'B_ALL_Cell_Line', 
       'B_ALL_Cell_Line', 'T_ALL_Cell_Line', 'B_ALL_Cell_Line', 'B_ALL_Cell_Line', 'T_ALL_Cell_Line', 'T_ALL_Cell_Line', 'B_ALL_Cell_Line',
       'B_ALL_Cell_Line', 'B_ALL_Cell_Line', 'B_ALL_Cell_Line', 'B_ALL_Cell_Line', 'T_ALL_Cell_Line', 'T_ALL_Cell_Line', 'B_ALL_Cell_Line',
       'T_ALL_Cell_Line', 'B_ALL_Patient', 'T_ALL_Cell_Line', 'AML_Cell_Line', 'Granulocytes', 'B_ALL_Patient', 'T_ALL_Cell_Line', 
       'AML_Cell_Line', 'Granulocytes', 'B_ALL_Patient', 'T_ALL_Cell_Line', 'AML_Cell_Line', 'B_ALL_Patient', 'AML_Cell_Line', 
       'AML_Cell_Line', 'B_ALL_Patient', 'AML_Cell_Line', 'AML_Cell_Line', 'B_ALL_Patient', 'AML_Cell_Line', 'AML_Cell_Line', 
       'B_ALL_Patient', 'B_ALL_Cell_Line', 'AML_Cell_Line', 'B_ALL_Patient', 'B_ALL_Cell_Line', 'AML_Cell_Line', 'B_ALL_Patient', 
       'B_ALL_Cell_Line', 'AML_Cell_Line', 'B_ALL_Patient', 'B_ALL_Cell_Line', 'B_Cells', 'B_ALL_Cell_Line', 'T_Cells', 'AML_Cell_Line', 
       'B_ALL_Patient', 'B_ALL_Cell_Line', 'Granulocytes', 'B_ALL_Cell_Line', 'Granulocytes', 'Monocytes', 'Monocytes', 
       'B_Cells', 'B_ALL_Cell_Line', 'T_Cells', 'AML_Cell_Line', 'B_ALL_Cell_Line', 'T_Cells', 'T_Cells', 'AML_Cell_Line', 'B_ALL_Cell_Line', 
       'T_Cells', 'T_Cells', 'AML_Cell_Line', 'B_Cells', 'B_ALL_Cell_Line', 'T_Cells', 'AML_Cell_Line', 'B_Cells', 'B_ALL_Cell_Line', 'T_Cells',
       'AML_Cell_Line', 'CD34', 'T_ALL_Cell_Line', 'T_ALL_Patient', 'AML_Cell_Line', 'CD34', 'T_ALL_Cell_Line', 'T_ALL_Patient', 'AML_Cell_Line', 
       'CD34', 'T_ALL_Cell_Line', 'T_ALL_Patient', 'AML_Cell_Line', 'T_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 
       'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'T_ALL_Patient', 'B_ALL_Patient', 
       'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'B_ALL_Patient', 'T_ALL_Patient', 'T_ALL_Patient', 
       'T_ALL_Patient', 'T_ALL_Patient', 'T_ALL_Patient', 'T_ALL_Patient', 'T_ALL_Patient', 'T_ALL_Patient', 'Granulocytes', 'Granulocytes', 
       'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'AML_Patient', 'AML_Patient', 'T_Cells', 
       'AML_Patient', 'AML_Patient', 'AML_Patient', 'B_Cells', 'B_Cells', 'B_Cells', 'B_Cells', 'B_Cells', 'B_Cells', 'B_Cells', 'T_Cells', 'Monocytes',
       'Monocytes', 'Monocytes', 'Monocytes', 'Granulocytes', 'T_Cells', 'T_Cells', 'T_Cells', 'T_Cells', 'T_Cells', 'T_Cells', 'T_Cells')

ex <- Biobase::exprs(gset)
# Log2 scale if it's required:
#     ex <- log2(ex + 1)
#     exprs(gset) <- ex


### Quality control:
##check if the datas are normalized:
pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()
# if we want to normalize the data: 
#     ex <- normalizeQuantiles(ex) 
#     exprs(gset) <- ex


## Correlation Heat map 
pdf("Results/CorrelationHeatmap.pdf", width = 20, height = 20)
pheatmap(cor(ex), labels_row = gr, labels_col = gr, 
         border_color = NA, bluered(256))
dev.off()


## principal component analysis
#importance of PCs
pc <- prcomp(ex)
pdf("Results/ImportanceOfPCs.pdf")
plot(pc)
dev.off()

# PC1 and PC2
pdf("Results/PC1andPC2.pdf")
plot(pc$x[,1:2])
dev.off()
# set the mean expression of all genes to zero
ex.scaled <- t(scale(t(ex), scale=FALSE))
pc <-prcomp(ex.scaled)
pdf("Results/Scaled_ImportanceOfPCs.pdf")
plot(pc)
dev.off()
pdf("Results/Scaled_PC1andPC2.pdf")
plot(pc$x[,1:2])
dev.off()


# PC of samples(we just need PC1 to PC3)
pcr <- data.frame(pc$rotation[,1:3], Group=gr)
pdf("Results/PCA samples.pdf")
ggplot(pcr, aes(x = PC1, y = PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()


#factor analysis
ex.fa.none <- factanal(ex, factors = 2, rotation = "none")
ex.fa.varimax <- factanal(ex, factors = 2, rotation = "varimax")
ex.fa.promax <- factanal(ex, factors = 2, rotation = "promax")
pdf("Results/factor_analysis.pdf")
par(mfrow = c(1,3))
plot(ex.fa.none$loadings[,1], 
     ex.fa.none$loadings[,2],
     xlab = "Factor 1", 
     ylab = "Factor 2", 
     ylim = c(-1,1),
     xlim = c(-1,1),
     main = "No rotation")
abline(h = 0, v = 0)

plot(ex.fa.varimax$loadings[,1], 
     ex.fa.varimax$loadings[,2],
     xlab = "Factor 1", 
     ylab = "Factor 2", 
     ylim = c(-1,1),
     xlim = c(-1,1),
     main = "Varimax rotation")

text(ex.fa.varimax$loadings[,1]-0.08, 
     ex.fa.varimax$loadings[,2]+0.08,
     colnames(ex),
     col="blue")
abline(h = 0, v = 0)

plot(ex.fa.promax$loadings[,1], 
     ex.fa.promax$loadings[,2],
     xlab = "Factor 1", 
     ylab = "Factor 2",
     ylim = c(-1,1),
     xlim = c(-1,1),
     main = "Promax rotation")
abline(h = 0, v = 0)
dev.off()


#SVD reduction
par(mar = rep(4, 4))
par(mfrow = c(1,3))
svdimg<-svd(ex)
U<-svdimg$u
d<-diag(svdimg$d)
V<-svdimg$v
pdf("Results/SVD_Decomposition.pdf")
plot(svdimg$d)
SVD <- svd(scale(ex))
plot(SVD$v[,1], xlab="Columns", main="Right Singular Vector (1st)", pch=19)
plot(SVD$v[,2], xlab="Columns", main="Right Singular Vector (2nd)", pch=19)
dev.off()


# t-SNE reduction
pdf("Results/t-SNE.pdf")
numTrain <- 5000
set.seed(1)
rows <- sample(1:nrow(ex), numTrain)
ex.train <- ex[rows,]
set.seed(1) 
tsne <- Rtsne(ex.train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
colors = rainbow(length(unique(gr)))
names(colors) = unique(gr)
par(mgp=c(2.5,1,0))
plot(tsne$Y, t='n', main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2, "cex.lab"=1.5)
text(tsne$Y, labels=gr, col=colors[gr])
dev.off()


### Differential expression analysis for all the normals and AML patients
gsms <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
               "XXXXXXXXXXXXXXXXXX1X1XXX1X1111X1XX11XX11X1X1X1X1X1",
               "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
               "11111111111111111111")
sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gset1 <- gset[ ,sel]

gr1 <- factor(sml)
groups <- make.names(c("AML","Healthy"))
levels(gr1) <- groups
gset1$group <- gr1
design <- model.matrix(~group + 0, gset1)
colnames(design) <- levels(gr1)
fit <- lmFit(gset1, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT, file="Results/AML_Normal.txt", row.names=F, sep="\t", quote=F)

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character((strsplit2(aml.up$Gene.symbol, "///"))))
write.table(aml.up.genes, file="Results/AML_Normal_Up.txt", row.names=F,
            col.names = F,sep="\t", quote=F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character((strsplit2(aml.down$Gene.symbol, "///"))))
write.table(aml.down.genes, file="Results/AML_Normal_Down.txt", row.names=F,
            col.names = F,sep="\t", quote=F)



### Differential expression analysis for AML and CD34
gsms <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXX22XXXXXXXXXXXXXXXXXXXX1",
               "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX00X000",
               "XXXXXXXX2222XXXXXXXX")

sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gset2 <- gset[ ,sel]

gr1 <- factor(sml)
groups <- make.names(c("AML","CD34","Monocytes"))
levels(gr1) <- groups
gset2$group <- gr1
design <- model.matrix(~group + 0, gset2)
colnames(design) <- levels(gr1)
fit <- lmFit(gset2, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT, file="Results/AML_CD34.txt", row.names=F, sep="\t", quote=F)

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character((strsplit2(aml.up$Gene.symbol, "///"))))
write.table(aml.up.genes, file="Results/AML_CD34_Up.txt", row.names=F,
            col.names = F,sep="\t", quote=F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character((strsplit2(aml.down$Gene.symbol, "///"))))
write.table(aml.down.genes, file="Results/AML_CD34_Down.txt", row.names=F,
            col.names = F,sep="\t", quote=F)


### Differential expression analysis for AML and Monocytes

cts <- paste(groups[1], groups[3], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT, file="Results/AML_Monocytes.txt", row.names=F, sep="\t", quote=F)

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character((strsplit2(aml.up$Gene.symbol, "///"))))
write.table(aml.up.genes, file="Results/AML_Monocytes_Up.txt", row.names=F,
            col.names = F,sep="\t", quote=F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character((strsplit2(aml.down$Gene.symbol, "///"))))
write.table(aml.down.genes, file="Results/AML_Monocytes_Down.txt", row.names=F,
            col.names = F,sep="\t", quote=F)