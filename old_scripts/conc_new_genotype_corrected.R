## 1-15-18 update 
## corrected genotype of mouse 42013 to Ddit3

## Define output file prefix and path to avoid overwriting during different test run
out_prefix<-"libby4"
out_analysis <- "./analysis/"
out_figure <- "./figure/"
out_data <- "./data_export/"

library(tidyverse)
source("~/Documents/code/01_function/my_edgeR.R")

## load rna-seq count file and design file
rna.file = "../conc_libby_new/data/Glaucoma_all_gene_counts.txt"
design.file = "./data/design_file.txt"

data.design = read.table(design.file,  sep="\t", head=T, quote="", check.names=F)

data.raw = read.table(rna.file,  sep="\t", head=T, quote="", check.names=F, row.names=1)
colnames(data.raw)= data.design$ID_simple

group=factor(data.design$Group)

data.design %>% group_by(Group) %>% summarise(N=n())

## 
all(colnames(data.raw) == data.design$ID_simple)


## DGEList, message=FALSE
library(edgeR)
y <- DGEList(data.raw, group=group, genes=row.names(data.raw)) # must specify
options(digits=3)
y$samples

## symbols, message=FALSE
library(org.Mm.eg.db)
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype="ENSEMBL", column="SYMBOL")  # keytype="ENSEMBL", attach gene Symbol from database to DGElist y 
y$genes$Entrezid <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype="ENSEMBL", column="ENTREZID")
y$genes$Genename <- mapIds(org.Mm.eg.db, rownames(y),
                           keytype="ENSEMBL", column="GENENAME")

  ### select(org.Hs.eg.db, keys=rownames(resultTable), columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
head(y$genes)

## dropNAsymbols
y <- y[!is.na(y$genes$Symbol),]
dim(y)

## keep
keep <- rowSums(cpm(y) > 1) >= 2
table(keep)

## ----filter--------------------------------------------------------------
y <- y[keep, , keep.lib.sizes=FALSE]

## ----norm----------------------------------------------------------------
y <- calcNormFactors(y)
y$samples

## ----mdsplot # what's difference to plot log.cpm or the whole object y? 
par(xpd = T, mar = par()$mar + c(0,0,0,7))  # to make legends outside of the plot
colors <- c("darkgreen","darkgreen","red","red", "blue","blue", "purple", "purple")
pch <- c(0,15, 1, 16, 2, 17, 5, 18)
plotMDS(y, top = 500, cex = 1, pch=pch[group], dim.plot = c(1,2), ndim = 3, gene.selection = "pairwise", col=colors[group])
legend(0.83, 0.025,levels(group), pch=pch, col=colors)
par(mar=c(5, 4, 4, 2.5) + 0.1)


par(xpd = T, mar = par()$mar + c(0,0,0,7))  # to make legends outside of the plot
colors <- c("darkgreen","darkgreen","red","red", "blue","blue", "purple", "purple")
pch <- c(0,15, 1, 16, 2, 17, 5, 18)
plotMDS(y, top = 500, cex = 1, pch=pch[group], dim.plot = c(2,3), ndim = 3, gene.selection = "pairwise", col=colors[group])
legend(0.83, 0.025,levels(group), pch=pch, col=colors)
par(mar=c(5, 4, 4, 2.5) + 0.1)

## ----mdplot, fig.cap="MD plot of log2-expression in sample x versus the average log2-expression across all other samples. Each point represents a gene, and the red line indicates a log-ratio of zero. The majority of points cluster around the red line."----
## to check individual sample library after normalization
as.list(seq(ncol(y$counts))) %>% walk(plotMD_All,object=y)

## ----design--------------------------------------------------------------
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

## ----estimateDisp--------------------------------------------------------
y <- estimateDisp(y, design, robust=TRUE)

## ----plotBCV, width="3.8in", fig.cap="Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene. The plot shows the square-root estimates of the common, trended and tagwise NB dispersions."----
yticks <- seq(0, 4, 0.1)
plotBCV(y, axes = FALSE)
axis(2, at = yticks)

## ----glmQLFit------------------------------------------------------------
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

## ----QLDisp, out.width="3.8in", fig.cap="A plot of the quarter-root QL dispersion against the average abundance of each gene. Estimates are shown for the raw (before EB moderation), trended and squeezed (after EB moderation) dispersions. Note that the QL dispersions and trend shown here are relative to the NB dispersion trend shown in Figure~\ref{fig:plotBCV}."----
plotQLDisp(fit)

## ----df.prior------------------------------------------------------------
summary(fit$df.prior)

## ----MakeContrasts--------------------------------------------------------------
con <- makeContrasts(
  a1_DJvsC = (Ddit3_Jun.CONC - Ddit3_Jun.DNT) - (Control.CONC - Control.DNT),
  a2_JvsC = (Jun.CONC - Jun.DNT) - (Control.CONC - Control.DNT),
  a3_DvsC = (Ddit3.CONC - Ddit3.DNT) - (Control.CONC - Control.DNT),
  a4_DJvsJ = (Ddit3_Jun.CONC - Ddit3_Jun.DNT) - (Jun.CONC - Jun.DNT),  
  a5_DJvsD = (Ddit3_Jun.CONC - Ddit3_Jun.DNT) - (Ddit3.CONC - Ddit3.DNT),
  a6_DvsJ = (Ddit3.CONC - Ddit3.DNT) - (Jun.CONC - Jun.DNT),
  a7_CONCvsDNT = Ddit3_Jun.CONC + Jun.CONC + Ddit3.CONC + Control.CONC - (Ddit3_Jun.DNT + Jun.DNT + Ddit3.DNT + Control.DNT),
  a8_DJvsC.CONC =  Ddit3_Jun.CONC - Control.CONC, 
  a9_JvsC.CONC = Jun.CONC - Control.CONC,
  a10_DvsC.CONC = Ddit3.CONC - Control.CONC,
  a11_CONCvsDNT.DJ = Ddit3_Jun.CONC - Ddit3_Jun.DNT,
  a12_CONCvsDNT.J = Jun.CONC - Jun.DNT,
  a13_CONCvsDNT.D = Ddit3.CONC - Ddit3.DNT,
  a14_CONCvsDNT.C = Control.CONC - Control.DNT,
  levels=design
)

## ----glmQLFTest----------------------------------------------------------
res <- glmQLFTest(fit, contrast=con[,7]) # just test CONC vs DNT in all background

## ----topTags-------------------------------------------------------------
topTags(res)

## ----decideTests---------------------------------------------------------
is.de <- decideTestsDGE(res)
summary(is.de)

## ----plotMDfit, fig.cap="MD plot showing the log-fold change and average abundance of each gene. Significantly up and down DE genes are highlighted in red and blue, respectively."----
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

## ----treat---------------------------------------------------------------
tr <- glmTreat(fit, contrast=con[,7], lfc=log2(1.5))
topTags(tr)

## ----treatdecideTests----------------------------------------------------
is.de <- decideTestsDGE(tr)
summary(is.de)

## ----plotMDtreat, fig.cap="MD plot showing the log-fold change and average abundance of each gene. Genes with fold-changes significantly greater than 1.5 are highlighted."----
plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

## ----cpm-----------------------------------------------------------------
logCPM <- cpm(y, prior.count=2, log=TRUE)
logCPM.PCA<-logCPM # save it later for PCA plot
rownames(logCPM) <- y$genes$Symbol
#colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
colnames(logCPM) <- data.design$ID_simple # get it into the y project

CPM <- cpm(y, normalized.lib.sizes = TRUE, log=FALSE)
rownames(CPM) <- y$genes$Symbol
colnames(CPM) <- data.design$ID_simple
CPM_exp <- as.data.frame(CPM)
CPM_exp <- rownames_to_column(CPM_exp, var="Symbol")
# export table
write.table(CPM_exp, paste(out_data, out_prefix, "_CPM.txt", sep=""), sep="\t", row.names = F, quote=F) #export


## ----order---------------------------------------------------------------
o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:100],]

## ----scale---------------------------------------------------------------
logCPM <- t(scale(t(logCPM)))

## ----heatmap, message=FALSE, fig.width=8, fig.height=12, fig.cap="Heat map across all the samples using the top 100 most DE genes between CONC and DNT"----
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", 
          trace="none", dendrogram="both", cexRow=0.5, cexCol=0.7, density.info="none",
          margin=c(10,9))

## ---Pincicpal component analysis ---
pca_original = prcomp(t(logCPM.PCA),scale=T, center=T)
pca_x <- pca_original$x
pca_table <- data.frame(pca_x , data.design)
x <- pca_original$sdev^2/sum(pca_original$sdev^2) # Proportion of Variance Explained for all components

## Scree plot
plot(x, xlab="Principal Component", ylab="Proportion of Variance Explained", type="b")
plot(cumsum(x), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type="b")

## PCA plot --PC1 vs PC2
PCA_plot <- function(pca_table, PC_x, PC_y, color, shape){
  #PC_x,PC_y are type of interger
  #color, shape, are type of string
  g <- ggplot(pca_table, aes_string(x=names(pca_table[PC_x]), y=names(pca_table[PC_y]), color=color, shape=shape)) 
  g <- g + geom_point(alpha=0.7, size=3) 
  # g <- g + labs(color = "Group", shape="Tissue")
  g <- g + theme_bw()
  g + labs(x = paste(names(pca_table[PC_x]), scales::percent(x[PC_x]),"variance explained", sep=" "), y=paste(names(pca_table[PC_y]), scales::percent(x[PC_y]),"variance explained", sep=" "))
  #filename <- paste()
  #ggsave(filename, width=7, height=7, units="in")
}

PCA_plot(pca_table, 1, 2, "Genotype", "Treatment")
PCA_plot(pca_table, 2, 3, "Genotype", "Treatment")
ggsave("figure/PCA2_3.png", dpi=600, width=5.3, height=4, units="in")
PCA_plot(pca_table, 1, 3, "Genotype", "Treatment")


PCA_plot(pca_table, 1, 2, "Sex", "Treatment")
PCA_plot(pca_table, 2, 3, "Sex", "Treatment")
PCA_plot(pca_table, 1, 3, "Sex", "Treatment")

PCA_plot(pca_table, 1, 2, "Gen", "Treatment")
PCA_plot(pca_table, 2, 3, "Gen", "Treatment")
PCA_plot(pca_table, 1, 3, "Gen", "Treatment")

PCA_plot(pca_table, 1, 2, "Pen", "Treatment")
PCA_plot(pca_table, 2, 3, "Pen", "Treatment")
PCA_plot(pca_table, 1, 3, "Pen", "Treatment")

## ----paired glmQLFtest---
print(con)
#as.list(seq(ncol(con)))%>% walk(my_glmTreat,fit=fit, Treat_FC=1)
#as.list(seq(ncol(con)))%>% walk(my_glmTreat,fit=fit, Treat_FC=1.2)
as.list(seq(ncol(con)))%>% walk(my_glmTreat,fit=fit, Treat_FC=1.5)

## ----anovaQLFtest--------------------------------------------------------
### Do not bother performing kegga and goana function, there multiple logFC made. It doesn't fit kegga and goana function which only takes one logFC.

contrast_i<-list(c(1:3), c(8,9,10), c(11,12,13,14))
contrast_i %>% walk(my_anovaQLF_FDR, fit=fit)

sessionInfo()
