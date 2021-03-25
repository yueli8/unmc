library(DESeq2)
library(data.table)
library(limma)
library(calibrate)
setwd("/media/hp/04c65089-71ff-4b33-9a30-c21b8c77eda2/li/unmc")
#hcq_vs_con
cts<-read.table("hcq_vs_con.txt",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("hcq",3),rep("con",3)), levels = c("hcq","con"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("hcq","con"))
dds$condition <- relevel(dds$condition, ref = "con")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"hcq_vs_con_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","hcq","con"))
resultsNames(dds2)
write.csv(res,"hcq_vs_con_deg.csv")


#pcq_vs_con
cts<-read.table("pcq_vs_con.txt",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("pcq",3),rep("con",3)), levels = c("pcq","con"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("pcq","con"))
dds$condition <- relevel(dds$condition, ref = "con")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"pcq_vs_con_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","pcq","con"))
resultsNames(dds2)
write.csv(res,"pcq_vs_con_deg.csv")


#phpma_vs_con
cts<-read.table("phpma_vs_con.txt",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("phpma",3),rep("con",3)), levels = c("phpma","con"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("phpma","con"))
dds$condition <- relevel(dds$condition, ref = "con")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"phpma_vs_con_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","phpma","con"))
resultsNames(dds2)
write.csv(res,"phpma_vs_con_deg.csv")

#volcano plot pcq_vs_con_DEG
res <- read.csv("pcq_vs_con_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

#volcano plot pcq_vs_con_DEG
res <- read.csv("pcq_vs_con_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6.5,5)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

#volcano plot phpma_vs_con_DEG
res <- read.csv("phpma_vs_con_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,5)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)


library(fgsea)
library(tidyverse)
library(data.table)
library(ggplot2)
#hcq_vs_con
res<-read.table("hcq_vs_con_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fwrite(fgseaRes_kegg, file="hcq_vs_con_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="hcq_vs_con_all.txt", sep="\t", sep2=c("", " ", ""))

#pcq_vs_con
res<-read.table("pcq_vs_con_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fwrite(fgseaRes_kegg, file="pcq_vs_con_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="pcq_vs_con_all.txt", sep="\t", sep2=c("", " ", ""))

#phpma_vs_con
res<-read.table("phpma_vs_con_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fwrite(fgseaRes_kegg, file="phpma_vs_con_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="phpma_vs_con_all.txt", sep="\t", sep2=c("", " ", ""))