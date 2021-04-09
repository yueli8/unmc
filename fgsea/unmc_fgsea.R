library(fgsea)
library(tidyverse)
library(data.table)
library(ggplot2)
setwd("/Users/yli1/unmc01")
#hcq_vs_con
res<-read.table("hcq_vs_con_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fwrite(fgseaRes_kegg, file="hcq_vs_con_kegg01.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="hcq_vs_con_all01.txt", sep="\t", sep2=c("", " ", ""))

plotEnrichment(pathways.kegg[["KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS"]],
               ranks) + labs(title="SYSTEMIC_LUPUS_ERYTHEMATOSUS")


#pcq_vs_con
res<-read.table("pcq_vs_con_fgsea.txt")
ranks <- deframe(res)
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fwrite(fgseaRes_kegg, file="pcq_vs_con_kegg01.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="pcq_vs_con_all01.txt", sep="\t", sep2=c("", " ", ""))


plotEnrichment(pathways.kegg[["KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS"]],
               ranks) + labs(title="SYSTEMIC_LUPUS_ERYTHEMATOSUS")

plotEnrichment(pathways.kegg[["KEGG_LEISHMANIA_INFECTION"]],
               ranks) + labs(title="LEISHMANIA_INFECTION")


#phpma_vs_con
res<-read.table("phpma_vs_con_fgsea.txt")
ranks <- deframe(res)
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fwrite(fgseaRes_kegg, file="phpma_vs_con_kegg01.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="phpma_vs_con_all01.txt", sep="\t", sep2=c("", " ", ""))

