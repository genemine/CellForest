# load R packages
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# GO enrichment analysis
g=genes
# ENSEMBL SYMBOL
g<- g %>%  bitr(fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
#GO
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = g$ENTREZID, 
                ont = "all", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
save(ego,file = "D:/enrich.Rdata")

head(ego)
dotplot(ego,showCategory=30,title="Top 30 Enriched GO Terms")
barplot(ego, showCategory=20,title="EnrichmentGO")
