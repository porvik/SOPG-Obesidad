# Análisis de enriquecimiento

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

## =========================
## MAPEAR SYMBOL -> ENTREZID
## =========================
res_df$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys = res_df$Gene.Name,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

## =========================
## FILTRAR GENES CON ENTREZID VÁLIDO
## =========================
res_df_filtered <- res_df[!is.na(res_df$ENTREZID), ]

## =========================
## PREPARAR RANKING PARA GSEA
## =========================
gene_list <- res_df_filtered$log2FoldChange
names(gene_list) <- res_df_filtered$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)  # orden descendente

## =========================
## GSEA: GO Biological Process
## =========================
gsea_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 1,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)

## =========================
## GSEA: Reactome
## =========================
gsea_reactome <- gsePathway(
  geneList = gene_list,
  organism = "human",
  minGSSize = 1,
  pvalueCutoff = 1,
  verbose = FALSE
)

## =========================
## MOSTRAR RESULTADOS EN CONSOLA
## =========================
cat("Top GO BP GSEA:\n")
print(head(as.data.frame(gsea_bp), 10))

cat("\nTop Reactome GSEA:\n")
print(head(as.data.frame(gsea_reactome), 10))

## =========================
## GUARDAR RESULTADOS EN CSV
## =========================
write.csv(as.data.frame(gsea_bp), "GSEA_GO_BP.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_reactome), "GSEA_Reactome.csv", row.names = FALSE)

## =========================
## VISUALIZACIONES
## =========================
print(dotplot(gsea_bp, showCategory = 20) +
        ggtitle("GO Biological Process GSEA") +
        theme(axis.text.y = element_text(size = 8)))

print(dotplot(gsea_reactome, showCategory = 20) +
        ggtitle("Reactome GSEA") +
        theme(axis.text.y = element_text(size = 7)))

## =========================
## CNETPLOT PARA TOP GO BP
## =========================
gsea_bp_symbol <- setReadable(gsea_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Solo mostrar cnetplot si hay términos enriquecidos
if(nrow(as.data.frame(gsea_bp_symbol)) > 0){
  print(cnetplot(gsea_bp_symbol, showCategory = 6))
} else {
  cat("\nNo hay términos enriquecidos suficientes para cnetplot.\n")
}
