
# Actividad 2
# DEG analysis using DSeq2

setwd("/mnt/Datos 1/BIOINFO/UNIR-Bioinfo/Secuenciación y ómicas/Actividades/Actividad2/TallerGrupal_Ficheros/")

if (!require("BiocManager", quietly = TRUE)) # Verificamos si el paquete está instalado sin cargarlo,
  install.packages("BiocManager") # Si no está, se instala.
BiocManager::install("DESeq2", force = TRUE) # Instalamos el programa que vamos a utilizar para realizar el análisis diferencial de genes
BiocManager::install("tximport")
library(DESeq2)
library(dplyr)
library(readr)
library(tximport)

sample_info = read_csv("Design.csv", locale = locale(encoding = "UTF-8"))
tx2gene = read_tsv("Transcrito_a_Gen.tsv", col_names = FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")

# Definimos rutas a los ficheros de cuantificación de Salmon
files = file.path("Salmon", sample_info$Sample, "quant.sf")
names(files) = sample_info$Sample

# Leemos los datos de expresión con tximport
txi = tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")

# Exploración de la matriz
conteos = round(txi$counts)
conteos

head(conteos)

# rownames(conteos) = conteos$Gene.ID
# head(conteos)

# genes = counts[, c("Gene.ID", "Gene.Name")]
conteos = conteos[, -c(3, 7, 8)]
head(conteos)

head(sample_info)
rownames(sample_info) = sample_info$Sample
head(sample_info)

sample_info$Condition[sample_info$Condition == 'Sobrepeso/Obeso1'] = 'Obeso1'
sample_info$Condition[sample_info$Condition == 'Sobrepeso/Obeso2'] = 'Obeso2'
sample_info

# Filtramos para quedarnos sólo con las condiciones a comparar
sample_info2 <- sample_info %>%
  filter(Condition %in% c("Normopeso", "Obeso1"))
sample_info2

# Definimos como factor a Condition y especificamos los niveles (en nuestro caso, sólo compararemos Normopeso y Obeso1)
sample_info2$Condition = factor(sample_info2$Condition, levels=c("Normopeso", "Obeso1"))
sample_info2$Condition


# DEG analysis

## No aplicamos filtros ya que trabajamos con valores de expresión muy bajos en este ejemplo

dds <- DESeqDataSetFromMatrix(countData=conteos, colData=sample_info2, design=~Condition)

dds <- DESeq(dds)

## Extraemos la información de los resultados específicos comparando los grupos Normopeso y Obeso1
res = results(dds, contrast=c("Condition", "Obeso1", "Normopeso"), alpha=1e-5) 
res


#Convertimos res en data frame:
res_df = as.data.frame(res)


# Visualización de los datos

## Es necesario crear una columna con el nombre de los genes (antes los teníamos como row names), la llamamos Gene.ID

res_df$Gene.ID <- rownames(res_df)
res_df

## MA plot
plotMA(res) # Los valores en azul corresponden a los genes que pasan el umbral que hemos especificado del alpha

## Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

EnhancedVolcano(
  res_df,
  lab=res_df$Gene.ID,
  x='log2FoldChange',
  y='pvalue',
  labSize = 3,
  axisLabSize = 10,
  drawConnectors = TRUE,   
  widthConnectors = 0.5,
  )

## HeatMap

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}

vsd <- varianceStabilizingTransformation(dds, blind = FALSE) ## Aplica la transformación variance-stabilizing transformation al objeto dds
mat <- assay(vsd)[(rownames(res_df)), ] ## Extrae la matriz de expresión transformada (filas = genes, columnas = muestras)
rownames(mat) <- res_df$Gene.ID ## Reemplaza los identificadores de fila por nombres de genes
mat_scaled <- t(scale(t(mat))) # Escala y transpone la matriz


### Extracción de lista top genes diferencialmente expresados

### Ordenar por significancia ajustada
res_ordered <- res_df[order(res_df$padj), ]
### Seleccionar los 10 genes más significativos
top10_genes <- head(res_ordered, 10)
top10_genes # aquí observamos que tansolo 8 genes cumplen valor p < 0.05

### Para resaltar los genes significativos:

genes_resaltar <- c("NTRK2", "LEP", "LEPR", "MC4R", "PCSK1", "SH2B1", "BDNF", "CADM2")

annotation_row <- data.frame(
  Significant = ifelse(rownames(mat_scaled) %in% genes_resaltar, "Sí", "No")
)

rownames(annotation_row) <- rownames(mat_scaled)

ann_colors <- list(
  Significant = c("Sí" = "yellow", "No" = "white")
)


### Heatmap incluyendo todos los genes para ver patrones (se resaltan en amarillo los significativos):

pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))


### Heatmap incluyendo únicamente los genes significativos

mat_subset <- mat_scaled[(head(order(res_df$padj), 8)), ] # en nuestro caso seleccionamos sólo los significativos padj<0.05, que son 8 como hemos visto anteriormente

pheatmap(mat_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50))



# Análisis de enriquecimiento

pkgs_bioc <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot")
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)
  library(p, character.only = TRUE)
}


## Creamos un subconjunto de datos con símbolos únicos y los resultados de significancia
genes_df <- unique(data.frame(symbol = res_df$Gene.ID,
                              log2FC = res_df$log2FoldChange,
                              padj = res_df$padj,
                              stringsAsFactors = FALSE))

## Eliminamos aquellos genes que no tienen nomenclatura, podríamos eliminar incluso con algún valor faltante
genes_df <- genes_df[!is.na(genes_df$symbol) & genes_df$symbol != "", ]

## Asociamos cada nombre de gen con el identificador de la base de datos ENTREZID
map <- bitr(genes_df$symbol,
            fromType = "SYMBOL",
            toType   = c("ENTREZID"),
            OrgDb    = org.Hs.eg.db)

## Unimos ambos datos
genes_mapped <- merge(genes_df, map, by.x = "symbol", by.y = "SYMBOL")
nrow(genes_df); nrow(genes_mapped) # Comprobamos duplicados: un nombre puede mapear a >1 ENTREZ (rara vez)
universe_entrez <- unique(map$ENTREZID) # Creamos el objeto de elementos únicos para los análisis

## Over-Representation Analysis (ORA)

sig_genes <- genes_mapped[genes_mapped$padj < 0.05 & abs(genes_mapped$log2FC) >= 0.25, ]
length(sig_genes$ENTREZID)

### Con threshold de 1 para log2FC, se obtienen 4 genes con expresión diferencialmente significativa
### Con threshold de 0.25 para log2FC, se obtienen 8 genes con expresión diferencialmente significativa
### Un log2FC de +/-0.25 implica un cambio de expresión de aproximadamente el 19% consideramos que biológicamente puede tener efecto para ciertos genes

# ORA: enrichGO (Biological Process)

ego_bp <- enrichGO(gene = sig_genes$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",                 # BP (proceso biológico), MF (función molecular), CC (componente celular)
                   pAdjustMethod = "BH",
                   universe = universe_entrez,
                   minGSSize = 3,
                   maxGSSize = 500,
                   pvalueCutoff = 0.2,
                   qvalueCutoff = 0.2)


dotplot(ego_bp, showCategory = 20) + ggtitle("GO:BP enrichment (ORA)") + theme(axis.text.y=element_text(size=8))


e_reactome <- enrichPathway(gene = sig_genes$ENTREZID,
                            organism = "human",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = universe_entrez,
                            minGSSize = 3,
                            qvalueCutoff = 0.2)

dotplot(e_reactome, showCategory = 20) + ggtitle("Reactome enrichment (ORA)") + theme(axis.text.y=element_text(size=6))


ego_bp_symbol <- setReadable(ego_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cnetplot(ego_bp_symbol, showCategory = 6)



# GSEA

gene_map <- bitr(unique(res_df$Gene.ID),
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

genes_mapped <- merge(res_df, gene_map, by.x = "Gene.ID", by.y = "SYMBOL")

# Eliminamos duplicados por ENTREZID (manteniendo el valor de mayor |log2FC|)
genes_mapped <- genes_mapped[order(abs(genes_mapped$log2FoldChange), decreasing = TRUE), ]
genes_mapped <- genes_mapped[!duplicated(genes_mapped$ENTREZID), ]

# Construimos geneList nombrado por ENTREZID
geneList <- genes_mapped$log2FoldChange
names(geneList) <- genes_mapped$ENTREZID

set.seed(123)
geneList <- geneList + rnorm(length(geneList), mean = 0, sd = 1e-6)
geneList <- sort(geneList, decreasing = TRUE)


## GSEA - GO

gsea_go <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  minGSSize = 4,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.8,
  eps = 0,
  nPermSimple = 10000,
  verbose = TRUE
)

head(as.data.frame(gsea_go))

library(enrichplot)

dotplot(gsea_go, showCategory = 20) +
  ggtitle("GSEA GO:CC") +
  theme(axis.text.y = element_text(size = 10))

ridgeplot(gsea_go) + ggtitle("GSEA GO:CC") + theme(axis.text.y = element_text(size = 8))


## GSEA - REACTOME

gsea_reactome <- gsePathway(
  geneList = geneList,
  organism = "human",
  minGSSize = 3,
  pvalueCutoff = 0.9,
  pAdjustMethod = "BH",
  eps = 0
)

dotplot(gsea_reactome, showCategory = 20) +
  ggtitle("GSEA Reactome") +
  theme(axis.text.y = element_text(size = 8))

ridgeplot(gsea_reactome) + ggtitle("GSEA Reactome") + theme(axis.text.y = element_text(size = 8))

