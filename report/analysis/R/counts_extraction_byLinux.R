library(readr)
library(tximport)

sample_info = read_csv("Design.csv", locale = locale(encoding = "UTF-8"))
tx2gene = read_tsv("Transcrito_a_Gen.tsv", col_names = FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")

# Definimos rutas a los ficheros de cuantificación de Salmon
files = file.path("Salmon", sample_info$Sample, "quant.sf")
names(files) = sample_info$Sample

# Leemos los datos de expresión con tximport
txi = tximport(files, type = "salmon", tx2gene = tx2gene)

# Exploración de la matriz
conteos = txi$counts
conteos
