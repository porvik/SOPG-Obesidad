# Instalar y cargar paquetes necesarios #####

# Paquetes de R
# install.packages(c("RColorBrewer","R.utils"))

# Paquetes de Bioconductor
# BiocManager::install(c(
#  "Rsubread",         # Control de calidad, alineamiento y cuantificación
#  "DESeq2",           # Análisis de expresión diferencial
#  "edgeR",            # Alternativa para análisis diferencial
#  "tximport",         # Importar datos de cuantificación
#  "GenomicFeatures",  # Trabajar con anotaciones genómicas
#  "biomaRt",          # Acceso a bases de datos genómicas
#  "rtracklayer",      # Leer archivos GTF/GFF
#  "clusterProfiler",  # Análisis de enriquecimiento funcional
#  "org.Hs.eg.db"      # Anotaciones del genoma humano
#  "AnnotationDbi"     # Trabajar con anotaciones
#   "EnhancedVolcano"  # Gráficos de tipo volcano plot
# ))

library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(EnhancedVolcano)
library(BiocManager)
library(Rsubread)
library(rtracklayer)
library(clusterProfiler)
library(AnnotationDbi)


# Establecer directorio de trabajo
setwd("/home/melissa/Documents/Máster bioinformática 2025-2027/Secuenciación y ómicas de próxima generación/Actividades/Actividades  2/Actividad 2")

# Creación de carpetas para mejor organización
# dir.create("fastq_files", showWarnings = FALSE)
# dir.create("referencias", showWarnings = FALSE)
# dir.create("indices", showWarnings = FALSE)
# dir.create("alineamientos", showWarnings = FALSE)
# dir.create("resultados", showWarnings = FALSE)
# dir.create("graficos", showWarnings = FALSE)








#Cargamos metadata de las muestras ####
#Esta información relaciona la información biológica con los archivos fastq de cada individuo.

Datos_muestras <- read.csv("design1.csv")
head(Datos_muestras)
str(Datos_muestras)
summary(Datos_muestras)


# Añadimos las rutas a los archivos FASTQ
Datos_muestras$fastq_file <- file.path(
  "fastq_files",
  paste0(Datos_muestras$sample, "_R1.fastq.gz")
)

# Verificamos que todos los archivos FASTQ existen
archivos_existen <- file.exists(Datos_muestras$fastq_file)
if (!all(archivos_existen)) {
  cat("\n Los siguientes archivos no se encontraron:\n")
  print(Datos_muestras$fastq_file[!archivos_existen])
  cat("\n Revisar que los archivos FASTQ estén en la carpeta 'fastq_files'\n")
} else {
  cat("\n Todos los archivos FASTQ encontrados\n")
}

# Guardar metadata con la nueva columna de los archivos fastq de cada individuo añadida.
write.csv(Datos_muestras, "Datos_muestras.csv", row.names = FALSE)








# Descargar y preparar los archivos de referencia ####

options(timeout = 3600) # 60 minutos, hay que ampliarlo porque por defecto el tiempo máx definido no es suficiente


# Genoma humano de referencia (sólo si no lo tenemos ya)  # Lleva unos 10-20 minutos
genome_url <- "http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
genome_file <- "referencias/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
genome_file_unzip <- "referencias/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

if (!file.exists(genome_file_unzip)) {
  cat("Descargando genoma de referencia...\n")
  download.file(genome_url, genome_file, mode = "wb")
  cat("Descomprimiendo genoma...\n")
  R.utils::gunzip(genome_file, remove = FALSE)
  cat("Genoma descargado y descomprimido\n")
} else {                                               #Evitamos que lo descargue de nuevo en caso de deternlo.
  cat("Genoma de referencia ya existe\n")
}

# Descargar anotación GTF
gtf_url <- "http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
gtf_file <- "referencias/Homo_sapiens.GRCh38.110.gtf.gz"
gtf_file_unzip <- "referencias/Homo_sapiens.GRCh38.110.gtf"

if (!file.exists(gtf_file_unzip)) {
  cat("Descargando archivo de anotación GTF...\n")
  download.file(gtf_url, gtf_file, mode = "wb")
  cat("Descomprimiendo GTF...\n")
  R.utils::gunzip(gtf_file, remove = FALSE)
  cat("GTF descargado y descomprimido\n")
} else {
  cat("Archivo GTF ya existe\n")
}










# Creación del índica del genoma con Rsubread - Me ha llevado una hora aprox. ####

index_base <- "indices/grch38_index"

if (!file.exists(paste0(index_base, ".00.b.tab"))) {
  cat("Construyendo índice del genoma...\n")

  buildindex(
    basename = index_base,
    reference = genome_file_unzip,
    indexSplit = TRUE, # Para genomas grandes
    memory = 8000 # 8GB de RAM (hay que ajustarlo según el ordenador con el que se ejecute)
  )

  cat(" Índice creado\n")
} else {
  cat("Índice del genoma ya existe\n")
}

View(Datos_muestras)







# Control de calidad con FastQC y MultiQC ####

# FASTQ hace el trabajo grande: analiza cada archivo fastQ individual, calcula todas
# las metricas de calidad, genera un reporte HTML por cada muestra. Es decir procesa los datos.

# Crear carpeta para resultados
dir.create("resultados/fastqc_reports", showWarnings = FALSE, recursive = TRUE)

# Verificar que FastQC y MultiQC están instalados. Si no están hay que descargarlos en la cmd.

fastqc_instalado <- system("which fastqc", ignore.stdout = TRUE) == 0
multiqc_instalado <- system("which multiqc", ignore.stdout = TRUE) == 0

if (!fastqc_instalado) {
  cat("ERROR: FastQC no está instalado\n")
  stop("FastQC requerido para continuar")
}

if (!multiqc_instalado) {
  cat("ADVERTENCIA: MultiQC no está instalado\n")
  cat("Continuaremos solo con FastQC...\n\n")
}

# Preparar lista de archivos FASTQ
archivos_fastq <- paste(Datos_muestras$fastq_file, collapse = " ")

# Ejecutar FastQC desde R

comando_fastqc <- sprintf(
  "fastqc %s -o resultados/fastqc_reports -t 4 --quiet",
  archivos_fastq
)

resultado_fastqc <- system(comando_fastqc)

if (resultado_fastqc == 0) {
  cat("✓ FastQC completado exitosamente\n\n")

# MULTIQC: no procesa datos, sino que agrega y visualiza. Lee los resultados de fastQC y combina los reportes en uno único.
# Crea gráficos comparativos entre muestras y genera tabla resumen.

?system
  
  if (multiqc_instalado) {
    cat("Generando reporte integrado con MultiQC...\n")

    comando_multiqc <- "multiqc resultados/fastqc_reports -o resultados/fastqc_reports --force --quiet --title 'QC RNA-seq Obesidad'"

    resultado_multiqc <- system(comando_multiqc)

    if (resultado_multiqc == 0) {   #Cuando devuelve 0 es porque no ha habido ningún error, se ha ejecutado con éxito.
      cat("MultiQC completado \n")  

      cat("Carpeta de resultados: resultados/fastqc_reports/\n\n")
      cat("REPORTE PRINCIPAL (abre este archivo):\n")
      cat("multiqc_report.html\n\n")
      cat("Reportes individuales por muestra:\n")
      for (muestra in Datos_muestras$sample_id) {
        cat(sprintf("   • %s_fastqc.html\n", muestra))
      }
      cat("Archiv 'resultados/fastqc_reports/multiqc_report.html' listo para visualizar en el navegador\n\n")
    } else {
      cat("Error al ejecutar MultiQC\n")
    }
  }
} else {
  cat("Error al ejecutar FastQC\n")
  cat("Revisar que los archivos FASTQ existen\n")
}

cat("\n Control de calidad finalizado\n")

#Valoración: Las librerías presentan alta calidad, baja duplicación, contenido GC consistentem con la especie humana
# no muestran señales de contaminación. Los datos parecen aptos para alineamiento y análisis de expresión diferencial.






# Alineamiento #### Tarda unos 15 min ####

# Crear los archivos de salida de los alineamientos
Datos_muestras$bam_file <- file.path(
  "alineamientos",
  paste0(Datos_muestras$sample, ".bam")
)

# Alinear cada muestra
for (i in 1:nrow(Datos_muestras)) {
  cat(paste0(
    "\n[", i, "/", nrow(Datos_muestras), "] Alineando ",
    Datos_muestras$sample_id[i], "...\n"
  ))

  if (!file.exists(Datos_muestras$bam_file[i])) {
    align(
      index = index_base,
      readfile1 = Datos_muestras$fastq_file[i],
      output_file = Datos_muestras$bam_file[i],
      type = "rna", 
      nthreads = 4, 
      unique = TRUE, # Solo lecturas únicas
      nBestLocations = 1 # Escoger el alineamiento por lectura
    )
    cat(paste0("OK ", Datos_muestras$sample[i], " completado\n"))
  } else {
    cat(paste0("OK ", Datos_muestras$sample[i], " ya existe, no se ejecuta\n"))
  }
}

cat("\n Todos los alineamientos completados\n")

# Verificamos que todos los archivos BAM existen
if (all(file.exists(Datos_muestras$bam_file))) {
  cat("Todos los archivos BAM generados \n")
} else {
  cat("Algunos archivos BAM no se generaron\n")
}






# Cuantificación ####

# Ejecutar featureCounts para contar lecturas por gen
fc <- featureCounts(
  files = Datos_muestras$bam_file,
  annot.ext = gtf_file_unzip,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE, # Agrupar exones por gen
  allowMultiOverlap = FALSE,
  isPairedEnd = FALSE, # Cambiar a TRUE si son paired-end
  nthreads = 4,
  verbose = TRUE
)

# Extraer matriz de cuentas
counts_matrix <- fc$counts
colnames(counts_matrix) <- Datos_muestras$sample

# Ver resumen de las cuentas
print(fc$stat)

cat("\nDimensiones de la matriz de cuentas:\n")
cat(paste0("Genes: ", nrow(counts_matrix), "\n"))
cat(paste0("Muestras: ", ncol(counts_matrix), "\n"))

# Ver primeras filas
cat("\nPrimeras filas de la matriz de cuentas:\n")
print(head(counts_matrix))

# Guardar matriz de cuentas
write.csv(counts_matrix,
  "resultados/matriz_cuentas_crudas.csv",
  row.names = TRUE
)






# Mapeo de los IDs de los genes a símbolos (Ids Ensembl) ####

# Obtener símbolos de genes desde IDs de Ensembl
library(org.Hs.eg.db)

gene_ids <- rownames(counts_matrix)

# Eliminar versiones de los IDs (ENSG00000000003.14 -> ENSG00000000003)
gene_ids_clean <- gsub("\\..*", "", gene_ids)

# Mapear a símbolos
gene_symbols <- mapIds(org.Hs.eg.db,
  keys = gene_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Crear tabla de mapeo
gene_mapping <- data.frame(
  ensembl_id = gene_ids,
  ensembl_id_clean = gene_ids_clean,
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)

# Para genes sin símbolo, se usa uel ID de Ensembl
gene_mapping$gene_symbol[is.na(gene_mapping$gene_symbol)] <-
  gene_mapping$ensembl_id_clean[is.na(gene_mapping$gene_symbol)]

# Guardar mapeo
write.csv(gene_mapping,
  "resultados/mapeo_genes.csv",
  row.names = FALSE
)


# Indicación de la comparativa y filtrado de las muestras ####

# Obeso1 vs Normopeso
Datos_muestras_filtrado <- Datos_muestras %>%
  filter(condition %in% c("Obeso1", "Normopeso"))
comparativa_nombre <- "Obeso1_vs_Normopeso"
grupo_ref <- "Normopeso"
grupo_test <- "Obeso1"

cat(paste0("Comparativa: ", comparativa_nombre, "\n"))
cat(paste0("Grupo de referencia: ", grupo_ref, "\n"))
cat(paste0("Grupo a comparar: ", grupo_test, "\n"))

# Revisar que las cuentas del archivo inicial de metadatos (Datos_muestras) coinide
# con el recuento en de las muestras presentes en el conteo.
counts_filtrado <- counts_matrix[, Datos_muestras_filtrado$sample]

cat(paste0("\nMuestras en el análisis: ", ncol(counts_filtrado), "\n"))
print(Datos_muestras_filtrado[, c("sample", "condition")])

# Filtrar matriz de cuentas.
keep_genes <- rowSums(counts_filtrado >= 2) >= 2
counts_filtrado_genes <- counts_filtrado[keep_genes, ]
(counts_filtrado_genes)

# Primero se han filtrado las muestras para asegurar la correspondencia entre la
# matriz de cuentas y los metadatos experimentales. Posteriormente, se han eliminado
# los genes con baja expresión para mejorar la robustez del análisis estadístico posterior.”
