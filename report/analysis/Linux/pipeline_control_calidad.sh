#!/bin/bash

# Crear entorno y activarlo
# conda create -n simpsons_rna seq fastqc multiqc fastp salmon -y
conda activate simpsons_rna

#Situarse en el directorio de fastqc para el control de calidad

mkdir -p Resultados/Quality/RAW Resultados/Quality/Filtered Resultados/Trimmed
echo "Empieza el analisis de calidad"

fastqc *fastq.gz -o Resultados/Quality/RAW
#Obtener los identificadores de cada muestra
ls *fastq.gz | cut -d _ -f 1 | sort -u > muestras.txt
#Realizar un bucle para filtrar
for i in $(cat muestras.txt); do fastp --in1 $i*R1* --in2 $i*R2* --out1 Resultados/Trimmed/$i"_R1_filtered.fastq.gz" --out2 Resultados/Trimmed/$i"_R2_filtered.fastq.gz" --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 12 --cut_mean_quality 30 --length_required 35 --json Resultados/Trimmed/$i.json --html Resultados/Trimmed/$i.html --thread 32; done
#Control de calidad de los ficheros filtrados: 
fastqc Resultados/Trimmed/*fastq.gz -o Resultados/Quality/Filtered/ --threads 32
# 8. Podemos generar un informe resumen de los pasos que hemos realizado con multiqc:
multiqc Resultados -o Resultados/Reporte_final