#!/bin/bash

# 1. Indexado (Ya lo has hecho, pero lo dejo aquí por si necesitas repetirlo)
#salmon index -t Referencia.fasta -i index -p 4

# 2. Cuantificación
# Creamos una lista con los nombres de todos los personajes
PERSONAJES="AbrahamSimpson BartSimpson HomerSimpson LisaSimpson MaggieSimpson MargeSimpson PattyBouvier SelmaBouvier"

# Iniciamos el bucle: Para cada personaje en la lista...
for PERSONAJE in $PERSONAJES
do
    echo "Procesando a: $PERSONAJE"
    
    salmon quant -i index -l A \
    -1 Fastqs/${PERSONAJE}_R1.fastq.gz \
    -2 Fastqs/${PERSONAJE}_R2.fastq.gz \
    -p 4 \
    -o Salmon/${PERSONAJE}
    
done

echo "¡Proceso terminado para todos los Simpsons!"