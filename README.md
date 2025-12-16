# Análisis de expresión diferencial de genes relacionados con la obesidad mediante RNA-seq

## Estructura del repositorio
El repositorio SOPG-Obesidad está organizado para facilitar un análisis reproducible de RNA-seq, separando claramente los datos de entrada, el análisis, los resultados y el material de entrega final.

- `analysis/`
Contiene el análisis reproducible en R Markdown (analysis.Rmd).
En este documento se implementa todo el pipeline del estudio, es decir importación de datos, cuantificación, normalización, análisis de expresión diferencial, visualización de resultados e interpretación biológica.

- `data/`
Incluye todos los datos necesarios para ejecutar el análisis:
    - `raw/`
    Datos originales de RNA-seq en formato FASTQ.
    - `genes/`
    Directorios correspondientes a los genes relacionados con obesidad incluidos en el estudio.
    - `Design.csv`
    Archivo de diseño experimental con la asignación de muestras a grupos metabólicos.
    - `Referencia.fasta`
    Secuencia de referencia utilizada para el mapeo o cuantificación.
    - `Transcrito_a_Gen.tsv`
    Tabla de correspondencia transcrito–gen para la agregación de conteos a nivel génico.

- `results/`
Almacena los resultados generados por el análisis:
    - `counts/`
    Matrices de conteos de expresión génica.
    - `tables/`
    Resultados estadísticos del análisis de expresión diferencial.
    - `plots/`
    Figuras generadas (volcano plots, heatmaps, etc.), listas para su uso en el póster.

- `poster/`
Contiene la plantilla oficial del póster científico utilizada para la presentación final de los resultados.

- `docs/`
Documentación adicional proporcionada en el contexto de la actividad.

- `public/`
Directorio reservado para la publicación web mediante GitLab Pages, generado automáticamente a partir del análisis con Quarto. Este directorio no se edita manualmente.

- `README.Rmd`
Documento principal del repositorio. Introduce el proyecto, describe los datos, el pipeline de análisis y la estructura del repositorio.