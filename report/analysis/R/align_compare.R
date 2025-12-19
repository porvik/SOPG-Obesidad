PROJECT_DIR <- "/home/vikto/work/bioinf/sopg/SOPG-Obesidad"

map <- read.csv(paste0(PROJECT_DIR, "/results/counts/using_R/mapeo_genes.csv"), row.names = 1)
df <- read.csv(paste0(PROJECT_DIR, "/results/counts/using_R/matriz_cuentas_crudas.csv"), row.names = 1)

df <- as.data.frame(df) |>
  tibble::rownames_to_column(var = "gene") |>
  dplyr::left_join(
    map |> dplyr::select(ensembl_id_clean, gene_symbol),
    by = c("gene" = "ensembl_id_clean")
  ) |>
  dplyr::mutate(
    label = ifelse(
      is.na(gene_symbol) | gene_symbol == gene,
      gene,
      gene_symbol
    )
  )
  
df_f <- df %>%
  filter(if_any(where(is.numeric), ~ . != 0))%>%
  select(-gene, -label)
df_f

df_2 <- read.csv(paste0(PROJECT_DIR, "/results/counts/using_linux/matriz_cuentas_crudas.csv"), row.names = 1)
df_2 <- rownames_to_column(df_2, var = "gene_symbol")
library(dplyr)

df_2_f <- df_2 %>%
  select(any_of(colnames(df_f)))

df_2_f
df_f

library(dplyr)
library(tibble)

rows <- intersect(rownames(df_f), rownames(df_2_f))
cols <- intersect(colnames(df_f), colnames(df_2_f))

df_diff <- df_f[rows, cols] - df_2_f[rows, cols]
df_diff
