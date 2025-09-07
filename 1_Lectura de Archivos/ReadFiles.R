library(tidyverse)

# Definir carpeta donde guardaste todos los TSV
my_wd <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables"

# Listar todos los .tsv
tsv_files <- list.files(my_wd, pattern = "\\.tsv$", full.names = TRUE)

# Función para leer y convertir a matriz con rownames
leer_tsv_abund <- function(filepath){
  read_tsv(filepath, col_types = cols()) %>%
    column_to_rownames(var = colnames(.)[1]) %>%  # primera columna como nombres de fila
    as.matrix()
}

# Leer todo en una lista nombrada con el nombre del proyecto (sin extensión)
lista_counts <- set_names(
  lapply(tsv_files, leer_tsv_abund),
  nm = tools::file_path_sans_ext(basename(tsv_files))
)

# Ejemplo: acceder al conteo del proyecto "106_infant_singapore"
counts_106 <- lista_counts[["106_infant_singapore"]]

# Ver dimensiones y un preview
dim(counts_106)
head(counts_106[, 1:5])
