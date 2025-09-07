library(readr)

path_tsv <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/23_obese.tsv"  # poné acá la ruta correcta

# 1) localizar la fila de cabecera real (la que empieza con "#OTU ID" o "OTU ID")
hdr_probe <- read_lines(path_tsv, n_max = 500)
hdr_idx <- which(grepl("^#?OTU[ _]?ID(\\t|,)", hdr_probe))[1]
stopifnot(!is.na(hdr_idx))

# 2) leer el TSV desde esa fila
M_raw <- read_tsv(path_tsv, skip = hdr_idx - 1, col_types = cols())

# 3) usar la 1ª columna como rownames (OTU IDs), quitando el "#"
names(M_raw)[1] <- sub("^#", "", names(M_raw)[1])
rownames(M_raw) <- M_raw[[1]]
M_raw <- as.data.frame(M_raw[ , -1, drop = FALSE])

# 4) asegurarnos de que todo sea numérico
M_raw[] <- lapply(M_raw, function(x) as.numeric(as.character(x)))

# 5) transponer: filas = muestras, columnas = OTUs
M <- t(as.matrix(M_raw))

# chequeo rápido
cat("Dimensiones M (muestras x OTUs):", nrow(M), "x", ncol(M), "\n")
cat("Ejemplo rownames (muestras):", paste(head(rownames(M)), collapse=", "), "\n")
cat("Ejemplo colnames (OTUs):", paste(head(colnames(M)), collapse=", "), "\n")
