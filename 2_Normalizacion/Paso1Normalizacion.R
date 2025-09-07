library(readr)
library(dplyr)
library(purrr)
library(tibble)

#Funcion que convierte archivo tsv a matriz.
leer_biom_tsv <- function(path) {
  # Leer, ignorando la línea de comentario del biom
  df <- read_tsv(path, comment = "#", col_types = cols(.default = col_guess()))
  stopifnot(ncol(df) >= 2)
  
  id_col <- names(df)[1]           # "OTU ID"
  ids <- as.character(df[[id_col]])
  if (anyDuplicated(ids) > 0) ids <- make.unique(ids)
  
  # Quitar la columna de IDs usando indexación base
  M <- as.data.frame(df[ , -1, drop = FALSE])
  
  # Coerción numérica segura
  M[] <- lapply(M, function(x) suppressWarnings(as.numeric(x)))
  if (any(!vapply(M, is.numeric, logical(1)))) {
    stop("Hay columnas no numéricas inesperadas.")
  }
  
  M[is.na(M)] <- 0
  rownames(M) <- ids
  as.matrix(M)
}

#Lo usamos para cada uno de los proyectos para obtener la matriz de interes.
path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/23_obese.tsv"   
obese_23 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/29_children_eu_africa.tsv"
children_eu_africa_29 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/34_L_H.tsv"
L_H_34 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/95_ovese.tsv"
ovese_95 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/98_combo_diet.tsv"
combo_diet_98 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/106_infant_singapore.tsv"
infant_singapore_106 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/126rural_af_urban_usa.tsv"
rural_af_urban_usa_126 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/147_cameroon.tsv"
cameroon_147 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/153_italian.tsv"
italian_153 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/161_lcarb_lfat.tsv"
lcarb_lfat_161 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/174_healthy_fibers.tsv"
healthy_fibers_174 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/192_infant.tsv"
infant_192 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/248_citizen.tsv"
citizen_248 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/346_asian.tsv"
asian_346 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/350_hadza.tsv"
hadza_350 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/367_healthy_japanese.tsv"
healthy_japanese_367 <- leer_biom_tsv(path)

path <- "/Users/sofioiene/Desktop/double_zero-main/projects_count_tables/471_swedish.tsv"
swedish_471 <- leer_biom_tsv(path)



# Paso todas las matrices a una lista
lista_counts <- list(
  obese_23 = obese_23,
  children_eu_africa_29 = children_eu_africa_29,
  L_H_34 = L_H_34,
  ovese_95 = ovese_95,
  combo_diet_98 = combo_diet_98,
  infant_singapore_106 = infant_singapore_106,
  rural_af_urban_usa_126 = rural_af_urban_usa_126,
  cameroon_147 = cameroon_147,
  italian_153 = italian_153,
  lcarb_lfat_161 = lcarb_lfat_161,
  healthy_fibers_174 = healthy_fibers_174,
  infant_192 = infant_192,
  citizen_248 = citizen_248,
  asian_346 = asian_346,
  hadza_350 = hadza_350,
  healthy_japanese_367 = healthy_japanese_367,
  swedish_471 = swedish_471
)
