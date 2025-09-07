library(vegan)
library(phyloseq)   # si estás usando objetos ps, si no lo omitimos

# 1. Tomar abundancias (POST_NZV que sabemos que funciona)
X <- abund_post_nzv_list[["95_ovese"]]   # tu tabla de OTUs filtrada
y <- groups_list[["95_ovese"]]           # los grupos

# 2. Armar metadata
metadata <- data.frame(Group = y)
rownames(metadata) <- rownames(X)  # asegurar matching

# 3. Calcular distancia Bray-Curtis
dist_bc <- vegdist(X, method = "bray")

# 4. PERMANOVA
perm <- adonis2(dist_bc ~ Group, data = metadata, permutations = 9999)
print(perm)

# 5. Homogeneidad de dispersión (beta-dispersion)
bd <- betadisper(dist_bc, metadata$Group)
anova(bd)




#Tengi que sacar las tablas de paso 1 normalizacion
library(vegan)


proj <- "95_ovese"

# --- Metadata (manifest) normalizada ---
md <- manifest_list[[proj]]
stopifnot(!is.null(md))
md$sample.id <- gsub("-", ".", md$sample.id)         # 'sample-19' -> 'sample.19'
metadata_95 <- data.frame(Group = factor(md$group),  # factor de grupos
                          row.names = md$sample.id)

# --- Matriz de abundancias tal cual la tenés ---
M <- as.data.frame(ovese_95)

# Normalizar nombres en matriz (tanto filas como columnas)
if (!is.null(rownames(M))) rownames(M) <- gsub("-", ".", rownames(M))
if (!is.null(colnames(M))) colnames(M) <- gsub("-", ".", colnames(M))

# Helper: pasar a relativas y quitar filas suma 0
to_rel <- function(A){
  A[is.na(A)] <- 0
  rs <- rowSums(A)
  A <- A[rs > 0, , drop = FALSE]
  sweep(A, 1, pmax(rowSums(A), 1e-12), "/")
}

# Detectar orientación: ¿muestras en filas o columnas?
ids_row <- intersect(rownames(M), rownames(metadata_95))
ids_col <- intersect(colnames(M), rownames(metadata_95))

if (length(ids_row) >= 3) {
  # Muestras en filas (raro en tu captura, pero lo contemplamos)
  A <- M[ids_row, , drop = FALSE]
  meta <- metadata_95[ids_row, , drop = FALSE]
} else if (length(ids_col) >= 3) {
  # Muestras en columnas -> transponer
  A <- t(M)
  # después de transponer, los rownames son los nombres de muestra (ya normalizados arriba)
  ids <- intersect(rownames(A), rownames(metadata_95))
  A <- A[ids, , drop = FALSE]
  meta <- metadata_95[ids, , drop = FALSE]
} else {
  stop("Sin IDs en común aún. Revisá que el manifest tenga IDs tipo 'sample-xx'/'sample.xx' como en la matriz.")
}

# Relativas y alineación final
A_rel <- to_rel(A)
meta  <- meta[rownames(A_rel), , drop = FALSE]

# Requisitos mínimos
if (nrow(A_rel) < 3 || length(unique(meta$Group)) < 2) {
  stop("Quedan <3 muestras o <2 grupos tras alinear; revisá el manifest y los IDs.")
}

# Bray–Curtis + PERMANOVA + betadisper
dist_bc <- vegdist(A_rel, method = "bray")
perm    <- adonis2(dist_bc ~ Group, data = meta, permutations = 9999)
bd      <- betadisper(dist_bc, meta$Group)

print(perm)
print(anova(bd))





library(vegan)

# 1) Matriz original
mat <- as.data.frame(ovese_95)

# 2) Filtrar a las 20 OTUs más abundantes
top20_taxa <- names(sort(colSums(mat), decreasing = TRUE))[1:20]
mat_top20 <- mat[, top20_taxa]

# 3) Abundancias relativas
mat_rel20 <- mat_top20 / rowSums(mat_top20)

# 4) Metadata (usamos la misma que antes: meta)
# asegurate que tenga rownames = nombres de las muestras
meta <- metadata_95[rownames(mat_rel20), , drop = FALSE]

# 5) Distancias Bray-Curtis
dist_bc20 <- vegdist(mat_rel20, method = "bray")

# 6) PERMANOVA
perm20 <- adonis2(dist_bc20 ~ Group, data = meta, permutations = 9999)
print(perm20)

# 7) Homogeneidad de dispersión
bd20 <- betadisper(dist_bc20, meta$Group)
print(anova(bd20))