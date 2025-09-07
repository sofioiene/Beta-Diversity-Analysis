library(vegan)
library(pheatmap)

# =========================
# 0) Configuración
# =========================
# Usamos raw + 10 normalizaciones (TSS incluida como "total")
norm_methods <- c(
  "raw",        # sin normalizar (RO)
  "total",      # TSS (ideal para Bray–Curtis)
  "max",
  "frequency",
  "normalize",
  "chi.square",
  "hellinger",
  "log",
  "range"
)

# Asegurar que M es numérica y sin NA
M <- as.matrix(M)
storage.mode(M) <- "double"
M[!is.finite(M)] <- 0

# =========================
# 1) Normalizaciones
# =========================
apply_norm <- function(mat, method) {
  if (method == "raw") {
    out <- mat
  } else if (method == "log") {
    out <- log1p(mat)                         # mantiene no-negatividad
  } else {
    out <- decostand(mat, method = method)    # total, max, frequency, pa, normalize, range, chi.square, hellinger
  }
  out[!is.finite(out)] <- 0
  # Por si alguna normalización deja un negativo numérico muy pequeño:
  mn <- suppressWarnings(min(out, na.rm = TRUE))
  if (is.finite(mn) && mn < 0) out <- out - mn + 1e-12
  out
}

# Pseudoconteo pequeño para evitar filas suma 0 tras normalizar
eps <- 1e-8
M_eps <- M + eps

norm_mats <- lapply(setNames(norm_methods, norm_methods), function(m) apply_norm(M_eps, m))

# Todas comparten EXACTAMENTE el mismo set y orden de muestras
common_ids <- Reduce(intersect, lapply(norm_mats, rownames))
if (length(common_ids) < 3) stop("Quedan <3 muestras comunes; revisá orientación/IDs.")
norm_mats <- lapply(norm_mats, function(A) A[common_ids, , drop = FALSE])

# =========================
# 2) Distancias Bray–Curtis
# =========================
dist_list <- lapply(norm_mats, function(A) vegdist(A, method = "bray"))

# Sanity check
if (any(vapply(dist_list, function(d) any(is.na(d)), logical(1)))) {
  stop("Alguna distancia tiene NA; revisá filas con suma 0 o valores no finitos.")
}

# =========================
# 3) Mantel all-vs-all
# =========================
mantel_mat <- matrix(NA_real_, length(norm_methods), length(norm_methods),
                     dimnames = list(norm_methods, norm_methods))

for (i in seq_along(dist_list)) {
  for (j in seq_along(dist_list)) {
    mantel_mat[i, j] <- if (i == j) 1 else unname(mantel(dist_list[[i]], dist_list[[j]], permutations = 999)$statistic)
  }
}

print(round(mantel_mat, 3))

# =========================
# 4) Heatmap
# =========================
pheatmap(mantel_mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "Mantel (Bray–Curtis) — 23_obese")
