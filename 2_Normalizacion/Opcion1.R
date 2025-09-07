suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(vegan)
})

# ---------- Normalizaciones ----------
norm_tss <- function(M) sweep(M, 2, colSums(M), FUN = "/")
norm_hellinger <- function(M) sqrt(norm_tss(M))
norm_clr <- function(M, pseudo = 1) {
  X <- M + pseudo
  P <- sweep(X, 2, colSums(X), "/")
  logP <- log(P)
  gm <- colMeans(logP)
  sweep(logP, 2, gm, "-")
}
to_presence_absence <- function(M) { (M > 0)*1 }

# ---------- Distancias fijas ----------
compute_distances <- function(M) {
  out <- list()
  # 1) Bray-Curtis @ TSS
  Mtss <- norm_tss(M)
  out$bray_tss <- as.matrix(vegdist(t(Mtss), method = "bray"))
  # 2) Euclídea @ Hellinger
  Mhel <- norm_hellinger(M)
  out$euclid_hellinger <- as.matrix(dist(t(Mhel), method = "euclidean"))
  # 3) Aitchison @ CLR
  Mclr <- norm_clr(M, pseudo = 1)
  out$aitchison_clr <- as.matrix(dist(t(Mclr), method = "euclidean"))
  # 4) Jaccard @ P/A
  Mpa <- to_presence_absence(M)
  out$jaccard_pa <- as.matrix(vegdist(t(Mpa), method = "jaccard", binary = TRUE))
  out
}

# ---------- Matriz de doble cero (Z) ----------
double_zero_matrix <- function(M) {
  B <- (M == 0)
  # proporción de features con 0 en ambas muestras
  denom <- nrow(M)
  Z <- matrix(0, ncol(M), ncol(M), dimnames = list(colnames(M), colnames(M)))
  for (i in seq_len(ncol(M))) {
    for (j in i:ncol(M)) {
      p <- sum(B[, i] & B[, j]) / denom
      Z[i, j] <- Z[j, i] <- p
    }
  }
  Z
}

# ---------- Wrapper para una matriz ----------
analyze_one <- function(M) {
  dists <- compute_distances(M)
  Z <- double_zero_matrix(M)
  list(distances = dists, double_zero = Z, libsize = colSums(M))
}

# ---------- Ejemplo rápido con un proyecto ----------
res <- analyze_one(lista_counts$obese_23)
str(res$distances)
range(res$double_zero)

# ---------- (opcional) correr sobre todos ----------
resultados <- imap(lista_counts, ~analyze_one(.x))






suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(purrr)
  library(tibble)
})

# --- compara una salida de analyze_one() contra Z ---
compare_dist_vs_double_zero <- function(res){
  Z <- res$double_zero
  upper_idx <- upper.tri(Z)
  map_dfr(names(res$distances), function(metodo){
    D <- res$distances[[metodo]]
    z_vec <- Z[upper_idx]
    d_vec <- D[upper_idx]
    tau <- suppressWarnings(cor(d_vec, z_vec, method = "kendall"))
    # Mantel: solo estadístico (sin permutaciones para que sea rápido)
    mant <- tryCatch(
      suppressWarnings(mantel(D, Z, method = "pearson", permutations = 0)$statistic),
      error = function(e) NA_real_
    )
    tibble(
      metodo = metodo,
      kendall_tau = round(tau, 3),
      mantel_r = round(mant, 3),
      mean_dist = round(mean(d_vec), 3),
      range_dist = paste0(round(min(d_vec),3), "–", round(max(d_vec),3))
    )
  })
}

# --- correr en toda la lista_counts ---
tabla_resumen_todos <- imap_dfr(lista_counts, function(M, nombre_proj){
  res <- analyze_one(M)  # usa tus funciones ya definidas
  compare_dist_vs_double_zero(res) %>%
    mutate(proyecto = nombre_proj, .before = 1)
})

# ordenar para leer fácil
tabla_resumen_todos <- tabla_resumen_todos %>%
  arrange(proyecto, metodo)

# ver resultados
print(tabla_resumen_todos, n = nrow(tabla_resumen_todos))

# (opcional) guardar a CSV
# write.csv(tabla_resumen_todos, file = "tabla_resumen_dist_vs_doble_cero.csv", row.names = FALSE)