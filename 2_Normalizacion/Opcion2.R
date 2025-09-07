suppressPackageStartupMessages({
  library(vegan)
  library(dplyr); library(purrr); library(tibble)
})

# --- Z rápida (doble cero) ---
double_zero_matrix_fast <- function(M){
  B <- (M == 0L)                  # lógico OK; crossprod cuenta TRUE como 1
  Z <- crossprod(B) / nrow(M)     # t(B) %*% B  -> muestras x muestras
  cn <- colnames(M); if (is.null(cn)) cn <- paste0("s", seq_len(ncol(M)))
  dimnames(Z) <- list(cn, cn)
  as.matrix(Z)
}

# --- Si querés, redefiní analyze_one para usar la Z rápida ---
analyze_one_fast <- function(M, pseudo_clr = 1){
  dists <- compute_distances(M)         # usa tus normalizaciones/distancias ya definidas
  Z <- double_zero_matrix_fast(M)
  list(distances = dists, double_zero = Z, libsize = colSums(M))
}

# --- Tabla resumen para un proyecto ---
compare_dist_vs_double_zero <- function(res){
  Z <- res$double_zero
  U <- upper.tri(Z)
  map_dfr(names(res$distances), function(metodo){
    D <- res$distances[[metodo]]
    z_vec <- Z[U]; d_vec <- D[U]
    tau <- suppressWarnings(cor(d_vec, z_vec, method = "kendall"))
    mant <- tryCatch(
      suppressWarnings(mantel(D, Z, method = "pearson", permutations = 0)$statistic),
      error = function(e) NA_real_
    )
    tibble(
      metodo = metodo,
      kendall_tau = round(tau, 3),
      mantel_r    = round(mant, 3),
      mean_dist   = round(mean(d_vec), 3),
      range_dist  = paste0(round(min(d_vec),3), "–", round(max(d_vec),3))
    )
  })
}

# ====== USO RÁPIDO (un proyecto) ======
res_obese23 <- analyze_one_fast(lista_counts$obese_23)
tabla_obese23 <- compare_dist_vs_double_zero(res_obese23)
print(tabla_obese23)

# ====== (Opcional) TODOS LOS PROYECTOS, con submuestreo defensivo ======
# Para evitar cuelgues en datasets con >500 muestras, submuestreamos a 400
proyectos <- names(lista_counts)
tabla_resumen_todos <- list()

for(nm in proyectos){
  M <- lista_counts[[nm]]
  if(ncol(M) > 500){
    set.seed(123)
    keep <- sample(colnames(M), 400)   # ajustá 400 según tu RAM
    M <- M[, keep, drop = FALSE]
  }
  res <- analyze_one_fast(M)
  tabla_resumen_todos[[nm]] <- compare_dist_vs_double_zero(res) %>%
    mutate(proyecto = nm, .before = 1)
  rm(res); gc()
}

tabla_resumen_todos <- bind_rows(tabla_resumen_todos) %>%
  arrange(proyecto, metodo)

# Ver todo
print(tabla_resumen_todos, n = nrow(tabla_resumen_todos))

# Guardar si querés
# write.csv(tabla_resumen_todos, "tabla_resumen_dist_vs_doble_cero.csv", row.names = FALSE)