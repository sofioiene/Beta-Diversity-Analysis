# Lista para resultados
tabla_resumen_todos1 <- list()

# Loop seguro
for(nm in names(lista_counts)){
  message("Procesando: ", nm)
  M <- lista_counts[[nm]]
  
  # Submuestreo defensivo
  if(ncol(M) > 400){
    set.seed(123)
    keep <- sample(colnames(M), 400)
    M <- M[, keep, drop = FALSE]
  }
  
  # Analiza y guarda
  res <- analyze_one_fast(M)  # usa la versión rápida con crossprod()
  tabla_resumen_todos1[[nm]] <- compare_dist_vs_double_zero(res) %>%
    mutate(proyecto = nm, .before = 1)
  
  rm(res); gc()
}

# Unir todo
tabla_resumen_todos1 <- bind_rows(tabla_resumen_todos1) %>%
  arrange(proyecto, metodo)

# Ver resumen
print(tabla_resumen_todos1, n = nrow(tabla_resumen_todos1))




library(ggplot2)

ggplot(tabela <- tabla_resumen_todos1,
       aes(x = metodo, y = proyecto, fill = kendall_tau)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  labs(title="Kendall τ: distancia vs % doble cero",
       x="", y="", fill="τ") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle=45, hjust=1))



library(dplyr); library(tidyr)

resumen_metodo <- tabla_resumen_todos1 %>%
  group_by(metodo) %>%
  summarise(
    n_proj = n(),
    med_tau = median(kendall_tau, na.rm=TRUE),
    IQR_tau = IQR(kendall_tau, na.rm=TRUE),
    med_mantel = median(mantel_r, na.rm=TRUE),
    IQR_mantel = IQR(mantel_r, na.rm=TRUE)
  ) %>%
  arrange(desc(med_tau))  # cambia a gusto

print(resumen_metodo)




library(purrr); library(tibble); library(dplyr)

pseudos <- c(0.5, 1, 5)
tabla_pseudo <- list()

for (ps in pseudos){
  tmp <- lapply(names(lista_counts), function(nm){
    M <- lista_counts[[nm]]
    if(ncol(M) > 400){ set.seed(123); M <- M[, sample(colnames(M), 400)] }
    dists <- list(aitchison_clr = as.matrix(dist(t(norm_clr(M, pseudo = ps)))))
    Z <- double_zero_matrix_fast(M)
    U <- upper.tri(Z)
    D <- dists$aitchison_clr
    tibble(
      proyecto = nm,
      metodo = paste0("aitchison_clr_p", ps),
      kendall_tau = suppressWarnings(cor(D[U], Z[U], method="kendall")),
      mantel_r = tryCatch(suppressWarnings(mantel(D, Z, method="pearson", permutations=0)$statistic), error=function(e) NA_real_)
    )
  }) %>% bind_rows()
  tabla_pseudo[[as.character(ps)]] <- tmp
  rm(tmp); gc()
}
tabla_pseudo <- bind_rows(tabla_pseudo)
print(tabla_pseudo, n = nrow(tabla_pseudo))
