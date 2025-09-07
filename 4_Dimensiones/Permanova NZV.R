library(vegan)
library(dplyr)
library(tibble)

# --- Helpers -------------------------------

to_numeric_df <- function(df) {
  df <- as.data.frame(df)
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  df[is.na(df)] <- 0
  df
}

align_groups <- function(X, y) {
  # y puede venir como vector sin nombres o data.frame
  if (is.data.frame(y)) {
    if ("group" %in% names(y)) y <- y$group else y <- unlist(y, use.names = FALSE)
  }
  y <- as.vector(y)
  
  # Si no tiene nombres pero el largo coincide, asumo el orden de X
  if (is.null(names(y)) || all(is.na(names(y)))) {
    if (length(y) == nrow(X)) {
      names(y) <- rownames(X)
    } else {
      stop("El vector de grupos no tiene nombres y su largo no coincide con nrow(X).")
    }
  }
  
  # Homogeneizar '-' vs '.'
  rn <- rownames(X)
  ny <- names(y)
  
  rn_dot <- gsub("-", ".", rn)
  ny_dot <- gsub("-", ".", ny)
  
  # Intenta alinear por nombres (con puntos)
  names(y) <- ny_dot
  rownames(X) <- rn_dot
  
  # Subset & reordenar
  y <- y[ rownames(X) ]
  
  list(X = X, y = factor(y))
}

safe_permanova <- function(dist, grp, perm = 9999) {
  out <- list(R2 = NA_real_, p = NA_real_, note = NA_character_)
  if (length(grp) < 3 || length(unique(grp)) < 2) {
    out$note <- "menos de 3 muestras o <2 grupos"
    return(out)
  }
  ad <- try(adonis2(dist ~ grp, permutations = perm), silent = TRUE)
  if (inherits(ad, "try-error")) { out$note <- "adonis2 error"; return(out) }
  out$R2 <- as.numeric(ad$R2[1])
  out$p  <- as.numeric(ad$`Pr(>F)`[1])
  out
}

safe_betadisp <- function(dist, grp, perm = 9999) {
  out <- list(p = NA_real_, note = NA_character_)
  if (length(grp) < 3 || length(unique(grp)) < 2) {
    out$note <- "menos de 3 muestras o <2 grupos"
    return(out)
  }
  bd <- try(betadisper(dist, grp), silent = TRUE)
  if (inherits(bd, "try-error")) { out$note <- "betadisper error"; return(out) }
  pt <- try(permutest(bd, permutations = perm), silent = TRUE)
  if (inherits(pt, "try-error")) { out$note <- "permutest error"; return(out) }
  out$p <- as.numeric(pt$tab[1, "Pr(>F)"])
  out
}

# --- Loop principal ------------------------

res <- list(); k <- 0

for (proj in names(abund_post_nzv_list)) {
  message("\n=== Proyecto: ", proj, " ===")
  
  X0 <- abund_post_nzv_list[[proj]]
  y0 <- groups_list[[proj]]
  
  if (is.null(X0) || is.null(y0)) {
    message("   -> faltan X o grupos; skip")
    next
  }
  
  # numérico y saneado
  X <- to_numeric_df(X0)
  
  # alinear grupos aun si y no tiene names()
  al <- try(align_groups(X, y0), silent = TRUE)
  if (inherits(al, "try-error")) {
    message("   -> no pude alinear grupos: ", al)
    next
  }
  X <- al$X; y <- al$y
  
  # eliminar filas con suma 0 o y NA
  rs <- rowSums(X, na.rm = TRUE)
  keep <- rs > 0 & !is.na(y)
  dropped <- sum(!keep)
  X <- X[keep, , drop = FALSE]
  y <- droplevels(y[keep])
  
  message("   -> muestras: ", nrow(X), " | OTUs: ", ncol(X), " | eliminadas: ", dropped)
  
  if (nrow(X) < 3 || length(unique(y)) < 2) {
    k <- k + 1
    res[[k]] <- tibble(
      proyecto = proj, n_muestras = nrow(X),
      grupos = paste(names(table(y)), table(y), collapse = "; "),
      R2 = NA_real_, p_permanova = NA_real_, p_dispersion = NA_real_,
      nota = "menos de 3 muestras o <2 grupos"
    )
    next
  }
  
  # Distancias Bray–Curtis (si da NA, probamos TSS rápido)
  d <- try(vegdist(X, method = "bray"), silent = TRUE)
  if (inherits(d, "try-error") || any(is.na(d))) {
    X_tss <- sweep(X, 1, pmax(rowSums(X), 1e-12), "/")
    d <- try(vegdist(X_tss, method = "bray"), silent = TRUE)
  }
  if (inherits(d, "try-error") || any(is.na(d))) {
    message("   -> distancias con NA; skip")
    k <- k + 1
    res[[k]] <- tibble(
      proyecto = proj, n_muestras = nrow(X),
      grupos = paste(names(table(y)), table(y), collapse = "; "),
      R2 = NA_real_, p_permanova = NA_real_, p_dispersion = NA_real_,
      nota = "distancias con NA"
    )
    next
  }
  
  per <- safe_permanova(d, y, perm = 9999)
  bd  <- safe_betadisp(d, y, perm = 9999)
  
  message(sprintf("   -> PERMANOVA: R2=%s  p=%s",
                  ifelse(is.na(per$R2), "NA", format(per$R2, digits=3)),
                  ifelse(is.na(per$p),  "NA", format(per$p,  digits=3))))
  message(sprintf("   -> betadisp p=%s",
                  ifelse(is.na(bd$p), "NA", format(bd$p, digits=3))))
  
  k <- k + 1
  res[[k]] <- tibble(
    proyecto = proj,
    n_muestras = nrow(X),
    grupos = paste(names(table(y)), table(y), collapse = "; "),
    R2 = per$R2, p_permanova = per$p, p_dispersion = bd$p,
    nota = dplyr::coalesce(per$note, bd$note, NA_character_)
  )
}

permanova_summary <- bind_rows(res) %>% arrange(proyecto)
print(permanova_summary)
# write.csv(permanova_summary, "permanova_betadisp_resumen.csv", row.names = FALSE)