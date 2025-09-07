# Paquetes
library(vegan)
library(dplyr)
library(tibble)
library(readr)   # opcional si querés exportar CSV

# =======================
# Helpers
# =======================

to_numeric_df <- function(df) {
  df <- as.data.frame(df)
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  df[is.na(df)] <- 0
  df
}

# Construye matriz ALL: filas = muestras, columnas = OTUs (sin NZV)
# original: data.frame con col "OTU ID" + columnas de samples
# manifest: data.frame con columnas sample.id y group
build_ALL <- function(original, manifest) {
  stopifnot("sample.id" %in% names(manifest), "group" %in% names(manifest))
  
  mat <- as.data.frame(original)
  # mover OTU IDs a rownames y quitar la columna
  otu_col <- which(names(mat) %in% c("OTU ID", "OTU_ID", "#OTU ID", "#OTU_ID"))
  if (length(otu_col) == 1) {
    rownames(mat) <- mat[[otu_col]]
    mat <- mat[, -otu_col, drop = FALSE]
  } else {
    stop("No encuentro columna 'OTU ID' en 'original'.")
  }
  
  mat <- to_numeric_df(mat)         # columnas = muestras
  A <- t(mat)                       # filas = muestras
  # homogeneizar IDs: '-' vs '.'
  rownames(A) <- gsub("-", ".", rownames(A))
  
  manifest <- manifest %>%
    mutate(sample.id = gsub("-", ".", sample.id)) %>%
    arrange(sample.id)
  
  # Subset y reordenar por sample.id del manifest
  keep_ids <- intersect(manifest$sample.id, rownames(A))
  A <- A[keep_ids, , drop = FALSE]
  rownames(A) <- keep_ids
  
  # devolver matriz y vector de grupos alineados
  grp <- droplevels(as.factor(manifest$group[match(rownames(A), manifest$sample.id)]))
  list(A = A, grp = grp)
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

# =======================
# Loop ALL-OTUs por proyecto
# =======================

perm <- 9999
res <- list(); i <- 0

for (proj in names(original_list)) {
  message("\n=== Proyecto: ", proj, " (ALL OTUs) ===")
  
  original <- original_list[[proj]]
  manifest <- manifest_list[[proj]]
  
  if (is.null(original) || is.null(manifest)) {
    message("   -> faltan 'original' o 'manifest'; skip")
    next
  }
  
  built <- try(build_ALL(original, manifest), silent = TRUE)
  if (inherits(built, "try-error")) {
    message("   -> build_ALL falló: ", built)
    next
  }
  
  A <- built$A
  grp <- built$grp
  
  # eliminar filas con suma 0 y grupos NA
  rs <- rowSums(A, na.rm = TRUE)
  keep <- rs > 0 & !is.na(grp)
  dropped <- sum(!keep)
  A <- A[keep, , drop = FALSE]
  grp <- droplevels(grp[keep])
  
  message("   -> muestras: ", nrow(A), " | OTUs: ", ncol(A), " | eliminadas: ", dropped)
  
  if (nrow(A) < 3 || length(unique(grp)) < 2) {
    i <- i + 1
    res[[i]] <- tibble(
      proyecto = proj,
      n_muestras = nrow(A),
      grupos = paste(names(table(grp)), table(grp), collapse = "; "),
      R2 = NA_real_, p_permanova = NA_real_, p_dispersion = NA_real_,
      nota = "menos de 3 muestras o <2 grupos"
    )
    next
  }
  
  # Distancia Bray–Curtis en ALL OTUs (si NA, fallback a TSS)
  d <- try(vegdist(A, method = "bray"), silent = TRUE)
  if (inherits(d, "try-error") || any(is.na(d))) {
    A_tss <- sweep(A, 1, pmax(rowSums(A), 1e-12), "/")
    d <- try(vegdist(A_tss, method = "bray"), silent = TRUE)
  }
  if (inherits(d, "try-error") || any(is.na(d))) {
    message("   -> distancias con NA; skip")
    i <- i + 1
    res[[i]] <- tibble(
      proyecto = proj,
      n_muestras = nrow(A),
      grupos = paste(names(table(grp)), table(grp), collapse = "; "),
      R2 = NA_real_, p_permanova = NA_real_, p_dispersion = NA_real_,
      nota = "distancias con NA"
    )
    next
  }
  
  per <- safe_permanova(d, grp, perm = perm)
  bd  <- safe_betadisp(d, grp, perm = perm)
  
  message(sprintf("   -> PERMANOVA: R2=%s  p=%s",
                  ifelse(is.na(per$R2), "NA", format(per$R2, digits = 3)),
                  ifelse(is.na(per$p),  "NA", format(per$p,  digits = 3))))
  message(sprintf("   -> betadisp p=%s",
                  ifelse(is.na(bd$p), "NA", format(bd$p, digits = 3))))
  
  i <- i + 1
  res[[i]] <- tibble(
    proyecto = proj,
    n_muestras = nrow(A),
    grupos = paste(names(table(grp)), table(grp), collapse = "; "),
    R2 = per$R2, p_permanova = per$p, p_dispersion = bd$p,
    nota = dplyr::coalesce(per$note, bd$note, NA_character_)
  )
}

permanova_ALL_summary <- bind_rows(res) %>% arrange(proyecto)
print(permanova_ALL_summary)

# (opcional) exportar:
# write_csv(permanova_ALL_summary, "PERMANOVA_ALL_OTUs_por_proyecto.csv")