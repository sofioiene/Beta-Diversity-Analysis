library(dplyr)
library(ggplot2)
library(readr)
library(grid)      # para unit()
# ---- Parámetros ----
TOP_N <- 10
MODO_SELECCION <- "correlacion"   # <<--- usar correlación

# helper: llevar a relativa y alinear filas con sample_ids del PCoA
to_relative <- function(M) {
  M <- as.data.frame(M)
  M[is.na(M)] <- 0
  rs <- rowSums(M)
  keep <- rs > 0
  if (!all(keep)) M <- M[keep, , drop = FALSE]
  rs <- pmax(rowSums(M), 1e-12)
  sweep(M, 1, rs, "/")
}

out_dir <- "/Users/sofioiene/Desktop/Dimensiones Codigo/PCOA con Vectores2/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_list_biplot <- list()

for (proj in names(pcoa_scores_list)) {
  scr <- pcoa_scores_list[[proj]]               # sample.id, group, PCoA1, PCoA2
  var_exp <- pcoa_var_list[[proj]]
  if (is.null(scr) || nrow(scr) < 3) next
  
  # matriz de abundancias relativas post-NZV (muestras en filas)
  A <- abund_post_nzv_list[[proj]] %>% to_relative()
  # alinear filas con el orden de scr$sample.id
  A <- A[match(scr$sample.id, rownames(A)), , drop = FALSE]
  
  # quitar columnas con varianza 0 (correlación indefinida)
  keep_cols <- apply(A, 2, sd, na.rm = TRUE) > 0
  A <- A[, keep_cols, drop = FALSE]
  if (ncol(A) == 0) next
  
  # --- vectores de taxa por correlación con ejes ---
  cors1 <- apply(A, 2, function(x) suppressWarnings(cor(x, scr$PCoA1, use="pairwise.complete.obs")))
  cors2 <- apply(A, 2, function(x) suppressWarnings(cor(x, scr$PCoA2, use="pairwise.complete.obs")))
  otu_scores <- data.frame(
    OTU   = colnames(A),
    Axis1 = cors1,
    Axis2 = cors2,
    stringsAsFactors = FALSE
  )
  otu_scores$length <- sqrt(otu_scores$Axis1^2 + otu_scores$Axis2^2)
  
  # seleccionar Top N por correlación (longitud del vector)
  top_otus <- otu_scores %>% arrange(desc(length)) %>% slice(1:min(TOP_N, n()))
  if (nrow(top_otus) == 0) next
  
  # escalar flechas para que queden dentro del panel
  rngx <- diff(range(scr$PCoA1))
  rngy <- diff(range(scr$PCoA2))
  max_arrow <- max(sqrt(top_otus$Axis1^2 + top_otus$Axis2^2))
  scale_factor <- 0.9 * min(rngx, rngy) / (2 * max_arrow + 1e-12)  # factor conservador
  
  top_otus$xend <- top_otus$Axis1 * scale_factor
  top_otus$yend <- top_otus$Axis2 * scale_factor
  
  # etiquetas de ejes (maneja casos sin nombres en var_exp)
  ax1_lab <- if (!is.null(var_exp) && !is.na(var_exp["Axis1"])) round(var_exp["Axis1"], 2) else NA
  ax2_lab <- if (!is.null(var_exp) && !is.na(var_exp["Axis2"])) round(var_exp["Axis2"], 2) else NA
  
  p <- ggplot() +
    geom_point(data = scr, aes(x = PCoA1, y = PCoA2, color = group), size = 2.5) +
    geom_segment(data = top_otus,
                 aes(x = 0, y = 0, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.4, color = "gray40") +
    geom_text(data = top_otus,
              aes(x = xend, y = yend, label = OTU),
              size = 3, vjust = -0.4, color = "black") +
    labs(
      title = paste0("Proyecto: ", proj),
      subtitle = paste0("PCoA Bray–Curtis (post nearZeroVar) – Top ", TOP_N, " por correlación"),
      x = if (!is.na(ax1_lab)) paste0("PCoA1 (", ax1_lab, "%)") else "PCoA1",
      y = if (!is.na(ax2_lab)) paste0("PCoA2 (", ax2_lab, "%)") else "PCoA2",
      color = "Condición"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  plot_list_biplot[[proj]] <- p
  
  ggsave(
    filename = paste0("PCoA_biplot_", proj, ".png"),
    plot     = p,
    path     = out_dir,
    width    = 8, height = 6, dpi = 300
  )
}
