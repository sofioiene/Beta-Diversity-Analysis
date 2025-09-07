library(ggplot2)
library(patchwork)

plot_list_pcoa <- list()

for (proj in names(pcoa_scores_list)) {
  df <- pcoa_scores_list[[proj]]
  var_exp <- pcoa_var_list[[proj]]
  
  # por si algún proyecto quedó sin tabla (skip en el loop anterior)
  if (is.null(df)) next
  
  p <- ggplot(df, aes(x = PCoA1, y = PCoA2, color = group)) +
    geom_point(size = 2.5) +
    labs(
      title = proj,
      subtitle = "PCoA (Bray-Curtis, post nearZeroVar)",
      x = paste0("PCoA1 (", round(var_exp["Axis1"], 2), "%)"),
      y = paste0("PCoA2 (", round(var_exp["Axis2"], 2), "%)"),
      color = "Grupo"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  plot_list_pcoa[[proj]] <- p
}

# 👉 grilla general (ejemplo 4x4)
grilla_pcoa <- wrap_plots(plotlist = plot_list_pcoa, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Mostrar en R
grilla_pcoa

# 👉 Guardar cada gráfico por separado
dir.create("", showWarnings = FALSE)
for (proj in names(plot_list_pcoa)) {
  ggsave(
    filename = paste0("PCoA_", proj, ".png"),
    plot     = plot_list_pcoa[[proj]],
    path     = "/Users/sofioiene/Desktop/Dimensiones Codigo/",
    width    = 7, height = 5, dpi = 300
  )
}
