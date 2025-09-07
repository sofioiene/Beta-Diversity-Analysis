my_wd <- "/home/biolab/Dropbox/research/beta_analysis/"

mantel_dir <- file.path(my_wd, "")

# Listar archivos que empiezan con "dist" y terminan con ".tsv"
rds_files <- list.files(paste0(my_wd, "mantel_test/"), pattern = "^dist.*\\.tsv$")

# Inicializamos lista para guardar los data.frames
all_results <- list()

for (myfile in rds_files) {
  # Leer archivo
  df <- read_tsv(paste0(my_wd, "mantel_test/", myfile))
  
  # Extraer nombre de la distancia desde el nombre del archivo
  dist_name <- myfile %>% 
    basename() %>% 
    str_remove("^dist_") %>% 
    str_remove("\\.tsv$")
  
  # Agregar columna con el nombre de la distancia
  df$distance <- dist_name
  
  # Guardar en lista
  all_results[[dist_name]] <- df
}

# Unir todos los data.frames en uno solo
mantel_combined <- bind_rows(all_results)

mantel_combined[mantel_combined$distance=="binomial", "distance"]<-"manhattan_log"

# Calcular resumen por distancia y threshold
summary_mantel <- mantel_combined %>%
  group_by(distance, threshold) %>%
  summarise(
    mean_mantel = mean(mantel_r_comunes, na.rm = TRUE),
    sd_mantel = sd(mantel_r_comunes, na.rm = TRUE),
    mean_nceros = mean(nceros, na.rm = TRUE),
    sd_nceros = sd(nceros, na.rm = TRUE),
    mean_taxacomunes = mean(ntaxacomunes, na.rm = TRUE),
    sd_taxacomunes = sd(ntaxacomunes, na.rm = TRUE),
    .groups = "drop"
  )

# Escalar dentro de cada distancia
summary_mantel <- summary_mantel %>%
  group_by(distance) %>%
  mutate(scale_factor = max(mean_mantel, na.rm = TRUE) / max(mean_nceros, na.rm = TRUE)) %>%
  ungroup()

# Graficar
p <- ggplot(summary_mantel, aes(x = threshold)) +
  # Línea y ribbon del test de Mantel
  geom_line(aes(y = mean_mantel, color = "Test de Mantel"), size = 1) +
  geom_ribbon(aes(ymin = mean_mantel - sd_mantel,
                  ymax = mean_mantel + sd_mantel,
                  fill = "Test de Mantel"),
              alpha = 0.2) +
  
  # Línea y ribbon de % de ceros (escalado)
  geom_line(aes(y = mean_nceros * scale_factor, color = "% de ceros"), linewidth = 1) +
  geom_ribbon(aes(ymin = (mean_nceros - sd_nceros) * scale_factor,
                  ymax = (mean_nceros + sd_nceros) * scale_factor,
                  fill = "% de ceros"),
              alpha = 0.2) +
  
  # Línea horizontal fija
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red", size = 0.5) +
  
  # Escala de ejes
  scale_y_continuous(
    name = "Media(Estadístico r del test de Mantel)",
    sec.axis = sec_axis(~ . / max(summary_mantel$scale_factor), name = "Porcentaje de ceros (%)")
  ) +
  scale_color_manual(
    name = "Curvas",
    values = c("Test de Mantel" = "red", "% de ceros" = "blue")
  ) +
  scale_fill_manual(
    name = "Curvas",
    values = c("Test de Mantel" = "red", "% de ceros" = "blue")
  ) +
  labs(
    title = "Comparación entre distancia y dispersión de ceros",
    subtitle = "Promedio del test de Mantel (comunes) vs. porcentaje de ceros",
    x = "Threshold de prevalencia"
  ) +
  facet_wrap(~ distance, nrow = 2, ncol = 4, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "red"),
    axis.title.y.right = element_text(color = "blue"),
    legend.position = "top"
  )

p

ggsave(
  filename = paste0(my_wd, "mantel_test/", "mantel_plot_raros.jpg"),
  plot = p,
  width = 16,       # en pulgadas
  height = 8,       # en pulgadas
  dpi = 300         # alta resolución
)


# 2do grafico 
# graficamos los raros vs los comunes! 

# Resumen por distancia y threshold
summary_mantel <- mantel_combined %>%
  group_by(distance, threshold) %>%
  summarise(
    mean_mantel_comunes = mean(mantel_r_comunes, na.rm = TRUE),
    sd_mantel_comunes = sd(mantel_r_comunes, na.rm = TRUE),
    mean_mantel_raros = mean(mantel_r_raros, na.rm = TRUE),
    sd_mantel_raros = sd(mantel_r_raros, na.rm = TRUE),
    .groups = "drop"
  )

# Graficar ambas curvas en el mismo eje Y
p <- ggplot(summary_mantel, aes(x = threshold)) +
  # Test de Mantel (comunes)
  geom_line(aes(y = mean_mantel_comunes, color = "Common"), size = 1) +
  geom_ribbon(aes(
    ymin = mean_mantel_comunes - sd_mantel_comunes,
    ymax = mean_mantel_comunes + sd_mantel_comunes,
    fill = "Common"
  ), alpha = 0.2) +
  
  # Test de Mantel (raros)
  geom_line(aes(y = mean_mantel_raros, color = "Low-prevalence"), size = 1) +
  geom_ribbon(aes(
    ymin = mean_mantel_raros - sd_mantel_raros,
    ymax = mean_mantel_raros + sd_mantel_raros,
    fill = "Low-prevalence"
  ), alpha = 0.2) +
  
  # Línea de referencia horizontal
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red", size = 0.5) +
  scale_y_continuous(name = "r statistic Mantel test (mean)") +
  scale_color_manual(
    name = "Taxa type",
    values = c("Common" = "darkred", "Low-prevalence" = "darkblue")
  ) +
  scale_fill_manual(
    name = "Taxa type",
    values = c("Common" = "red", "Low-prevalence" = "blue")
  ) +
  labs(
    title = "Test de Mantel: comunes vs. raros",
    subtitle = "Promedio ± desviación estándar por distancia y threshold",
    x = "Prevalence threshold"
  ) +
  facet_wrap(~ distance, nrow = 2, ncol = 4, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(color = "black"),
    legend.position = "top"
  )

p

## 3er grafico 

# Función auxiliar para encontrar intersección con y = 0.75
get_threshold_at_075 <- function(df, column_name) {
  df %>%
    arrange(threshold) %>%
    group_by(distance) %>%
    summarise(
      cross = list({
        x = threshold
        y = .data[[column_name]]
        idx = which(diff(sign(y - 0.75)) != 0)
        if (length(idx) > 0) {
          # interpolación lineal
          i = idx[1]
          x0 = x[i]; x1 = x[i+1]
          y0 = y[i]; y1 = y[i+1]
          x_int = x0 + (0.75 - y0) * (x1 - x0) / (y1 - y0)
          tibble(threshold = x_int)
        } else {
          tibble()
        }
      }),
      .groups = "drop"
    ) %>%
    unnest(cross)
}

# Calcular intersecciones
cross_common <- get_threshold_at_075(summary_mantel, "mean_mantel_comunes") %>%
  mutate(type = "Common")

cross_raros <- get_threshold_at_075(summary_mantel, "mean_mantel_raros") %>%
  mutate(type = "Low-prevalence")

cross_labels <- cross_labels %>%
  mutate(
    x_label = threshold + ifelse(type == "Common", -0.07, 0.07),
    y_label = 0.75 - 0.4
  )

# Graficar
p <- ggplot(summary_mantel, aes(x = threshold)) +
  # Líneas y áreas
  geom_line(aes(y = mean_mantel_comunes, color = "Common"), linewidth = 1) +
  geom_ribbon(aes(ymin = mean_mantel_comunes - sd_mantel_comunes,
                  ymax = mean_mantel_comunes + sd_mantel_comunes,
                  fill = "Common"), alpha = 0.2) +
  
  geom_line(aes(y = mean_mantel_raros, color = "Low-prevalence"), linewidth = 1) +
  geom_ribbon(aes(ymin = mean_mantel_raros - sd_mantel_raros,
                  ymax = mean_mantel_raros + sd_mantel_raros,
                  fill = "Low-prevalence"), alpha = 0.2) +
  
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "black") +
  
  # Puntos de intersección
  geom_point(data = cross_labels, aes(x = threshold, y = 0.75, color = type), size = 2) +
  
  # Segmentos que conectan punto con texto
  geom_segment(data = cross_labels,
               aes(x = threshold, xend = x_label, y = 0.75, yend = y_label, color = type),
               arrow = arrow(length = unit(0.02, "npc")), linewidth = 0.4, show.legend = FALSE) +
  
  # Texto desplazado
  geom_label(data = cross_labels,
             aes(x = x_label, y = y_label, label = round(threshold, 2)),
             fill = "lightgrey",             # fondo gris claro
             color = "black",                # texto negro
             label.size = 0,
             alpha = 0.6,
             size = 3.5, show.legend = FALSE) +  # Escalas
  scale_color_manual(
    name = "Taxa type",
    values = c("Common" = "darkred", "Low-prevalence" = "darkblue")
  ) +
  scale_fill_manual(
    name = "Taxa type",
    values = c("Common" = "red", "Low-prevalence" = "blue")
  ) +
  
  # Etiquetas y tema
  labs(
    title = "Mantel test by taxa prevalence",
    subtitle = "Comparison of common vs. low-prevalence taxa",
    x = "Prevalence threshold",
    y = "Mean r statistic test Mantel"
  ) +
  facet_wrap(~ distance, nrow = 2, ncol = 4, scales = "free_y") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title.y = element_text(color = "black")
  )

# Mostrar gráfico
p
