asvtable <- read_tsv(paste0(my_wd,"2. files tsv/106_infant_singapore.tsv"), comment = "#", skip=1)

archivos <- list.files(path = paste0(my_wd, "kendall_test/"), pattern = "^kendall_.*106_infant_singapore\\.tsv$", full.names = TRUE)

nombres_extraidos <- sub("^kendall_(.*)106_infant_singapore\\.tsv$", "\\1", basename(archivos))


# definimos la línea de corte entre discordantes y concordantes: 
threshold<-0

asv_matriz <- asvtable %>%
  select(where(is.numeric))  # mantiene solo columnas numéricas = samples

totales_kendall <- map2_dfr(archivos, nombres_extraidos, ~ {
  read_tsv(.x, col_types = cols(
    sample = col_character(),
    value = col_double()
  )) %>%
    mutate(
      distances = .y,
      dist_1 = sub("_vs_.*", "", .y),
      dist_2 = sub(".*_vs_", "", .y)
    )
})

totales_kendall <- totales_kendall %>%
  mutate(grupo_concordancia = case_when(
    value < 0 ~ "Discordant",
    value >= 0 & value < 0.3 ~ "Low concordance",
    value >= 0.3 & value < 0.6 ~ "Moderate concordance",
    value >= 0.6 ~ "High concordance"
  )) %>%
  mutate(grupo_concordancia = factor(grupo_concordancia, levels = c(
    "Discordant", "Low concordance", "Moderate concordance", "High concordance"
  )))

# 2. Distancias disponibles
todas_las_distancias <- unique(totales_kendall$distances)

# 3. Lista para guardar resultados
lista_resultados <- list()

for (dist in todas_las_distancias) {
  
  mypru <- totales_kendall %>% filter(distances == dist)
  cambios_fuertes <- mypru %>% filter(value < threshold) %>% pull(sample)
  
  dist_1 <- unique(mypru$dist_1)
  dist_2 <- unique(mypru$dist_2)
  
  asv_cambia <- asv_matriz[, colnames(asv_matriz) %in% cambios_fuertes, drop = FALSE]
  asv_nocambia <- asv_matriz[, !(colnames(asv_matriz) %in% cambios_fuertes), drop = FALSE]
  
  if (length(cambios_fuertes) == 0) {
    pruCambia <- data.frame(
      porcentaje = rep(NA, ncol(asv_nocambia)),
      categ = "Change",
      distancia = dist,
      dist_1 = dist_1,
      dist_2 = dist_2
    )
  } else {
    pruCambia <- data.frame(
      porcentaje = colMeans(asv_cambia == 0),
      categ = "Change",
      distancia = dist,
      dist_1 = dist_1,
      dist_2 = dist_2
    )
  }
  
  pruNoCambia <- data.frame(
    porcentaje = colMeans(asv_nocambia == 0),
    categ = "No change",
    distancia = dist,
    dist_1 = dist_1,
    dist_2 = dist_2
  )
  
  lista_resultados[[dist]] <- rbind(pruCambia, pruNoCambia)
}

# 8. Combinar todo
datos_completos <- bind_rows(lista_resultados)

# 9. Etiquetas más legibles
datos_completos$categ <- factor(datos_completos$categ, 
                                levels = c("Change", "No change"),
                                labels = c("Discordant", "Concordant"))

orden_columnas <- c("gower", "raup", "manhattan_raw", "cao", "chao")
orden_filas <- c("bray_edger", "chao", "cao", "manhattan_raw", "raup")

# Asegurar que los niveles estén correctamente ordenados
datos_completos <- datos_completos %>%
  mutate(
    dist_1 = factor(dist_1, levels = orden_columnas),
    dist_2 = factor(dist_2, levels = rev(orden_filas))
  )

# Crear tabla con combinaciones reales de dist_1 y dist_2
#combinaciones <- datos_completos %>%
#  dplyr::distinct(dist_1, dist_2) %>%
#  dplyr::mutate(x = 0.3, xend = 0.3, y = 0, yend = Inf)

ggplot(datos_completos, aes(x = porcentaje, fill = categ)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(values = c("Discordant" = "red", "Concordant" = "blue")) +
  labs(
    title = "Percentage of zeros per sample",
    subtitle = "Grouped by distance metric and tau threshold",
    x = "Percentage of zeros per sample",
    y = "Frequency",
    fill = "Sample group"
  ) +
  facet_grid(dist_1 ~ dist_2, scales = "free_y") +
  theme_minimal() +  theme(
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 7)
  ) 
