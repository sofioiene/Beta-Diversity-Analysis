library(reshape2)
library(gridExtra)
library(tidyverse)
#library(ape)
library(magrittr)
library(plotly)
library(kableExtra)
#library(Polychrome)
library(ggplot2)
library(ggpmisc)
library(vegan)
library(cluster)
library(DescTools)

my_wd<-"/home/biolab/Dropbox/research/beta_analysis/"

archivos <- list.files(path = paste0(my_wd, "kendall_test/"), pattern = "^kendall_.*106_infant_singapore\\.tsv$", full.names = TRUE)

nombres_extraidos <- sub("^kendall_(.*)106_infant_singapore\\.tsv$", "\\1", basename(archivos))

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


orden_columnas <- c("gower", "raup", "manhattan_raw", "cao", "chao")
orden_filas <- c("bray_edger", "chao", "cao", "manhattan_raw", "raup")

totales_kendall <- totales_kendall %>%
  mutate(
    dist_1 = factor(dist_1, levels = orden_columnas),
    dist_2 = factor(dist_2, levels = rev(orden_filas))
  )

# Crear tabla con combinaciones reales de dist_1 y dist_2
combinaciones <- totales_kendall %>%
  dplyr::distinct(dist_1, dist_2) %>%
  dplyr::mutate(x = 0, xend = 0, y = 0, yend = Inf)

# Graficar con líneas solo en las combinaciones reales
ggplot(totales_kendall, aes(x = value, fill = grupo_concordancia)) +
  geom_histogram(binwidth = 0.05, color = "black") +
  facet_grid(dist_1 ~ dist_2, drop = TRUE) +
  scale_fill_manual(values = c(
    "Discordant" = "red",
    "Low concordance" = "orange",
    "Moderate concordance" = "skyblue",
    "High concordance" = "forestgreen"
  )) +
  labs(
    title = "Distribución de coeficientes de Kendall",
    subtitle = "Faceteado por combinaciones de distancias",
    x = "Kendall Tau",
    y = "Frequency",
    fill = "Concordance level"
  ) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  geom_segment(
    data = combinaciones,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  )

