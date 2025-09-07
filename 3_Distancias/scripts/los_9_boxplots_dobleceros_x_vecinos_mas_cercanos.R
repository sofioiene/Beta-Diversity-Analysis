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
library(ggpubr)


my_wd<-"/home/biolab/Dropbox/research/beta_analysis/"

archivos <- list.files(path = paste0(my_wd, "kendall_test/"), pattern = "^df_doble_zero_.*\\.tsv$", full.names = TRUE)

nombres_extraidos <- sub("^df_doble_zero_(.*)\\.tsv$", "\\1", basename(archivos))

asvtable <- read_tsv(paste0(my_wd,"kendall_test/",archivos[1]), comment = "#")

lista_df_largos <- map2(archivos, nombres_extraidos, ~ {
  read_tsv(.x) %>%
    pivot_longer(
      cols = starts_with("neighbor_"),
      names_to = "neighbor_pos",
      names_prefix = "neighbor_",
      values_to = "value"
    ) %>%
    mutate(
      neighbor_pos = as.integer(neighbor_pos),
      archivo_id = .y
    )
})

vecinos_all_dist <- bind_rows(lista_df_largos)

grupo<-40

# Crear nueva columna para agrupar en bins de 10
vecinos_all_dist <- vecinos_all_dist %>%
  mutate(grupo_vecinos = paste0(
    floor((neighbor_pos - 1) / grupo) * grupo + 1, "-", 
    ceiling((neighbor_pos) / grupo) * grupo
  ))

# Ordenar el factor para que el eje x esté bien ordenado
vecinos_all_dist$grupo_vecinos <- factor(
  vecinos_all_dist$grupo_vecinos,
  levels = unique(vecinos_all_dist$grupo_vecinos[order(as.numeric(gsub("-.*", "", vecinos_all_dist$grupo_vecinos)))])
)

# Calcular la mediana por grupo
medianas <- vecinos_all_dist %>%
  group_by(archivo_id, grupo_vecinos) %>%
  summarise(mediana = median(value, na.rm = TRUE), .groups = "drop")


ggboxplot(
  vecinos_all_dist, 
  x = "grupo_vecinos", 
  y = "value",
  outlier.shape = NA,
  color = "grupo_vecinos",
  notch = TRUE,
  facet.by = "archivo_id",
  scales = "free_y"
) +
  geom_point(
    data = medianas,
    aes(x = grupo_vecinos, y = mediana, group = 1),
    color = "steelblue", size = 1.8
  ) +
  geom_line(
    data = medianas,
    aes(x = grupo_vecinos, y = mediana, group = 1),
    color = "steelblue", linewidth = 0.8
  ) +
  labs(
    x = "Groups of neighbor proximity (in blocks of 40)",
    y = "Number of double zero",
    title = "Distribución de ceros dobles por grupo de vecinos",
    subtitle = "Agrupado por distancias (archivo_id)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = "none"
  )
