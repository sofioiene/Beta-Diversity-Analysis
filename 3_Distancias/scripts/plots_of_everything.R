library(vegan)
library(tidyverse)
library(magrittr)
library(kableExtra)
library(gridExtra)
library(diverse)
library(reshape2)
library(ggpmisc)  
library(igraph)
library(broom)

my_wd<-"/Users/sofioiene/Desktop/double_zero-main/scripts"

rds_files<-list.files(paste0(my_wd), pattern = "tsv")

i=0
for(RDS_FILE_NAME in rds_files){ 
  myrds_file<-RDS_FILE_NAME
  mifile <- read_tsv(paste0(my_wd,myrds_file), comment = "#") 
  
  filemelted<- mifile %>%
    pivot_longer(
      cols = -c(Var1, Var2, dzero, grupoVar1, grupoVar2),  # dejamos estas columnas quietas
      names_to = "distancia",
      values_to = "valor"
    )
  
  filemelted$project<-str_remove(RDS_FILE_NAME, ".tsv")
  
  if(i==0) {
    totales<-filemelted
    i=1
  } else {
    totales <- rbind(totales, filemelted)
  }
  totales=totales[!is.na(totales$grupoVar1),]
  totales=totales[!is.na(totales$grupoVar2),]
  
}

## grafico 1

## chequeamos correlaciones entre todas las distancias: 
## graficar las correlaciones entre distancias: 

totales_wide <- totales %>%
  pivot_wider(
    id_cols = c(Var1, Var2, grupoVar1, grupoVar2, project, dzero),
    names_from = distancia,
    values_from = valor
  )

matriz_datos <- totales_wide[, 6:31]

correlaciones <- cor(matriz_datos, use = "pairwise.complete.obs")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(matriz_datos)

corrplot::corrplot(correlaciones[-1,-1],
                   diag = FALSE,
                   method = "color",
                   type = "upper",
                   tl.col = "black",
                   addCoef.col = "black",  # añade coeficientes de correlación
                   p.mat = p.mat[-1,-1],          # añade nivel de significancia
                   sig.level = 0.05,       # umbral de significancia
                   insig = "blank",        # oculta las no significativas
                   number.cex = 0.7)       # tamaño de los números

# fin grafico de correlaciones 

# vamos a hacer los grupos de distancias mayores a 0.75 


# Convertimos a matriz binaria de adyacencia: 1 si correlación > 0.75 y no es la diagonal
adj_matrix <- (correlaciones > 0.85) & upper.tri(correlaciones)
adj_matrix[is.na(adj_matrix)] <- FALSE

# Creamos un grafo no dirigido a partir de la matriz de adyacencia
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Obtenemos los componentes conectados (es decir, grupos de distancias correlacionadas)
grupos <- components(g)

# Mostramos los grupos
grupos$membership  # cada distancia con su grupo
split(names(grupos$membership), grupos$membership)  # distancias por grupo

# graficos de dispersion de las distancias seleccionadas (las 8)

# distancias elegidas que no correlacionan entre si: 
grouped<- c("chisq", "manhattan_log", "bray_edger", "raup", "mountford" ,  "manhattan_raw", "cao", "chao")

totales2<-totales[totales$distancia %in% grouped, ]

my_colors <- c(
  "af_urban_usa"       = "#5C2121",
  "asian"              = "#66CC66",
  "cameroon"           = "#0066FF",
  "children_eu_africa" = "#0000CC",
  "citizen"            = "#FF3333",
  "hadza"              = "#FFD700",
  "healthy_fibers"     = "#228B22",
  "healthy_japanese"   = "#99FF99",
  "infant"             = "#FFB6C1",
  "infant_singapore"   = "#00FFFF",
  "italian"            = "#FF99CC",
  "L_H"                = "#00CCCC",
  "lcarb_lfat"         = "#FF00FF",
  "obese"              = "#FFCC00",
  "ovese"              = "#0033CC",
  "swedish"            = "#CC66FF"
)

totales2$project <- as.factor(totales2$project)
names(my_colors) <- levels(totales2$project)

# coloreado por proyecto
grafico <- ggplot(totales2, aes(x = dzero, y = valor, color = project)) +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_wrap(~distancia, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Double Zero Frequency",
    y = "Distance value",
    title = "Double zero freq vs distance values",
    color = "Project"
  ) +
  scale_color_manual(values = my_colors) +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right"
  )

ggsave(paste0(my_wd, "dispersión_8_distancias.png"), plot = grafico, device = "png", width = 15,   # tamaño más razonable
       height = 8,   # tamaño más razonable
       units = "in", bg = "white") # usa pulgadas para el tamaño

## esta tabla muestra la correlación del grupo de distancias manhattan_log

totales_filtrado <- totales2 %>%
  filter(distancia == "manhattan_log")

tabla_modelos <- totales_filtrado %>%
  group_by(project) %>%
  summarise(
    modelo = list(lm(valor ~ dzero)),
    cor_pearson = cor(dzero, valor, method = "pearson"),
    cor_spearman = cor(dzero, valor, method = "spearman"),
    tidy_model = list(tidy(modelo[[1]])),
    glance_model = list(glance(modelo[[1]])),
    .groups = "drop"
  ) %>%
  mutate(
    b0 = map_dbl(tidy_model, ~ .x %>% filter(term == "(Intercept)") %>% pull(estimate)),
    b1 = map_dbl(tidy_model, ~ .x %>% filter(term == "dzero") %>% pull(estimate)),
    r2 = map_dbl(glance_model, "r.squared")
  ) %>%
  select(project, b0, b1, r2, cor_pearson, cor_spearman)


n <- nrow(tabla_modelos)
y_max <- max(totales_filtrado$valor)
y_min <- min(totales_filtrado$valor)
espaciado <- (y_max - y_min) / (n + 1)

tabla_modelos_ordenada <- tabla_modelos %>%
  mutate(
    y_pos = y_max - row_number() * espaciado,
    dzero = min(totales_filtrado$dzero),  # x fijo a la izquierda
    label = paste0("rho = ", round(cor_spearman, 2))
  )

grafico2 <- ggplot(totales_filtrado, aes(x = dzero, y = valor, color = project)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_text(data = tabla_modelos_ordenada,
            aes(x = dzero, y = y_pos, label = label, color = project),
            hjust = -0.05,
            size = 3,
            show.legend = FALSE) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Double Zero Frequency",
    y = "Distance value",
    title = "Double zero freq vs Manhattan distance",
    color = "Project"
  ) +
  scale_color_manual(values = my_colors) +
  theme(
    legend.position = "right"
  ) +
  xlim(min(totales_filtrado$dzero) - 0.03, NA)  # para que no se corte

grafico2

# tabla de correlaciones de las distancias seleccionadas: 
# por proyecto
# (no tiene sentido)

cor_matrix <- totales2 %>%
  group_by(project, distancia) %>%                     # agrupás por proyecto y distancia
  summarize(correlacion = cor(dzero, valor, method = "pearson"), .groups = "drop") %>%  # calculás correlación
  pivot_wider(
    names_from = distancia, 
    values_from = correlacion
  )

cor_matrix_prom <- cor_matrix %>%
  bind_rows(
    cor_matrix %>%
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
      mutate(project = "Average")  # nombre de la fila de promedio
  )

formatted_table <- cor_matrix_prom %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  mutate(across(where(is.numeric), ~cell_spec(
    .x,
    color = ifelse(abs(.x) > 0.75, "white", "black"),
    background = ifelse(abs(.x) > 0.75, "red", "")
  ))) %>%
  mutate(project = ifelse(project == "Average",
                          cell_spec(project, bold = TRUE),
                          project)) %>%
  kable(booktabs = TRUE, linesep = "", format = "html", escape = FALSE) %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(nrow(cor_matrix_prom), bold = TRUE, extra_css = "border-top: double;")

formatted_table

### boxplot que discrimina entre cohortes
###

totales2$same_group <- totales2$grupoVar1 == totales2$grupoVar2

p <- ggplot(totales2, aes(x = same_group, y = valor)) +
  geom_jitter(
    aes(color = dzero),  # ahora el color solo es para los puntos
    size = 0.4,
    width = 0.2
  ) +
  geom_boxplot(
    aes(fill = same_group),  # ¡usar fill acá!
    color = "black",          # borde negro fijo
    alpha = 0.3, 
    outlier.shape = NA,
    notch = TRUE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "yellow", "FALSE" = "red"), 
    labels = c("TRUE" = "SAME", "FALSE" = "DIFF")
  ) +
  scale_color_gradient(
    low = "#004d00",   # verde clarito (inicio de la escala)
    high = "#E0ECD2"   # verde oscuro (final de la escala)
  ) +
  scale_x_discrete(
    labels = NULL
  ) +
  labs(
    color = "Double zero \n frequency", 
    fill = "Sample cohort", 
    title = "All projects"  # Título agregado aquí
  ) +
  theme_minimal() +
  facet_wrap(~distancia, scales = "free_y", ncol = 8) +
  xlab("Distances") +  # Cambiar el nombre del eje X
  ylab("Value")  # Cambiar el nombre del eje Y

pp<-p + stat_compare_means(aes(group = same_group), label = "p.signif")

pp

ggsave(paste0(my_wd,"/boxplot_projs_InterIntraCohort_notch", ".png"), plot = pp, width = 15,   # tamaño más razonable
height = 8,   # tamaño más razonable
units = "in", bg = "white") # usa pulgadas para el tamaño


### grafico de dispersión de
### prop de doble cero vs las diffs de alfa

# distancias elegidas que no correlacionan entre si: 
grouped<- c("shannon_diff", "bp_diff", "faith_diff", "simpson_diff", "singletons_diff" ,  "robbins_diffs", "richness_diff")

totales2<-totales[totales$distancia %in% grouped, ]

my_colors <- c(
  "af_urban_usa"       = "#5C2121",
  "asian"              = "#66CC66",
  "cameroon"           = "#0066FF",
  "children_eu_africa" = "#0000CC",
  "citizen"            = "#FF3333",
  "hadza"              = "#FFD700",
  "healthy_fibers"     = "#228B22",
  "healthy_japanese"   = "#99FF99",
  "infant"             = "#FFB6C1",
  "infant_singapore"   = "#00FFFF",
  "italian"            = "#FF99CC",
  "L_H"                = "#00CCCC",
  "lcarb_lfat"         = "#FF00FF",
  "obese"              = "#FFCC00",
  "ovese"              = "#0033CC",
  "swedish"            = "#CC66FF"
)

totales2$project <- as.factor(totales2$project)
names(my_colors) <- levels(totales2$project)

# coloreado por proyecto
ggplot(totales2, aes(x = dzero, y = valor, color = project)) +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_wrap(~distancia, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Double Zero Frequency",
    y = "Diff alpha matrics",
    title = "Double zero freq vs alpha metrics difference",
    color = "Project"
  ) +
  scale_color_manual(values = my_colors) +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right"
  )


### grafico de dispersión de
### prop de doble cero vs las diffs de alfa

# distancias elegidas que no correlacionan entre si: 
grouped<- c("shannon_diff", "bp_diff", "faith_diff", "simpson_diff", "singletons_diff" ,  "robbins_diffs", "richness_diff")

totales2<-totales[totales$distancia %in% grouped, ]

my_colors <- c(
  "af_urban_usa"       = "#5C2121",
  "asian"              = "#66CC66",
  "cameroon"           = "#0066FF",
  "children_eu_africa" = "#0000CC",
  "citizen"            = "#FF3333",
  "hadza"              = "#FFD700",
  "healthy_fibers"     = "#228B22",
  "healthy_japanese"   = "#99FF99",
  "infant"             = "#FFB6C1",
  "infant_singapore"   = "#00FFFF",
  "italian"            = "#FF99CC",
  "L_H"                = "#00CCCC",
  "lcarb_lfat"         = "#FF00FF",
  "obese"              = "#FFCC00",
  "ovese"              = "#0033CC",
  "swedish"            = "#CC66FF"
)

totales2$project <- as.factor(totales2$project)
names(my_colors) <- levels(totales2$project)

# coloreado por proyecto
ggplot(totales2, aes(x = dzero, y = valor, color = project)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Double Zero Frequency",
    y = "Diff alpha matrics",
    title = "Double zero freq vs alpha metrics difference",
    color = "Project"
  ) +
  scale_color_manual(values = my_colors) +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right"
  )


