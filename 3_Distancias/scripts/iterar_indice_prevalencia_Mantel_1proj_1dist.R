library(tidyverse)
library(edgeR)
library(vegan)

## prueba de sensibilidad del thershold de prevalencia
## en 1 proyecto y una distancia

# PARAMETROS ----
my_wd<-"/home/ignacio/Dropbox/research/beta_analysis/"

rds_files<-list.files(paste0(my_wd, "2. files tsv"), pattern = "tsv")

dist_list <- list(
  #chisq = vegdist(decostand(asv_t, "total"), "chisq"),
  #binomial = vegdist(decostand(asv_t, "pa"), "binomial"),
  #raup = vegdist(decostand(asv_t, "pa"), "raup"),
  #mountford = vegdist(decostand(asv_t, "pa"), "mountford"),
  #manhattan_raw = vegdist(asv_t, "manhattan"),
  #chao = vegdist(asv_t, "chao"),
  cao = vegdist(asv_t, "cao"),
  #bray_edger = vegdist(norm_counts_t, method = "bray")
)

dist_name <- "binomial"
thresholds <- seq(0.05, 0.95, by = 0.1)

# Lista para guardar resultados
mantel_results <- data.frame(threshold = numeric(),
                             mantel_r_comunes = numeric(),
                             ntaxacomunes = numeric(),
                             nceros = numeric(), 
                             mantel_r_raros = numeric())

for (name in rds_files) {
  start_time <- Sys.time()
  myrds_file<-name
  cat(paste0("empiezo con ", myrds_file, start_time, "\n"))

  suppressMessages({
    asvtable <- read_tsv(paste0(my_wd,"2. files tsv/",myrds_file), comment = "#", skip=1) %>%  column_to_rownames(var = "OTU ID")
  })
  
  suppressMessages({  
  latabla<-read_delim(paste0(my_wd, "7.manifest_consolidated/", myrds_file), 
                      delim = NULL, escape_double = FALSE, 
                      trim_ws = TRUE)
  })
    
  asvtable<-asvtable[,colnames(asvtable)%in%latabla$`sample-id`]

  ntaxa<-dim(asvtable)[1]
  porcentaje_ceros <- sum(asvtable == 0, na.rm = TRUE) / (nrow(asvtable) * ncol(asvtable)) * 100
  asv_t <- t(asvtable)

  if(dist_name=="bray_edger"){
    # MATRIZ NORMALIZADA (TMM) ----
    dge <- DGEList(counts = asvtable)
    dge <- calcNormFactors(dge, method = "TMM")
    norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
    norm_counts_t <- t(norm_counts)
    dist_total <- vegdist(norm_counts_t, method = "bray")
  } else {
    dist_total <- vegdist(decostand(asv_t, "pa"), "binomial")
  }
  
# Iteramos sobre thresholds de prevalencia
for (thresh in thresholds) {
  cat(paste0("iteracion ", thresh, " del proyecto ", myrds_file, "\n"))
  
  prevalencia <- rowSums(asvtable > 0) / ncol(asvtable)
  raros <- prevalencia < thresh
  comunes <- !raros
  
  norm_raros <- asvtable[raros, , drop = FALSE]
  norm_comunes <- asvtable[comunes, , drop = FALSE]
  
  # Quitamos muestras con suma 0
  quitar_samples_comunes <- colnames(norm_comunes)[colSums(norm_comunes) == 0]
  quitar_samples_raros <- colnames(norm_raros)[colSums(norm_raros) == 0]
  
  samples_validos <- setdiff(colnames(asvtable), union(quitar_samples_comunes, quitar_samples_raros))
  
  if (length(samples_validos) < 3) next  # Skip si no hay suficientes muestras
  
  norm_comunes <- norm_comunes[, samples_validos]

  taxacomunes<-round(dim(norm_comunes)[1]/ntaxa,3)
  porcentaje_ceros_comunes <- round(sum(norm_comunes == 0, na.rm = TRUE) / (nrow(norm_comunes) * ncol(norm_comunes)) * 100, 2)
  
  norm_raros <- norm_raros[, samples_validos]
  asvtable_filt <- asvtable[, samples_validos]

  if(dist_name=="bray_edger"){
    # Normalizamos con TMM
    dge <- DGEList(counts = asvtable_filt)
    dge <- calcNormFactors(dge)
    norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
  
    # Submatrices raros y comunes
    dge_raros <- DGEList(counts = norm_raros)
    dge_raros <- calcNormFactors(dge_raros)
    norm_raros_cpm <- cpm(dge_raros, normalized.lib.sizes = TRUE)
    
    dge_comunes <- DGEList(counts = norm_comunes)
    dge_comunes <- calcNormFactors(dge_comunes)
    norm_comunes_cpm <- cpm(dge_comunes, normalized.lib.sizes = TRUE)
    
    
    dist_total <- vegdist(t(norm_counts), "manhattan")
    dist_raros <- vegdist(t(norm_raros_cpm), "manhattan")
    dist_comunes <- vegdist(t(norm_comunes_cpm), "manhattan")
  }else {

    dist_total <- vegdist(decostand(asv_t, "pa"), "binomial")
    dist_raros <- vegdist(decostand(t(norm_raros),"pa"), "binomial")
    dist_comunes <- vegdist(decostand(t(norm_comunes),"pa"), "binomial")
  }
    cat(paste0("calcule distancias", myrds_file, "\n" ))

  # Validamos que todas las distancias tengan las mismas muestras
  ids <- Reduce(intersect, list(rownames(as.matrix(dist_total)), 
                                rownames(as.matrix(dist_raros)), 
                                rownames(as.matrix(dist_comunes))))
  if (length(ids) < 3) next
  
  if (anyNA(as.matrix(dist_comunes))) next
  if (anyNA(as.matrix(dist_raros))) next
  
  m1 <- mantel(as.matrix(dist_total)[ids, ids], as.matrix(dist_raros)[ids, ids])
  m2 <- mantel(as.matrix(dist_total)[ids, ids], as.matrix(dist_comunes)[ids, ids])
  cat(paste0("calcule los dos mantel ", myrds_file, "\n"))
  
  mantel_results <- rbind(mantel_results,
                          data.frame(threshold = thresh,
                                     mantel_r_comunes = m2$statistic,
                                     mantel_r_raros = m1$statistic, 
                                     ntaxacomunes = taxacomunes,
                                     nceros = porcentaje_ceros_comunes, 
                                     project = str_remove(myrds_file, ".tsv")))
}
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(paste0("Terminé con ", myrds_file, " en ", round(elapsed_time, 2), " minutos\n"))
}

mantel_results<-mantel_results[!mantel_results$project=="98_combo_diet", ]

write.table(mantel_results, 
            file = paste0(my_wd, "mantel_test/", dist_name, ".tsv"), 
            sep = "\t", row.names = FALSE, quote = TRUE)


## grafico todos los proyectos en dos ejes Y
# mantel_comunes y nceros

# 1. Definir un factor de escala para alinear visualmente los dos valores
# Puede ajustarse según los rangos que tengas
scale_factor <- max(mantel_results$mantel_r_comunes, na.rm = TRUE) / max(mantel_results$ntaxacomunes, na.rm = TRUE)

# 2. Crear el gráfico
ggplot(mantel_results, aes(x = threshold)) +
  geom_line(aes(y = mantel_r_comunes, color = "Mantel R (comunes)"), size = 1) +
  geom_line(aes(y = ntaxacomunes * scale_factor, color = "Porcentaje de ceros"), size = 1) +
  facet_wrap(~ project) +
  scale_y_continuous(
    name = "Mantel R (comunes)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Porcentaje de ceros")
  ) +
  scale_color_manual(
    name = "Variable",
    values = c("Mantel R (comunes)" = "red", "Porcentaje de ceros" = "blue")
  ) +
  labs(
    title = paste0("Mantel R y % de ceros vs Threshold ", dist_name),
    x = "Threshold de prevalencia"
  ) +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "red"),
    axis.title.y.right = element_text(color = "blue"),
    legend.position = "bottom"
  )


summary_mantel <- mantel_results %>%
  group_by(threshold) %>%
  summarise(
    mean_mantel = mean(mantel_r_comunes, na.rm = TRUE),
    sd_mantel = sd(mantel_r_comunes, na.rm = TRUE),
    mean_nceros = mean(nceros, na.rm = TRUE),
    sd_nceros = sd(nceros, na.rm = TRUE),
    mean_taxacomunes = mean(ntaxacomunes, na.rm = TRUE),
    sd_taxacomunes = sd(ntaxacomunes, na.rm = TRUE)
    
  )

scale_factor <- max(summary_mantel$mean_mantel) / max(summary_mantel$mean_nceros)

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
  
  # Línea horizontal en 0.75
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red", size = 0.5) +
  
  # Ejes y etiquetas
  scale_y_continuous(
    name = "Media(Estadístico r del test de Mantel)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Porcentaje de ceros (%)")
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
    title = paste0("Comparación de distancia ", dist_name, " y dispersión de ceros"),
    subtitle = "Promedio del test de Mantel (comunes) vs. porcentaje de ceros",
    x = "Threshold de prevalencia"
  ) +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "red"),
    axis.title.y.right = element_text(color = "blue"),
    legend.position = "top")

p

ggsave(paste0(my_wd, "mantel_test/grafico_", dist_name, "_mantel_vs_nceros.jpg"), plot = p, width = 10, height = 6, dpi = 300)


##### taxacomunes vs mantel

scale_factor <- max(summary_mantel$mean_mantel) / max(summary_mantel$mean_taxacomunes)

p <- ggplot(summary_mantel, aes(x = threshold)) +
  # Línea y ribbon del test de Mantel
  geom_line(aes(y = mean_mantel, color = "Test de Mantel"), size = 1) +
  geom_ribbon(aes(ymin = mean_mantel - sd_mantel,
                  ymax = mean_mantel + sd_mantel,
                  fill = "Test de Mantel"),
              alpha = 0.2) +
  
  # Línea y ribbon de % de ceros (escalado)
  geom_line(aes(y = mean_taxacomunes * scale_factor, color = "% de ceros"), linewidth = 1) +
  geom_ribbon(aes(ymin = (mean_taxacomunes - sd_taxacomunes) * scale_factor,
                  ymax = (mean_taxacomunes + sd_taxacomunes) * scale_factor,
                  fill = "% de ceros"),
              alpha = 0.2) +
  
  # Línea horizontal en 0.75
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red", size = 0.5) +
  
  # Ejes y etiquetas
  scale_y_continuous(
    name = "Media(Estadístico r del test de Mantel)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Taxa comunes (%)")
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
    title = paste0("Comparación de distancia ", dist_name, " y dispersión de ceros"),
    subtitle = "Promedio del test de Mantel (comunes) vs. porcentaje de taxa comunes",
    x = "Threshold de prevalencia"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

p

ggsave(paste0(my_wd, "mantel_test/grafico_", dist_name, "_mantel_vs_ntaxa.jpg"), plot = p, width = 10, height = 6, dpi = 300)


