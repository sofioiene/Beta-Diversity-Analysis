library(reshape2)
library(gridExtra)
library(tidyverse)
#library(ape)
library(magrittr)
#library(plotly)
library(kableExtra)
#library(Polychrome)
library(ggplot2)
library(ggpmisc)
library(vegan)
library(cluster)
library(DescTools)

contarDobleZero <- function(asv_table, s1, s2) {
  len <- asv_table %>% filter(.data[[s1]] == 0 & .data[[s2]] == 0) %>% count()
  return(len[[1]])
}

my_wd<-"/Users/ignaciocassol/Dropbox/research/beta_analysis/"

rds_files<-list.files(paste0(my_wd, "2. files tsv"), pattern = "tsv")

## Función 1/3

getAllNeighbors <- function(dist_matrix) {
  # Función interna para obtener las muestras más cercanas
  get_nearest_samples <- function(dist_matrix, sample_name) {
    # Excluir la diagonal (la distancia a sí mismo)
    dist_to_sample <- dist_matrix[sample_name, ]
    dist_to_sample <- dist_to_sample[!names(dist_to_sample) %in% sample_name]  # Eliminar la distancia a sí mismo
    
    # Ordenar las distancias de menor a mayor (más cercanas primero)
    sorted_samples <- sort(dist_to_sample)
    
    return(sorted_samples)
  }
  
  # Crear una lista con las muestras más cercanas para cada muestra
  all_nearest_samples <- list()
  for (sample in rownames(dist_matrix)) {
    all_nearest_samples[[sample]] <- get_nearest_samples(dist_matrix, sample)
  }
  
  return(all_nearest_samples)
}

## Función 2/3

rank_samples_from_list <- function(neighbor_list) {
  lapply(neighbor_list, function(x) {
    names(sort(x))  # Ordena por valor (menor a mayor), devuelve los nombres
  })
}

## Función 3/3

compare_kendall_for_all_fast <- function(ranking1, ranking2) {
  samples <- names(ranking1)
  pb <- progress_bar$new(
    format = "Comparando [:bar] :percent en :elapsed",
    total = length(samples), clear = FALSE, width = 60
  )
  
  res <- numeric(length(samples))
  names(res) <- samples
  
  for (i in seq_along(samples)) {
    pb$tick()
    sample <- samples[i]
    
    x <- ranking1[[sample]]
    y <- ranking2[[sample]]
    
    # Verificar que hay valores en común
    common <- intersect(x, y)
    if (length(common) < 2) {
      res[i] <- NA
      next
    }
    
    x_ranks <- match(common, x)
    y_ranks <- match(common, y)
    
    # Calcular Kendall Tau A
    res[i] <- KendallTauA(x_ranks, y_ranks)
  }
  
  return(res)
}


#for (RDS_FILE_NAME in rev(rds_files)){

  RDS_FILE_NAME<-(rds_files)[1]
  
  myrds_file<-RDS_FILE_NAME

  print(paste0("proyecto: ", RDS_FILE_NAME))
  asvtable <- read_tsv(paste0(my_wd,"2. files tsv/",myrds_file), comment = "#", skip=1)
  i=0
  
  latabla<-read_delim(paste0(my_wd, "7.manifest_consolidated/", myrds_file), 
                      delim = NULL, escape_double = FALSE, 
                      trim_ws = TRUE)
  
  asvtable<-asvtable[,colnames(asvtable)%in%latabla$`sample-id`]
  
  asvtable<-asvtable
  mimuestra<-asvtable
  
  dist_gower <- daisy(decostand(t(asvtable), "pa"), metric = "gower")
  dist_chisq = vegdist(decostand(t(asvtable), "total"), "chisq")
  dist_manhattan_log = vegdist(decostand(t(asvtable), "log", logbase = 2), "manhattan")
  dist_raup = vegdist(decostand(t(asvtable), "pa"), "raup")
  dist_mountford = vegdist(decostand(t(asvtable), "pa"), "mountford")
  dist_manhattan_raw = vegdist(t(asvtable), "manhattan")
  dist_cao = vegdist(t(asvtable), "cao")
  dist_chao = vegdist(t(asvtable), "chao")
  dist_bray_edger = vegdist(t(asvtable), method = "bray")
    
  neighbors_gower <- getAllNeighbors(as.matrix(dist_gower))
  neighbors_chisq = getAllNeighbors(as.matrix(dist_chisq))
  neighbors_manhattan_log = getAllNeighbors(as.matrix(dist_manhattan_log))
  neighbors_raup = getAllNeighbors(as.matrix(dist_raup))
  neighbors_mountford = getAllNeighbors(as.matrix(dist_mountford))
  neighbors_manhattan_raw = getAllNeighbors(as.matrix(dist_manhattan_raw))
  neighbors_cao = getAllNeighbors(as.matrix(dist_cao))
  neighbors_chao = getAllNeighbors(as.matrix(dist_chao))
  neighbors_bray_edger = getAllNeighbors(as.matrix(dist_bray_edger))
  
  rankings_gower <- rank_samples_from_list(neighbors_gower)
  rankings_chisq = rank_samples_from_list(neighbors_chisq)
  rankings_manhattan_log = rank_samples_from_list(neighbors_manhattan_log)
  rankings_raup = rank_samples_from_list(neighbors_raup)
  rankings_mountford = rank_samples_from_list(neighbors_mountford)
  rankings_manhattan_raw = rank_samples_from_list(neighbors_manhattan_raw)
  rankings_cao = rank_samples_from_list(neighbors_cao)
  rankings_chao = rank_samples_from_list(neighbors_chao)
  rankings_bray_edger = rank_samples_from_list(neighbors_bray_edger)

  df_gower <- enframe(rankings_gower, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_chisq <- enframe(rankings_chisq, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_manhattan_log <- enframe(rankings_manhattan_log, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_raup <- enframe(rankings_raup, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_mountford <- enframe(rankings_mountford, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_manhattan_raw <- enframe(rankings_manhattan_raw, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_cao <- enframe(rankings_cao, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_chao <- enframe(rankings_chao, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  df_bray_edger <- enframe(rankings_bray_edger, name = "sample", value = "neighbors") %>%
    unnest(neighbors) %>%
    group_by(sample) %>%
    mutate(rank = paste0("neighbor_", row_number())) %>%
    pivot_wider(names_from = rank, values_from = neighbors) %>%
    ungroup()
  
  neighbor_cols_gower <- names(df_gower)[-1]
  neighbor_cols_chisq <- names(df_chisq)[-1]
  neighbor_cols_manhanttan_log <- names(df_manhattan_log)[-1]
  neighbor_cols_raup <- names(df_raup)[-1]
  neighbor_cols_mountford <- names(df_mountford)[-1]
  neighbor_cols_manhattan_raw <- names(df_manhattan_raw)[-1]
  neighbor_cols_cao <- names(df_cao)[-1]
  neighbor_cols_chao <- names(df_chao)[-1]
  neighbor_cols_bray_edger <- names(df_bray_edger)[-1]
  
  df_doble_zero_gower <- df_gower
  df_doble_zero_chisq <- df_chisq
  df_doble_zero_manhattan_log <- df_manhattan_log
  df_doble_zero_raup <- df_raup
  df_doble_zero_mountford <- df_mountford
  df_doble_zero_manhattan_raw <- df_manhattan_raw
  df_doble_zero_cao <- df_cao
  df_doble_zero_chao <- df_chao
  df_doble_zero_bray_edger <- df_bray_edger
  
  # para Gower:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_gower))) {
    sample_i <- df_gower$sample[i]
    
    for (col in neighbor_cols_gower) {
      neighbor_sample <- df_gower[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_gower[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_gower, file = paste0(my_wd, "kendall_test/df_doble_zero_gower.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # para Chisq:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_chisq))) {
    sample_i <- df_chisq$sample[i]
    
    for (col in neighbor_cols_chisq) {
      neighbor_sample <- df_chisq[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_chisq[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_chisq, file = paste0(my_wd, "kendall_test/df_doble_zero_chisq.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # para manhattan_log:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_manhattan_log))) {
    sample_i <- df_manhattan_log$sample[i]
    
    for (col in neighbor_cols_manhanttan_log) {
      neighbor_sample <- df_manhattan_log[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_manhattan_log[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_manhattan_log, file = paste0(my_wd, "kendall_test/df_doble_zero_manhattan_log.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # para raup:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_raup))) {
    sample_i <- df_raup$sample[i]
    
    for (col in neighbor_cols_raup) {
      neighbor_sample <- df_raup[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_raup[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_raup, file = paste0(my_wd, "kendall_test/df_doble_zero_raup.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # para mountford:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_mountford))) {
    sample_i <- df_mountford$sample[i]
    
    for (col in neighbor_cols_mountford) {
      neighbor_sample <- df_mountford[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_mountford[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_mountford, file = paste0(my_wd, "kendall_test/df_doble_zero_mountford.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # para manhattan_raw:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_manhattan_raw))) {
    sample_i <- df_manhattan_raw$sample[i]
    
    for (col in neighbor_cols_manhattan_raw) {
      neighbor_sample <- df_manhattan_raw[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_manhattan_raw[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_manhattan_raw, file = paste0(my_wd, "kendall_test/df_doble_zero_manhattan_raw.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # para cao:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_cao))) {
    sample_i <- df_cao$sample[i]
    
    for (col in neighbor_cols_cao) {
      neighbor_sample <- df_cao[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_cao[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_cao, file = paste0(my_wd, "kendall_test/df_doble_zero_cao.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # para chao:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_chao))) {
    sample_i <- df_chao$sample[i]
    
    for (col in neighbor_cols_chao) {
      neighbor_sample <- df_chao[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_chao[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_chao, file = paste0(my_wd, "kendall_test/df_doble_zero_chao.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # para bray_edger:
  # Iteramos sobre cada fila y columna vecina
  for (i in seq_len(nrow(df_bray_edger))) {
    sample_i <- df_bray_edger$sample[i]
    
    for (col in neighbor_cols_bray_edger) {
      neighbor_sample <- df_bray_edger[[col]][i]
      
      # Calculamos proporción de doble cero y asignamos
      df_doble_zero_bray_edger[[col]][i] <- contarDobleZero(mimuestra, sample_i, neighbor_sample)
    }
  }
  write.table(df_doble_zero_bray_edger, file = paste0(my_wd, "kendall_test/df_doble_zero_bray_edger.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  
  # me voy a quedar solamente con estos rankings para comparar: 

  # Lista de rankings y sus nombres
  ranking_list <- list(
    gower = rankings_gower,
    raup = rankings_raup,
    manhattan_raw = rankings_manhattan_raw,
    cao = rankings_cao,
    chao = rankings_chao,
    bray_edger = rankings_bray_edger
  )
  
  # Comparar todos los pares
  for (i in 1:(length(ranking_list) - 1)) {
    for (j in (i + 1):length(ranking_list)) {
      
      # Obtener nombres y objetos
      r1_name <- names(ranking_list)[i]
      r2_name <- names(ranking_list)[j]
      r1 <- ranking_list[[i]]
      r2 <- ranking_list[[j]]
      
      # Calcular correlación de Kendall
      kendall_scores <- compare_kendall_for_all_fast(r1, r2)
      
      # Generar nombre de archivo
      output_file <- paste0(
        my_wd, "kendall_test/kendall_", r1_name, "_vs_", r2_name, RDS_FILE_NAME
      )
      
      # Guardar
      write.table(kendall_scores, file = output_file, sep = "\t", row.names = TRUE, quote = FALSE)
    }
  }
  
  #kendall_scores <- compare_kendall_for_all_fast(jaccard_rankings, bray_rankings)
  #print(kendall_scores)
  
  #write.table(kendall_scores, file = paste0(my_wd, "kendall_test/kendall_Bray_vs_Jacc",RDS_FILE_NAME), sep = "\t", row.names = TRUE, quote = FALSE)
  

  ## para plotear:
    
  kendall_scores <- read_delim(paste0(my_wd, "kendall_test/kendall_106_infant_singapore_Jacc_vs_Gower.tsv"), 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)
  
  kendall_scores <- setNames(kendall_scores$tau, kendall_scores$sample)
  
  # Gráfico 1/2
  kendall_df <- data.frame(
    sample = names(kendall_scores),
    kendall_tau = as.numeric(kendall_scores)
  )
  
  kendall_df <- kendall_df %>%
    arrange(desc(kendall_tau)) %>%
    mutate(sample = factor(sample, levels = sample))  # ordena las muestras
  
  kendall_df$sample <- factor(kendall_df$sample, levels = rev(levels(kendall_df$sample)))
  
  ggplot(kendall_df, aes(x = sample, y = 1, fill = kendall_tau)) +
    geom_tile(color = "white") +
    #geom_text(aes(label = round(kendall_tau, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                         name = "Tau") +
    labs(title = "Jaccard vs Bray-Curtis",
         subtitle = "Concordance rankings", 
         x = "", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    coord_flip()
  
  
  # Gráfico 2/2
  
  ## graficamos el histogramas de frecuencias de concordancia: 

  df <- data.frame(
    sample = names(kendall_scores),
    kendall_scores = as.numeric(kendall_scores)
  )
  
  # Crear una variable categórica para colorear
  df <- df %>%
    mutate(grupo_concordancia = case_when(
      kendall_scores < 0 ~ "Discordant",
      kendall_scores >= 0 & kendall_scores < 0.3 ~ "Low concordance",
      kendall_scores >= 0.3 & kendall_scores < 0.6 ~ "Moderate concordance",
      kendall_scores >= 0.6 ~ "High concordance"
    ))
  
  # Convertimos en factor para que respete el orden
  df$grupo_concordancia <- factor(df$grupo_concordancia, levels = c("Discordant", "Low concordance", "Moderate concordance", "High concordance"))
  
  # Graficamos
  ggplot(df, aes(x = kendall_scores, fill = grupo_concordancia)) +
    geom_histogram(binwidth = 0.05) +
    scale_fill_manual(values = c("red", "orange", "skyblue", "forestgreen")) +
    labs(title= "Distribution of Kendall Tau Scores",
         subtitle = "by Concordance Level", 
         x = "Kendall Tau", y = "Frecuency", fill = "Concordance level") +
    theme_minimal()
  
#}

##################
  
  mysamples<-names(kendall_scores[kendall_scores>0])
  
  #tableaux<-asvtable[, colnames(asvtable)%in%mysamples]
  
  #dim(tableaux)
  
  #tableaux <- as.matrix(sapply(tableaux, as.numeric))
  
  cambios_fuertes <- names(kendall_scores)[kendall_scores < 0]
  
  alpha_values_good_projects <- read_csv(paste0(my_wd,"alpha_values_good_projects"))
  alpha_values<-alpha_values_good_projects[alpha_values_good_projects$proj=="106_infant_singapore",]

  alpha_values <- alpha_values %>%
    mutate(categ = ifelse(sample %in% cambios_fuertes, "Change", "No Change"))
  
  # Pasamos a formato largo
  alpha_long <- alpha_values %>%
    select(sample, categ, shannon, bp, richness, faith, simpson) %>%
    pivot_longer(
      cols = -c(sample, categ),
      names_to = "metric",
      values_to = "value"
    )
  
  # Graficamos
  ggplot(alpha_long, aes(x = value, fill = categ)) +
    geom_histogram(position = "dodge", alpha = 0.5, bins = 30) +
    scale_fill_manual(
      values = c("Change" = "red", "No Change" = "blue"),
      labels = c("Change" = "Discordant", "No Change" = "Concordant")
    ) +
    facet_wrap(~ metric, scales = "free", ncol = 2) +
    labs(
      title = "Alpha diversity metrics by group",
      subtitle = "Samples grouped according to tau < 0",
      x = "Metric value",
      y = "Frequency",
      fill = "Sample group"
    ) +
    theme_minimal()
  
  
  # Extraer filas (muestras) de la tabla transpuesta
  asv_transpuesta <- t(asvtable)
  asv_muestras_cambio <- asv_transpuesta[rownames(asv_transpuesta) %in% cambios_fuertes, ]
  asv_muestras_cambioB <- asv_transpuesta[!rownames(asv_transpuesta) %in% cambios_fuertes, ]
  
  # Riqueza (número de ASVs distintos de 0)
  
  richnessCambia <- (cbind(as.data.frame(rowSums(asv_muestras_cambio > 0)), "Change"))
  richnessNoCambia <- (cbind(as.data.frame(rowSums(asv_muestras_cambioB > 0)), "No Change"))
  names(richnessCambia)<-c("porcentaje", "categ")
  names(richnessNoCambia)<-c("porcentaje", "categ")
  nueRichness<-rbind(richnessCambia, richnessNoCambia)
  
  
  readsCambia<-(cbind(as.data.frame(rowSums(asv_muestras_cambio)), "Change"))
  readsNoCambia<-(cbind(as.data.frame(rowSums(asv_muestras_cambioB)), "No change"))
  names(readsCambia)<-c("porcentaje", "categ")
  names(readsNoCambia)<-c("porcentaje", "categ")
  nueReads<-rbind(readsCambia, readsNoCambia)
  
  pruCambia<-(cbind(as.data.frame(rowMeans(asv_muestras_cambio == 0)), "Change"))
  pruNoCambia<-(cbind(as.data.frame(rowMeans(asv_muestras_cambioB == 0)), "No change"))
  names(pruCambia)<-c("porcentaje", "categ")
  names(pruNoCambia)<-c("porcentaje", "categ")
  nuew<-rbind(pruCambia, pruNoCambia)
  
  ggplot(nuew, aes(x = porcentaje, fill = categ)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Change" = "red", "No change" = "blue")) +
    labs(
      title = "Percentage of zeros per sample", 
      subtitle = "Samples grouped according to tau < 0", 
      x = "Percentage of zeros per sample", 
      y = "Frequency", 
      fill = "Sample group"
    ) +
    theme_minimal()
  
  ggplot(nuew, aes(x = porcentaje, fill = categ)) +
    geom_histogram(position = "dodge", alpha = 0.5) +
    scale_fill_manual(
      values = c("Change" = "red", "No change" = "blue"),
      labels = c("Change" = "Discordant", "No change" = "Concordant")
    ) +
    labs(
      title = "Percentage of zeros per sample", 
      subtitle = "Samples grouped according to tau < 0", 
      x = "Percentage of zeros per sample", 
      y = "Frequency", 
      fill = "Sample group"
    ) +
    theme_minimal()
  
  
  ggplot(nueReads, aes(x = porcentaje, fill = categ)) +
    geom_density(alpha = 0.5) +
    labs(title= "Jaccard vs Gower - Read size per sample", 
         subtitle = "samples groups according to tau<0", 
         x = "Read size per sample", y = "Frecuency", fill = "Sample group") 
  
  ggplot(nueReads, aes(x = porcentaje, fill = categ)) +
    geom_histogram(position = "dodge", alpha = 0.5) +
    scale_fill_manual(
      values = c("Change" = "red", "No change" = "blue"),
      labels = c("Change" = "Discordant", "No change" = "Concordant")
    ) +
    labs(
      title = "Number or reads", 
      subtitle = "Samples grouped according to tau < 0", 
      x = "Number of reads", 
      y = "Frequency", 
      fill = "Sample group"
    ) +
    theme_minimal()
  
  
  ggplot(nueRichness, aes(x = porcentaje, fill = categ)) +
    geom_density(alpha = 0.5) +
    labs(title= "Jaccard vs Gower - Richness per sample", 
         subtitle = "samples groups according to tau<0", 
         x = "Richness per sample", y = "Frecuency", fill = "Sample group") 

  ggplot(nueRichness, aes(x = porcentaje, fill = categ)) +
    geom_histogram(position = "dodge", alpha = 0.5) +
    scale_fill_manual(
      values = c("Change" = "red", "No Change" = "blue"),
      labels = c("Change" = "Discordant", "No change" = "Concordant")
    ) +
    labs(
      title = "Richness", 
      subtitle = "Samples grouped according to tau < 0", 
      x = "Richness", 
      y = "Frequency", 
      fill = "Sample group"
    ) +
    theme_minimal()  
  