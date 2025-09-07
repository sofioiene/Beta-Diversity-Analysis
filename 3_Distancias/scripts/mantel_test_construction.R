library(vegan)
library(tidyverse)
library(magrittr)
library(kableExtra)
library(gridExtra)
library(diverse)
library(reshape2)
library(ggpmisc)
library(edgeR)

my_wd<-"/Users/ignaciocassol/Dropbox/research/beta_analysis/"

rds_files<-list.files(paste0(my_wd, "2. files tsv"), pattern = "tsv")

metrics_all_dist_mantels<-{}
all_projs_mantels_R<-{}
all_projs_mantels_C<-{}

for (name in rds_files[1:16]) {
  myrds_file<-name
  print(paste0("proyecto: ", name))
  asvtable <- read_tsv(paste0(my_wd,"2. files tsv/",myrds_file), comment = "#", skip=1) %>%  column_to_rownames(var = "OTU ID")
  i=0
  
  latabla<-read_delim(paste0(my_wd, "7.manifest_consolidated/", myrds_file), 
                      delim = NULL, escape_double = FALSE, 
                      trim_ws = TRUE)
  
  asvtable<-asvtable[,colnames(asvtable)%in%latabla$`sample-id`]
  asv_t <- t(asvtable)
  prevalencia <- rowSums(asvtable > 0) / ncol(asvtable)
  raros <- prevalencia < .1
  table(raros)
  comunes <- !raros
  
  norm_raros <- asvtable[raros, ]
  norm_comunes <- asvtable[comunes, ]
  # quitamos columnas (muestras) vacías:
  
  quitar_samples_comunes<-colnames(norm_comunes)[colSums(norm_comunes)==0]
  quitar_samples_raros<-colnames(norm_raros)[colSums(norm_raros)==0]
  samples_validos <- colnames(asvtable)[!colnames(asvtable) %in% c(quitar_samples_comunes, quitar_samples_raros)]

  norm_raros <- norm_raros[, samples_validos]
  norm_comunes <- norm_comunes[, samples_validos]
  asv_t<-asv_t[samples_validos,]
  
  all_dist_mantels_R<-{}
  all_dist_mantels_C<-{}

    # Crear objeto DGEList y normalizar con TMM (edgeR)
    dge <- DGEList(counts = asvtable[,samples_validos])
    dge <- calcNormFactors(dge, method = "TMM")
    norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
    
    # Transponer porque vegdist espera muestras en filas
    norm_counts_t <- t(norm_counts)
    
    dist_list <- list(
      chisq = vegdist(decostand(asv_t, "total"), "chisq"),
      binomial = vegdist(decostand(asv_t, "pa"), "binomial"),
      bray_edger = vegdist(norm_counts_t, method = "bray"),
      raup = vegdist(decostand(asv_t, "pa"), "raup"),
      mountford = vegdist(decostand(asv_t, "pa"), "mountford"),
      manhattan_raw = vegdist(asv_t, "manhattan"),
      cao = vegdist(asv_t, "cao"),
      chao = vegdist(asv_t, "chao")
    )
    
    # Transponer tablas para vegdist
    norm_raros_t <- t(norm_raros)
    norm_comunes_t <- t(norm_comunes)
    
    # Crear listas de distancias específicas para raros y comunes
    dist_list_raros <- list(
      chisq = vegdist(decostand(norm_raros_t, "total"), "chisq"),
      binomial = vegdist(decostand(norm_raros_t, "pa"), "binomial"),
      bray_edger = vegdist(cpm(DGEList(counts = norm_raros), normalized.lib.sizes = TRUE) %>% t(), method = "bray"),
      raup = vegdist(decostand(norm_raros_t, "pa"), "raup"),
      mountford = vegdist(decostand(norm_raros_t, "pa"), "mountford"),
      manhattan_raw = vegdist(norm_raros_t, "manhattan"),
      cao = vegdist(norm_raros_t, "cao"),
      chao = vegdist(norm_raros_t, "chao")
    )
    
    dist_list_comunes <- list(
      chisq = vegdist(decostand(norm_comunes_t, "total"), "chisq"),
      binomial = vegdist(decostand(norm_comunes_t, "pa"), "binomial"),
      bray_edger = vegdist(cpm(DGEList(counts = norm_comunes), normalized.lib.sizes = TRUE) %>% t(), method = "bray"),
      raup = vegdist(decostand(norm_comunes_t, "pa"), "raup"),
      mountford = vegdist(decostand(norm_comunes_t, "pa"), "mountford"),
      manhattan_raw = vegdist(norm_comunes_t, "manhattan"),
      cao = vegdist(norm_comunes_t, "cao"),
      chao = vegdist(norm_comunes_t, "chao")
    )
  #i <- 0
  for (distancz in names(dist_list)) {
      
    dist_total <- as.matrix(dist_list[[distancz]])
    dist_raros <- as.matrix(dist_list_raros[[distancz]])
    dist_comunes <- as.matrix(dist_list_comunes[[distancz]])

    # Igualamos el orden de muestras en las 3 matrices
    comunes_ids <- Reduce(intersect, list(rownames(dist_total), rownames(dist_raros), rownames(dist_comunes)))
    
    pru1 <- mantel(dist_total[comunes_ids, comunes_ids], dist_raros[comunes_ids, comunes_ids])
    pru2 <- mantel(dist_total[comunes_ids, comunes_ids], dist_comunes[comunes_ids, comunes_ids])
    
    all_dist_mantels_R <- rbind(all_dist_mantels_R, round(pru1$statistic, 2))
    all_dist_mantels_C <- rbind(all_dist_mantels_C, round(pru2$statistic, 2))
  }
    
    
  names(all_dist_mantels_R)<-str_remove(myrds_file, ".tsv")
  names(all_dist_mantels_C)<-str_remove(myrds_file, ".tsv")
  all_projs_mantels_R<-cbind(all_projs_mantels_R, all_dist_mantels_R)
  all_projs_mantels_C<-cbind(all_projs_mantels_C, all_dist_mantels_C)
  
  metrics_all_dist_mantels<-rbind(metrics_all_dist_mantels, c(table(raros)[2]/nrow(asvtable),length(c(quitar_samples_comunes, quitar_samples_raros)),dim(asvtable)[1],dim(asvtable)[2]))
  
}

colnames(all_projs_mantels_C)<-str_remove(rds_files[1:16],".tsv")
colnames(all_projs_mantels_R)<-str_remove(rds_files[1:16],".tsv")
rownames(all_projs_mantels_C)<-names(dist_list)
rownames(all_projs_mantels_R)<-names(dist_list)

View(all_projs_mantels_R)
colnames(metrics_all_dist_mantels)<-c("prop_rare", "sample_removed", "ASVs", "samples")
View(metrics_all_dist_mantels)

write.table(all_projs_mantels_R, file = paste0(my_wd, "mantel_test/", "mantel_of_rare_samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all_projs_mantels_C, file = paste0(my_wd, "mantel_test/", "mantel_of_common_samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(metrics_all_dist_mantels, file = paste0(my_wd, "mantel_test/", "experiment_metrics_mantel_test.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# tabla 1

df <- as.data.frame(all_projs_mantels_R)

# Convertimos los valores mayores a 0.75 en celdas rojas
df_coloreado <- df %>%
  mutate(across(
    everything(),
    ~ cell_spec(round(.x, 2), background = ifelse(.x > 0.75, "red", "white"),
                color = ifelse(.x > 0.75, "white", "black"))
  ))

# Añadir nombres de fila si los tenés
rownames(df_coloreado) <- rownames(all_projs_mantels_R)

# Imprimir con kable
kable(df_coloreado, escape = FALSE, format = "html") %>%
  kable_styling("striped", full_width = F)


# tabla 2

df <- as.data.frame(all_projs_mantels_C)

# Convertimos los valores mayores a 0.75 en celdas rojas
df_coloreado <- df %>%
  mutate(across(
    everything(),
    ~ cell_spec(round(.x, 2), background = ifelse(.x > 0.75, "red", "white"),
                color = ifelse(.x > 0.75, "white", "black"))
  ))

# Añadir nombres de fila si los tenés
rownames(df_coloreado) <- rownames(all_projs_mantels_R)

# Imprimir con kable
kable(df_coloreado, escape = FALSE, format = "html") %>%
  kable_styling("striped", full_width = F)
