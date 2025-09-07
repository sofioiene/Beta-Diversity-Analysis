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

contarDobleZero <- function(asv_table, s1, s2) {
  len <- asv_table %>% filter(.data[[s1]] == 0 & .data[[s2]] == 0) %>% count()
  return(len[[1]]/dim(asv_table)[1])
}

j=0

for (name in rds_files) {
  myrds_file<-name
  print(paste0("proyecto: ", name))
  asvtable <- read_tsv(paste0(my_wd,"2. files tsv/",myrds_file), comment = "#", skip=1) %>%  column_to_rownames(var = "OTU ID")
  i=0
  
  latabla<-read_delim(paste0(my_wd, "7.manifest_consolidated/", myrds_file), 
                      delim = NULL, escape_double = FALSE, 
                      trim_ws = TRUE)
  
  asvtable<-asvtable[,colnames(asvtable)%in%latabla$`sample-id`]
  asv_t <- t(asvtable)

  # Crear objeto DGEList y normalizar con TMM (edgeR)
  dge <- DGEList(counts = asvtable)
  dge <- calcNormFactors(dge, method = "TMM")
  norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
  
  # Transponer porque vegdist espera muestras en filas
  norm_counts_t <- t(norm_counts)
  
  dist_list <- list(
    chord = vegdist(decostand(asv_t, "normalize"), "chord"),
    chisq = vegdist(decostand(asv_t, "total"), "chisq"),
    bray_log = vegdist(decostand(asv_t, "log", logbase = 2), "bray"),
    bray_total = vegdist(decostand(asv_t, "total"), "bray"),
    bray_hell = vegdist(decostand(asv_t, "hellinger"), "bray"),
    canberra_log = vegdist(decostand(asv_t, "log", logbase = 2), "canberra"),
    canberra_total = vegdist(decostand(asv_t, "total"), "canberra"),
    clark_log = vegdist(decostand(asv_t, "log", logbase = 2), "clark"),
    clark_total = vegdist(decostand(asv_t, "total"), "clark"),
    eucl_log = vegdist(decostand(asv_t, "log", logbase = 2), "euclidean"),
    eucl_hell = vegdist(decostand(asv_t, "hellinger"), "euclidean"),
    binomial = vegdist(decostand(asv_t, "pa"), "binomial"),
    raup = vegdist(decostand(asv_t, "pa"), "raup"),
    mountford = vegdist(decostand(asv_t, "pa"), "mountford"),
    jaccard = vegdist(decostand(asv_t, "pa"), "jaccard"),
    manhattan_log = vegdist(decostand(asv_t, "log", logbase = 2), "manhattan"),
    manhattan_total = vegdist(decostand(asv_t, "total"), "manhattan"),
    manhattan_raw = vegdist(asv_t, "manhattan"),
    kulc_log = vegdist(decostand(asv_t, "log", logbase = 2), "kulczynski"),
    kulc_total = vegdist(decostand(asv_t, "total"), "kulczynski"),
    kulc_raw = vegdist(asv_t, "kulczynski"),
    cao = vegdist(asv_t, "cao"),
    chao = vegdist(asv_t, "chao"),
    horn = vegdist(decostand(asv_t, "total"), "horn"),
    bray_edger = vegdist(norm_counts_t, method = "bray")
  )
  
  i <- 0
  for (distancz in names(dist_list)) {
    # Ejecutar la función de distancia sobre asv_t
    dist_matrix <- as.matrix(dist_list[[distancz]])
    
    # Convertir a upper triangular
    upper_matrix <- dist_matrix
    upper_matrix[!upper.tri(upper_matrix)] <- NA
    mydatatoplot1 <- melt(upper_matrix, na.rm = TRUE)
    colnames(mydatatoplot1)[3] <- distancz
    
    if (i == 0) {
      mydatatoplot <- mydatatoplot1
      i <- 1
    } else {
      mydatatoplot <- merge(mydatatoplot, mydatatoplot1, by = c("Var1", "Var2"))
    }
  }
  
  # Proporción de doble cero
  mydatatoplot$dzero <- mapply(
    function(v1, v2) contarDobleZero(asvtable, as.character(v1), as.character(v2)),
    mydatatoplot$Var1, mydatatoplot$Var2
  )
  
  # Agregar grupo para Var1
  mydatatoplot <- mydatatoplot %>%
    left_join(latabla %>% dplyr::select(`sample-id`, group), by = c("Var1" = "sample-id")) %>%
    rename(grupoVar1 = group)

    
  # Agregar grupo para Var2
  mydatatoplot <- mydatatoplot %>%
    left_join(latabla %>% dplyr::select(`sample-id`, group), by = c("Var2" = "sample-id")) %>%
    rename(grupoVar2 = group)
  
  # Guardar
  write.table(mydatatoplot, 
              file = paste0(my_wd, "double_zero_all_dist_with_all_norms/", myrds_file), 
              sep = "\t", row.names = FALSE, quote = FALSE)
}
