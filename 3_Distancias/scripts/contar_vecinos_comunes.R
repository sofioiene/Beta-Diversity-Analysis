# datos requeridos
# df_qgower, df_jacc y df_doble_zero_jacc

contar_vecinos_comunes <- function(df1, df2, k = 10) {
  # Verificamos que las muestras estén en el mismo orden
  if (!all(df1$sample == df2$sample)) {
    stop("Los data frames no tienen las muestras en el mismo orden.")
  }
  
  # Para cada fila (muestra), contar cuántos vecinos se repiten
  vecinos_comunes <- sapply(1:nrow(df1), function(i) {
    vecinos_1 <- unlist(df1[i, 2:(k + 1)])
    vecinos_2 <- unlist(df2[i, 2:(k + 1)])
    
    length(intersect(vecinos_1, vecinos_2))
  })
  
  # Devolver el resultado como data frame
  return(data.frame(
    sample = df1$sample,
    vecinos_comunes = vecinos_comunes
  ))
}

resultado_comunes <- contar_vecinos_comunes(df_qgower, df_qjacc, k = 50)
hist(resultado_comunes$vecinos_comunes)

resultado_comunes <- contar_vecinos_comunes(df_qbray, df_qjacc, k = 50)
hist(resultado_comunes$vecinos_comunes)

mean(resultado_comunes$vecinos_comunes)



analizar_exclusiones_por_doble_cero_posicional <- function(df_vecinos_ref, df_vecinos_comp, df_doble_cero, k = 10) {
  if (!all(df_vecinos_ref$sample == df_vecinos_comp$sample)) {
    stop("Los data frames no tienen las muestras en el mismo orden.")
  }
  
  resultados <- lapply(1:nrow(df_vecinos_ref), function(i) {
    sample_i <- df_vecinos_ref$sample[i]
    
    vecinos_ref <- unlist(df_vecinos_ref[i, 2:(k + 1)])
    vecinos_comp <- unlist(df_vecinos_comp[i, 2:(k + 1)])
    
    # Vecinos de referencia que NO están en el otro método
    vecinos_excluidos <- setdiff(vecinos_ref, vecinos_comp)
    
    # Encontrar en qué columnas están esos vecinos excluidos (posición en vecinos_ref)
    posiciones_excluidas <- which(vecinos_ref %in% vecinos_excluidos)
    
    # Usar esas posiciones para obtener los valores de doble cero desde df_doble_cero
    valores_doble_cero <- unlist(df_doble_cero[i, posiciones_excluidas + 1])  # +1 por la columna 'sample'
    
    data.frame(
      sample = sample_i,
      vecinos_excluidos = length(posiciones_excluidas),
      suma_doble_cero = sum(valores_doble_cero),
      promedio_doble_cero = if(length(valores_doble_cero) > 0) mean(valores_doble_cero) else NA
    )
  })
  
  do.call(rbind, resultados)
}

resultado_exclusiones_bray <- analizar_exclusiones_por_doble_cero_posicional(
  df_vecinos_ref = df_qbray,
  df_vecinos_comp = df_qjacc,
  df_doble_cero = df_doble_zero_bray,
  k = 10
)
