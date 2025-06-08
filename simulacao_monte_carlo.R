# --- FUNÇÃO BASE DE SIMULAÇÃO COM DADOS FALTANTES ----
fn_simula_colada <- function(n_original, tau, m, phi = 0.7, sigma = 1, seed = NULL) {
  if (n_original <= 0 | tau <= 0 | m <= 0) {
    stop('n_original, tau e m devem ser positivos.')
  }
  if (n_original <= tau + m) {
    stop('n_original deve ser maior que tau + m.')
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  gera_ar1 <- function(n, phi, sigma) {
    e <- rnorm(n, mean = 0, sd = sigma)
    x <- numeric(n)
    for (t in 2:n) {
      x[t] <- phi * x[t-1] + e[t]
    }
    return(x)
  }
  
  encontra_maior_segmento <- function(x) {
    is_na <- is.na(x)
    rle_na <- rle(!is_na)
    comprimentos <- rle_na$lengths
    valores <- rle_na$values
    if (all(!valores)) return(numeric(0))
    idx_max <- which.max(comprimentos * valores)
    posicoes <- cumsum(comprimentos)
    start <- ifelse(idx_max == 1, 1, posicoes[idx_max - 1] + 1)
    end <- posicoes[idx_max]
    return(x[start:end])
  }
  
  X <- gera_ar1(n_original, phi, sigma)
  indices_obs <- c(1:tau, (tau + m + 1):n_original)
  X_NA <- X
  X_NA[-indices_obs] <- NA
  Z <- X_NA[!is.na(X_NA)]
  media_obs <- mean(X_NA, na.rm = TRUE)
  X_imputado_media <- X_NA
  X_imputado_media[is.na(X_imputado_media)] <- media_obs
  dados_observados <- X_NA[!is.na(X_NA)]
  X_imputado_hotdeck <- X_NA
  X_imputado_hotdeck[is.na(X_imputado_hotdeck)] <- sample(dados_observados, sum(is.na(X_NA)), replace = TRUE)
  maior_parte <- encontra_maior_segmento(X_NA)
  
  return(list(
    serie_original = X,
    serie_com_NA = X_NA,
    serie_colada = Z,
    serie_imputada_media = X_imputado_media,
    serie_imputada_hotdeck = X_imputado_hotdeck,
    maior_parte = maior_parte,
    indices_observados = indices_obs,
    tau = tau,
    m = m
  ))
}

# --- FUNÇÃO DE MONTE CARLO ----
fn_monte_carlo <- function(n_simulacoes, n_original, tau, m, phi = 0.7, sigma = 1,
                           seed_inicial = NULL,
                           metodo = c("original", "colada", "imputada_media", "imputada_hotdeck", "maior_parte")) {
  metodo <- match.arg(metodo)
  medias <- numeric(n_simulacoes)
  variancias <- numeric(n_simulacoes)
  acf1 <- numeric(n_simulacoes)
  acf2 <- numeric(n_simulacoes)
  acf3 <- numeric(n_simulacoes)
  
  for (i in 1:n_simulacoes) {
    seed_atual <- ifelse(is.null(seed_inicial), NULL, seed_inicial + i)
    resultado <- fn_simula_colada(n_original, tau, m, phi, sigma, seed = seed_atual)
    
    serie_escolhida <- switch(metodo,
                              original = resultado$serie_original,
                              colada = resultado$serie_colada,
                              imputada_media = resultado$serie_imputada_media,
                              imputada_hotdeck = resultado$serie_imputada_hotdeck,
                              maior_parte = resultado$maior_parte)
    
    # Média
    medias[i] <- mean(serie_escolhida, na.rm = TRUE)
    
    # Variância
    variancias[i] <- var(serie_escolhida, na.rm = TRUE)
    
    # Autocovariâncias
    acf1[i] <- mean((serie_escolhida[2:n_original] - mean(serie_escolhida, na.rm = TRUE)) *
                      (serie_escolhida[1:(n_original - 1)] - mean(serie_escolhida, na.rm = TRUE)), na.rm = TRUE)
    acf2[i] <- mean((serie_escolhida[3:n_original] - mean(serie_escolhida, na.rm = TRUE)) *
                      (serie_escolhida[1:(n_original - 2)] - mean(serie_escolhida, na.rm = TRUE)), na.rm = TRUE)
    acf3[i] <- mean((serie_escolhida[4:n_original] - mean(serie_escolhida, na.rm = TRUE)) *
                      (serie_escolhida[1:(n_original - 3)] - mean(serie_escolhida, na.rm = TRUE)), na.rm = TRUE)
  }
  
  return(list(medias = medias, variancias = variancias, acf1 = acf1, acf2 = acf2, acf3 = acf3))
}

# --- FUNÇÃO DE CÁLCULO DE MÉTRICAS ----
calcular_metricas <- function(m, phi, n_original, tau) {
  cat("Executando para m =", m, "tau =", tau, "phi =", phi, "\n")
  
  metodos <- c("original", "colada", "imputada_media", "imputada_hotdeck", "maior_parte")
  resultados <- data.frame(Metodo = character(),
                           Media = numeric(),
                           Variancia = numeric(),
                           ACF1 = numeric(),
                           ACF2 = numeric(),
                           ACF3 = numeric(),
                           Vies_Media = numeric(),
                           Vies_Variancia = numeric(),
                           Vies_ACF1 = numeric(),
                           Vies_ACF2 = numeric(),
                           Vies_ACF3 = numeric(),
                           EQM_Media = numeric(),
                           EQM_Variancia = numeric(),
                           EQM_ACF1 = numeric(),
                           EQM_ACF2 = numeric(),
                           EQM_ACF3 = numeric(),
                           Phi = numeric(),
                           N = integer(),
                           Tau = integer(),
                           m = integer(),
                           stringsAsFactors = FALSE)
  
  # Valor real para comparação
  medias_original <- fn_monte_carlo(n_simulacoes, n_original, tau, m, phi, sigma, seed_inicial, metodo = "original")
  valor_real_media <- mean(medias_original$medias)
  valor_real_variancia <- mean(medias_original$variancias)
  valor_real_acf1 <- mean(medias_original$acf1)
  valor_real_acf2 <- mean(medias_original$acf2)
  valor_real_acf3 <- mean(medias_original$acf3)
  
  for (metodo in metodos) {
    resultado_simulacao <- fn_monte_carlo(n_simulacoes, n_original, tau, m, phi, sigma, seed_inicial, metodo = metodo)
    
    # Média, variância e autocovariâncias
    media_m <- mean(resultado_simulacao$medias)
    variancia_m <- mean(resultado_simulacao$variancias)
    acf1_m <- mean(resultado_simulacao$acf1)
    acf2_m <- mean(resultado_simulacao$acf2)
    acf3_m <- mean(resultado_simulacao$acf3)
    
    # Cálculo do viés e EQM para cada métrica
    vies_media <- media_m - valor_real_media
    vies_variancia <- variancia_m - valor_real_variancia
    vies_acf1 <- acf1_m - valor_real_acf1
    vies_acf2 <- acf2_m - valor_real_acf2
    vies_acf3 <- acf3_m - valor_real_acf3
    
    eqm_media <- mean((resultado_simulacao$medias - valor_real_media)^2)
    eqm_variancia <- mean((resultado_simulacao$variancias - valor_real_variancia)^2)
    eqm_acf1 <- mean((resultado_simulacao$acf1 - valor_real_acf1)^2)
    eqm_acf2 <- mean((resultado_simulacao$acf2 - valor_real_acf2)^2)
    eqm_acf3 <- mean((resultado_simulacao$acf3 - valor_real_acf3)^2)
    
    # Armazenando resultados
    resultados <- rbind(resultados, data.frame(
      Metodo = metodo,
      Media = media_m,
      Variancia = variancia_m,
      ACF1 = acf1_m,
      ACF2 = acf2_m,
      ACF3 = acf3_m,
      Vies_Media = vies_media,
      Vies_Variancia = vies_variancia,
      Vies_ACF1 = vies_acf1,
      Vies_ACF2 = vies_acf2,
      Vies_ACF3 = vies_acf3,
      EQM_Media = eqm_media,
      EQM_Variancia = eqm_variancia,
      EQM_ACF1 = eqm_acf1,
      EQM_ACF2 = eqm_acf2,
      EQM_ACF3 = eqm_acf3,
      Phi = phi,
      N = n_original,
      Tau = tau,
      m = m,
      stringsAsFactors = FALSE
    ))
  }
  return(resultados)
}

# --- CONFIGURAÇÕES GERAIS ----
n_simulacoes <- 500
sigma <- 1
seed_inicial <- 123
# --- NOVOS CENÁRIOS --- 
n_values <- c(25, 50, 75, 100, 150, 200)
perc_faltantes <- c(0.10, 0.30, 0.50)

# --- EXECUTA TODAS AS COMBINAÇÕES COM OS NOVOS PARÂMETROS ---
tabela_resultados_novos <- do.call(rbind, lapply(n_values, function(n_original) {
  tau <- floor(n_original * 0.30)  # Tau é 30% de n
  do.call(rbind, lapply(perc_faltantes, function(p) {
    m <- floor(n_original * p)  # m é 10%, 30%, 50% de n
    if (n_original > tau + m) {
      calcular_metricas(m = m, phi = 0.7, n_original = n_original, tau = tau)
    } else {
      warning(paste("Cenário inválido: n =", n_original, "tau =", tau, "m =", m))
      NULL
    }
  }))
}))

# --- FILTRAR APENAS INFORMAÇÕES DE MÉDIA ---- 
tabela_resultados_media_novos <- tabela_resultados_novos[, c("Metodo", "Media", "Vies_Media", "EQM_Media", "Phi", "N", "Tau", "m")]

# --- VISUALIZAÇÃO FINAL ---
print(tabela_resultados_media_novos)


# Reorganiza os dados para formato wide (cada método em uma coluna)
tabela_wide <- tabela_resultados_media_novos %>%
  select(Metodo, Media, Vies_Media, EQM_Media, N, Tau, m) %>%
  pivot_wider(
    names_from = Metodo,
    values_from = c(Media, Vies_Media, EQM_Media)
  )

# Pegando apenas as infos da média 

tabela_media <- tabela_wide %>%
  select(N, Tau, m, starts_with("Media_"))


tabela_vies <- tabela_wide %>%
  select(N, Tau, m, starts_with("Vies_Media"))

tabela_eqm <- tabela_wide %>%
  select(N, Tau, m, starts_with("EQM_Media"))


tabela_eqm
