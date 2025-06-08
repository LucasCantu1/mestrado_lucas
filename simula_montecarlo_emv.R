library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(tibble)
library(stringr)

# Função aux ---- 

fn_gera_serie_colada <- function(n, phi, sigma, tau, m, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Verificações básicas
  if (tau + m - 1 > n) {
    stop("O intervalo de NA excede o tamanho da série (n).")
  }
  
  # Inicializa a série
  y <- numeric(n)
  y[1] <- rnorm(1, mean = 0, sd = sigma / sqrt(1 - phi^2))  # valor inicial em equilíbrio
  
  # Gera a série AR(1)
  for (t in 2:n) {
    y[t] <- phi * y[t - 1] + rnorm(1, mean = 0, sd = sigma)
  }
  
  # Cria a versão incompleta
  y_incompleta <- y
  y_incompleta[(tau+1):(tau + m)] <- NA
  
  # Imputação pela média
  media_obs <- mean(y_incompleta, na.rm = TRUE)
  y_media <- y_incompleta
  y_media[is.na(y_media)] <- media_obs
  
  # Imputação hot-deck (valor aleatório observado)
  y_hotdeck <- y_incompleta
  valores_observados <- y_hotdeck[!is.na(y_hotdeck)]
  y_hotdeck[is.na(y_hotdeck)] <- sample(valores_observados, sum(is.na(y_hotdeck)), replace = TRUE)
  
  # Série colada = remove os NAs e "gruda" os pedaços
  y_colada_curta <- y_incompleta[!is.na(y_incompleta)]
  y_colada <- rep(NA, n)
  y_colada[1:length(y_colada_curta)] <- y_colada_curta
  
  # Determinar os dois trechos possíveis
  idx_antes <- 1:(tau)
  idx_depois <- (tau + m + 1):n
  
  # Inicializa coluna maior_trecho como NA
  maior_trecho <- rep(NA, n)
  
  # Seleciona e preenche o maior trecho observado
  if (length(idx_antes) >= length(idx_depois)) {
    maior_trecho[idx_antes] <- y[idx_antes]
  } else {
    maior_trecho[idx_depois] <- y[idx_depois]
  }
  
  # Dataframe final com maior_trecho
  df <- data.frame(
    tempo = 1:n,
    completa = y,
    incompleta = y_incompleta,
    imput_media = y_media,
    imput_hotdeck = y_hotdeck,
    serie_colada = y_colada,
    maior_trecho = maior_trecho
  )
  
  return(df)
}



# Estimação EMV (definida antes)
estima_emv_ar1 <- function(y) {
  y <- y[!is.na(y)]
  y_lag <- y[-length(y)]
  y_curr <- y[-1]
  phi_hat <- sum(y_lag * y_curr) / sum(y_lag^2)
  resid <- y_curr - phi_hat * y_lag
  sigma_hat <- sqrt(mean(resid^2))
  return(list(phi = phi_hat, sigma = sigma_hat))
}

# Monte Carlo com EMV
fn_monte_carlo <- function(n_sim, n = 100, phi, sigma, tau = 40, m = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  resultado_total <- tibble()
  
  for (i in 1:n_sim) {
    df <- fn_gera_serie_colada(n = n, phi = phi, sigma = sigma, tau = tau, m = m)
    
    completa      <- df$completa
    imput_media   <- df$imput_media
    imput_hotdeck <- df$imput_hotdeck
    colada        <- df %>% filter(!is.na(serie_colada)) %>% pull(serie_colada)
    maior_trecho  <- df %>% filter(!is.na(maior_trecho)) %>% pull(maior_trecho)
    
    est_completa     <- estima_emv_ar1(completa)
    est_media        <- estima_emv_ar1(imput_media)
    est_hotdeck      <- estima_emv_ar1(imput_hotdeck)
    est_colada       <- estima_emv_ar1(colada)
    est_maior_trecho <- estima_emv_ar1(maior_trecho)
    
    resultado_total <- bind_rows(resultado_total, tibble(
      sim = i,
      phi_real = phi,
      sigma_real = sigma,
      
      phi_completa = est_completa$phi,
      sigma_completa = est_completa$sigma,
      
      phi_media = est_media$phi,
      sigma_media = est_media$sigma,
      
      phi_hotdeck = est_hotdeck$phi,
      sigma_hotdeck = est_hotdeck$sigma,
      
      phi_colada = est_colada$phi,
      sigma_colada = est_colada$sigma,
      
      phi_maior_trecho = est_maior_trecho$phi,
      sigma_maior_trecho = est_maior_trecho$sigma
    ))
  }
  
  return(resultado_total)
}






# Parâmetros dos cenários
phis <- c(0.2, 0.5, 0.8)
sigmas <- c(0.5, 1.0)
cenarios <- expand.grid(phi = phis, sigma = sigmas)

# Roda simulações
resultados_todos <- map2_dfr(
  cenarios$phi, cenarios$sigma,
  ~fn_monte_carlo(n_sim = 500, phi = .x, sigma = .y)
)






# Atualização com EQM
tabela_final <- resultados_todos %>%
  group_by(phi_real, sigma_real) %>%
  summarise(
    # completa
    phi_completa_medio = mean(phi_completa),
    phi_completa_sd = sd(phi_completa),
    phi_completa_eqm = mean((phi_completa - phi_real)^2),
    
    sigma_completa_medio = mean(sigma_completa),
    sigma_completa_sd = sd(sigma_completa),
    sigma_completa_eqm = mean((sigma_completa - sigma_real)^2),
    
    # média
    phi_media_medio = mean(phi_media),
    phi_media_sd = sd(phi_media),
    phi_media_eqm = mean((phi_media - phi_real)^2),
    
    sigma_media_medio = mean(sigma_media),
    sigma_media_sd = sd(sigma_media),
    sigma_media_eqm = mean((sigma_media - sigma_real)^2),
    
    # hotdeck
    phi_hotdeck_medio = mean(phi_hotdeck),
    phi_hotdeck_sd = sd(phi_hotdeck),
    phi_hotdeck_eqm = mean((phi_hotdeck - phi_real)^2),
    
    sigma_hotdeck_medio = mean(sigma_hotdeck),
    sigma_hotdeck_sd = sd(sigma_hotdeck),
    sigma_hotdeck_eqm = mean((sigma_hotdeck - sigma_real)^2),
    
    # colada
    phi_colada_medio = mean(phi_colada),
    phi_colada_sd = sd(phi_colada),
    phi_colada_eqm = mean((phi_colada - phi_real)^2),
    
    sigma_colada_medio = mean(sigma_colada),
    sigma_colada_sd = sd(sigma_colada),
    sigma_colada_eqm = mean((sigma_colada - sigma_real)^2),
    
    # maior trecho
    phi_maior_medio = mean(phi_maior_trecho),
    phi_maior_sd = sd(phi_maior_trecho),
    phi_maior_eqm = mean((phi_maior_trecho - phi_real)^2),
    
    sigma_maior_medio = mean(sigma_maior_trecho),
    sigma_maior_sd = sd(sigma_maior_trecho),
    sigma_maior_eqm = mean((sigma_maior_trecho - sigma_real)^2),
    
    .groups = "drop"
  )



tabela_phi_por_metodo <- resultados_todos %>%
  select(phi_real, starts_with("phi_")) %>%  # manter apenas colunas de phi
  pivot_longer(
    cols = -phi_real,
    names_to = "metodo",
    values_to = "phi_estimado"
  ) %>%
  mutate(
    metodo = str_replace(metodo, "phi_", "")  # extrai o nome do método
  ) %>%
  group_by(phi_real, metodo) %>%
  summarise(
    media = mean(phi_estimado),
    desvio_padrao = sd(phi_estimado),
    eqm = mean((phi_estimado - phi_real)^2),
    .groups = "drop"
  )

tabela_phi_por_metodo
