library(dplyr)
library(tidyr)
library(purrr)
library(knitr)
library(kableExtra)
library(ggplot2)


arma_with_imputations <- function(n, p, q, ar_coefs = NULL, ma_coefs = NULL, sigma = 1, tau, m, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Verificações
  if (tau + m > n) {
    stop("tau + m não pode exceder n")
  }
  
  if (is.null(ar_coefs)) {
    ar_coefs <- runif(p, -0.8, 0.8)
  }
  
  if (is.null(ma_coefs)) {
    ma_coefs <- runif(q, -0.8, 0.8)
  }
  
  # Série verdadeira
  y_true <- arima.sim(model = list(ar = ar_coefs, ma = ma_coefs), n = n, sd = sigma)
  
  # Inserir missings
  y_missing <- y_true
  missing_idx <- seq(tau + 1, tau + m)
  y_missing[missing_idx] <- NA
  
  # Imputação da média
  mean_value <- mean(y_missing, na.rm = TRUE)
  y_mean_imp <- ifelse(is.na(y_missing), mean_value, y_missing)
  
  # Hotdeck simples 
  observed_values <- y_missing[!is.na(y_missing)]
  y_hotdeck <- y_missing
  y_hotdeck[is.na(y_hotdeck)] <- sample(observed_values, sum(is.na(y_hotdeck)), replace = TRUE)
  
  # Função maior bloco
  encontra_maior_segmento <- function(x) {
    is_na <- is.na(x)
    rle_na <- rle(!is_na)
    comprimentos <- rle_na$lengths
    valores <- rle_na$values
    
    if (all(!valores)) {
      return(numeric(0))
    }
    
    idx_max <- which.max(comprimentos * valores)
    posicoes <- cumsum(comprimentos)
    
    start <- ifelse(idx_max == 1, 1, posicoes[idx_max - 1] + 1)
    end <- posicoes[idx_max]
    
    return(x[start:end])
  }
  
  maior_bloco_valores <- encontra_maior_segmento(y_missing)
  
  y_largest_block <- rep(NA, n)
  if (length(maior_bloco_valores) > 0) {
    y_largest_block[1:length(maior_bloco_valores)] <- maior_bloco_valores
  }
  
  y_collapsed_full <- rep(NA, n)
  collapsed_values <- y_missing[!is.na(y_missing)]
  y_collapsed_full[seq_along(collapsed_values)] <- collapsed_values
  
  df <- data.frame(
    t = 1:n,
    y_true = y_true,
    y_missing = y_missing,
    y_mean_imp = y_mean_imp,
    y_hotdeck = y_hotdeck,
    y_largest_block = y_largest_block,
    y_collapsed = y_collapsed_full
  )
  
  return(df)
}




simula_monte_carlo <- function(n_sim = 1000,
                               n = 200,
                               p = 1,
                               q = 1,
                               tau = 50,
                               m = 20,
                               sigma = 1) {
  
  mc_results <- purrr::map_dfr(1:n_sim, function(i) {
    
    df <- arma_with_imputations(
      n = n,
      p = p,
      q = q,
      tau = tau,
      m = m,
      sigma = sigma,
      seed = i  # <----- seed única por simulação
    )
    
    tibble(
      sim = i,
      completa = mean(df$y_true),
      media = mean(df$y_mean_imp),
      colada = mean(df$y_collapsed, na.rm = TRUE),
      hotdeck = mean(df$y_hotdeck),
      maior_bloco = mean(df$y_largest_block, na.rm = TRUE),
      n = n,
      p = p,
      q = q,
      tau = tau,
      m = m,
      sigma = sigma
    )
  })
  
  mc_results_long <- mc_results %>%
    tidyr::pivot_longer(
      cols = c(completa, media, colada, hotdeck, maior_bloco),
      names_to = "method",
      values_to = "mean_value"
    )
  
  return(mc_results_long)
}



cenarios_arma <- tibble(
  p = c(1, 0, 1),
  q = c(0, 1, 1),
  modelo = c("AR(1)", "MA(1)", "ARMA(1,1)")
)

# Definir os parâmetros (n, tau, m) conforme sua lista
parametros <- tribble(
  ~n, ~tau, ~m,
  25, 5, 5,
  25, 5, 10,
  25, 5, 15,
  25, 13, 3,
  25, 13, 6,
  25, 13, 9,
  25, 18, 2,
  25, 18, 4,
  25, 18, 5,
  50, 10, 10,
  50, 10, 20,
  50, 10, 30,
  50, 25, 6,
  50, 25, 13,
  50, 25, 19,
  50, 35, 4,
  50, 35, 8,
  50, 35, 13,
  75, 15, 15,
  75, 15, 30,
  75, 15, 45,
  75, 38, 9,
  75, 38, 19,
  75, 38, 28,
  75, 53, 5,
  75, 53, 11,
  75, 53, 16,
  100, 20, 20,
  100, 20, 40,
  100, 20, 60,
  100, 50, 13,
  100, 50, 25,
  100, 50, 38,
  100, 70, 8,
  100, 70, 15,
  100, 70, 23,
  150, 30, 30,
  150, 30, 60,
  150, 30, 90,
  150, 75, 19,
  150, 75, 38,
  150, 75, 56,
  150, 105, 11,
  150, 105, 22,
  150, 105, 33,
  200, 40, 40,
  200, 40, 80,
  200, 40, 120,
  200, 100, 25,
  200, 100, 50,
  200, 100, 75,
  200, 140, 15,
  200, 140, 30,
  200, 140, 45
)

# Número de simulações por cenário (pode ajustar)
n_sim <- 1000

# Função para rodar para um cenário ARMA + conjunto (n,tau,m)
rodar_simulacao <- function(p, q, modelo, n, tau, m, n_sim, sigma = 1) {
  sim_res <- simula_monte_carlo(
    n_sim = n_sim,
    n = n,
    p = p,
    q = q,
    tau = tau,
    m = m,
    sigma = sigma
  )
  # Não adicionar colunas que já existem no df externo!
  # sim_res <- sim_res %>%
  #   mutate(modelo = modelo, n = n, tau = tau, m = m)
  
  return(sim_res)
}

start_time <- Sys.time()

resultado_final <- parametros %>%
  crossing(cenarios_arma) %>%
  mutate(
    resultados = pmap(
      list(p, q, modelo, n, tau, m),
      rodar_simulacao,
      n_sim = n_sim
    )
  ) %>%
  select(
    modelo, resultados
  ) %>%
  unnest(resultados)

end_time <- Sys.time()

print(end_time - start_time)


# Salvando CSV ---- 

write.csv(resultado_final, 'db_sims_media.csv')


# Começando tabela do latex ---- 

db_medias <- 
  resultado_final %>%
  group_by(
    modelo,
    n,
    p,
    q,
    tau,
    m,
    sigma,
    metodo = method
  ) %>%
  summarise(
    media_simulacao = mean(mean_value),
    var_simulacao = var(mean_value),
    vies_simulacao = mean(mean_value) - 0
  ) %>%
  ungroup() %>%
  mutate(
    eqm_simulacao = var_simulacao + vies_simulacao^2
  )


# Tabela para o AR(1) ---- 

db_ar1_media <- 
  db_medias %>% 
  filter(
    modelo == 'AR(1)'
  ) %>%
  select(
    n,
    tau,
    m,
    sigma,
    metodo,
    media_simulacao,
    eqm_simulacao
  ) %>%
  pivot_wider(
    names_from = metodo,
    values_from = c(media_simulacao, eqm_simulacao)
  )


db_ar1_media %>%
  select(
    n, tau, m, starts_with("media")
  ) %>%
  mutate(across(starts_with("media"), ~ round(.x, 4))) %>%
  View()


db_ar1_media %>%
  select(
    n, tau, m, starts_with("eqm")
  ) %>%
  mutate(across(starts_with("eqm"), ~ round(.x, 4))) %>%
  View()


# Visualizações ---- 

library(ggplot2)

db_medias %>%
  glimpse()
