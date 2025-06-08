# Gerando as simulações para a média 


library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

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



# Valores verdadeiros
media_true <- 0
phi <- 0.5
sigma <- 1

# Função Monte Carlo adaptada apenas para a média
fn_monte_carlo_media <- function(n_sim, n = 25, phi = 0.5, sigma = 1, tau = 10, m, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  resultado_total <- tibble()
  
  for (i in 1:n_sim) {
    
    df <- fn_gera_serie_colada(n = n, phi = phi, sigma = sigma, tau = tau, m = m)
    
    completa <- df$completa
    imput_media <- df$imput_media
    imput_hotdeck <- df$imput_hotdeck
    colada <- df %>% filter(!is.na(serie_colada)) %>% pull(serie_colada)
    maior_trecho <- df %>% filter(!is.na(maior_trecho)) %>% pull(maior_trecho)
    
    medias <- tibble(
      serie = c("completa", "imput_media", "imput_hotdeck", "colada", "maior_trecho"),
      media = c(mean(completa),
                mean(imput_media),
                mean(imput_hotdeck),
                mean(colada),
                mean(maior_trecho)),
      sim = i
    )
    
    resultado_total <- bind_rows(resultado_total, medias)
  }
  
  return(resultado_total)
}




# Aqui eu rodei alguns pontos diferentes para fazer a tabela, mas a função parece
# estar fazendo o esperado, podemos definir melhor os pontos 

# Rodar para diferentes valores de m
valores_m <- c(10, 50, 75)

# Valores verdadeiros
media_true <- 0
phi <- 0.5
sigma <- 1


resultados_m <- map_dfr(valores_m, ~{
  df <- fn_monte_carlo_media(n_sim = 1000, n = 100, phi = 0.5, sigma = 1, tau = 20, m = .x)
  df$m <- .x
  df
})

# Tabela final resumida para as médias
tabela_resumo_medias <- resultados_m %>%
  group_by(m, serie) %>%
  summarise(
    media_media = mean(media),
    var_media = var(media),
    bias_media = media_media - media_true,
    eqm_media = bias_media^2 + var_media,
    .groups = "drop"
  )

print(tabela_resumo_medias)
