# Pacotes ---- 

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


# Autocovariâncias ---- 

fn_monte_carlo_autocov_resumo_mvar <- function(n_sim, n = 100, phi = 0.5, sigma = 1, tau = 40, m = c(20), seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  resultado_total <- tibble()
  
  # Valores teóricos das autocovariâncias
  acv1_true <- sigma^2 * phi^1 / (1 - phi^2)
  acv2_true <- sigma^2 * phi^2 / (1 - phi^2)
  acv3_true <- sigma^2 * phi^3 / (1 - phi^2)
  
  for (current_m in m) {
    for (i in 1:n_sim) {
      
      df <- fn_gera_serie_colada(n = n, phi = phi, sigma = sigma, tau = tau, m = current_m)
      
      completa <- df$completa
      imput_media <- df$imput_media
      imput_hotdeck <- df$imput_hotdeck
      
      colada <- df %>% filter(!is.na(serie_colada)) %>% pull(serie_colada)
      maior_trecho <- df %>% filter(!is.na(maior_trecho)) %>% pull(maior_trecho)
      
      auto_completa <- acf(completa, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
      auto_imput_media <- acf(imput_media, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
      auto_imput_hotdeck <- acf(imput_hotdeck, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
      auto_colada <- acf(colada, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
      auto_maior_trecho <- acf(maior_trecho, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
      
      tabela_sim <- tibble(
        sim = i,
        m = current_m,
        serie = rep(c("completa", "imput_media", "imput_hotdeck", "colada", "maior_trecho"), each = 3),
        lag = rep(1:3, times = 5),
        autocov = c(auto_completa[2:4],
                    auto_imput_media[2:4],
                    auto_imput_hotdeck[2:4],
                    auto_colada[2:4],
                    auto_maior_trecho[2:4])
      )
      
      resultado_total <- bind_rows(resultado_total, tabela_sim)
    }
  }
  
  resumo_estatisticas <- resultado_total %>%
    group_by(m, serie, lag) %>%
    summarise(
      media_autocov = mean(autocov, na.rm = TRUE),
      var_autocov = var(autocov, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      acv_true = case_when(
        lag == 1 ~ acv1_true,
        lag == 2 ~ acv2_true,
        lag == 3 ~ acv3_true
      ),
      bias_autocov = media_autocov - acv_true,
      eqm_autocov = bias_autocov^2 + var_autocov
    ) %>%
    select(m, serie, lag, media_autocov, var_autocov, acv_true, bias_autocov, eqm_autocov)
  
  return(list(
    resultado_total = resultado_total,
    resumo_estatisticas = resumo_estatisticas
  ))
}



# Rodando cenários ---- 



# Cenário n = 25
res_n25 <- fn_monte_carlo_autocov_resumo_mvar(
  n_sim = 1000,
  n = 25,
  phi = 0.5,
  sigma = 1,
  tau = 10,
  m = c(7, 10, 13),
  seed = 123
)

# Cenário n = 50
res_n50 <- fn_monte_carlo_autocov_resumo_mvar(
  n_sim = 1000,
  n = 50,
  phi = 0.5,
  sigma = 1,
  tau = 20,
  m = c(14, 20, 26),
  seed = 123
)

# Cenário n = 75
res_n75 <- fn_monte_carlo_autocov_resumo_mvar(
  n_sim = 1000,
  n = 75,
  phi = 0.5,
  sigma = 1,
  tau = 30,
  m = c(21, 30, 39),
  seed = 123
)

# Cenário n = 100
res_n100 <- fn_monte_carlo_autocov_resumo_mvar(
  n_sim = 1000,
  n = 100,
  phi = 0.5,
  sigma = 1,
  tau = 40,
  m = c(28, 40, 52),
  seed = 123
)

# Cenário n = 250
res_n250 <- fn_monte_carlo_autocov_resumo_mvar(
  n_sim = 1000,
  n = 250,
  phi = 0.5,
  sigma = 1,
  tau = 100,
  m = c(70, 100, 130),
  seed = 123
)


todos_resumos <- bind_rows(
  res_n25$resumo_estatisticas %>% mutate(n = 25),
  res_n50$resumo_estatisticas %>% mutate(n = 50),
  res_n75$resumo_estatisticas %>% mutate(n = 75),
  res_n100$resumo_estatisticas %>% mutate(n = 100),
  res_n250$resumo_estatisticas %>% mutate(n = 250)
) %>%
  select(n, everything())  


View(todos_resumos)


write.csv2(todos_resumos, 'resumos_autocov.csv')
