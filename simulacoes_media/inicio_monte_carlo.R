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



fn_monte_carlo <- function(n_sim, n = 100, phi = 0.5, sigma = 1, tau = 40, m = 20, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  resultado_total <- tibble()
  
  for (i in 1:n_sim) {
    
    # 1. Gera o df 
    df <- fn_gera_serie_colada(n = n, phi = phi, sigma = sigma, tau = tau, m = m)
    
    # 2. Extrai vetores
    completa <- df$completa
    imput_media <- df$imput_media
    imput_hotdeck <- df$imput_hotdeck
    
    colada <- df %>% filter(!is.na(serie_colada)) %>% pull(serie_colada)
    maior_trecho <- df %>% filter(!is.na(maior_trecho)) %>% pull(maior_trecho)
    
    # 3. Estatísticas
    medias <- tibble(
      serie = c("completa", "imput_media", "imput_hotdeck", "colada", "maior_trecho"),
      media = c(mean(completa, na.rm = TRUE),
                mean(imput_media, na.rm = TRUE),
                mean(imput_hotdeck, na.rm = TRUE),
                mean(colada),
                mean(maior_trecho))
    )
    
    variancias <- tibble(
      serie = c("completa", "imput_media", "imput_hotdeck", "colada", "maior_trecho"),
      variancia = c(var(completa, na.rm = TRUE),
                    var(imput_media, na.rm = TRUE),
                    var(imput_hotdeck, na.rm = TRUE),
                    var(colada),
                    var(maior_trecho))
    )
    
    # Autocovariâncias para lag 1, 2 e 3
    auto_completa <- acf(completa, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
    auto_imput_media <- acf(imput_media, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
    auto_imput_hotdeck <- acf(imput_hotdeck, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
    auto_colada <- acf(colada, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
    auto_maior_trecho <- acf(maior_trecho, lag.max = 3, plot = FALSE, type = "covariance")$acf[,1,1]
    
    autocovs1 <- tibble(
      serie = c("completa", "imput_media", "imput_hotdeck", "colada", "maior_trecho"),
      autocov_lag1 = c(auto_completa[2], auto_imput_media[2], auto_imput_hotdeck[2], auto_colada[2], auto_maior_trecho[2])
    )
    
    autocovs2 <- tibble(
      serie = c("completa", "imput_media", "imput_hotdeck", "colada", "maior_trecho"),
      autocov_lag2 = c(auto_completa[3], auto_imput_media[3], auto_imput_hotdeck[3], auto_colada[3], auto_maior_trecho[3])
    )
    
    autocovs3 <- tibble(
      serie = c("completa", "imput_media", "imput_hotdeck", "colada", "maior_trecho"),
      autocov_lag3 = c(auto_completa[4], auto_imput_media[4], auto_imput_hotdeck[4], auto_colada[4], auto_maior_trecho[4])
    )
    
    # Junta tudo
    tabela_final <- medias %>%
      left_join(variancias, by = "serie") %>%
      left_join(autocovs1, by = "serie") %>%
      left_join(autocovs2, by = "serie") %>%
      left_join(autocovs3, by = "serie") %>%
      mutate(sim = i) %>%
      select(sim, everything())
    
    # Acumula o resultado
    resultado_total <- bind_rows(resultado_total, tabela_final)
  }
  
  return(resultado_total)
}

tabela_mc <- fn_monte_carlo(n_sim = 1000)


# Valores verdadeiros do modelo (ajuste conforme necessário)
media_true <- 0
phi <- 0.5
variancia_true <- 1/(1-(0.5)^2)
sigma <- 1
acv1_true <- sigma^2 * phi^1 / (1 - phi^2)
acv2_true <- sigma^2 * phi^2 / (1 - phi^2)
acv3_true <- sigma^2 * phi^3 / (1 - phi^2)

# Estatísticas Monte Carlo
teste_completo <- 
  tabela_mc %>%
  group_by(serie) %>%
  summarise(
    # Média
    media_media = mean(media, na.rm = TRUE),
    var_media = var(media, na.rm = TRUE),
    bias_media = media_media - media_true,
    eqm_media = bias_media^2 + var_media,
    
    # Variância
    media_variancia = mean(variancia, na.rm = TRUE),
    var_variancia = var(variancia, na.rm = TRUE),
    bias_variancia = media_variancia - variancia_true,
    eqm_variancia = bias_variancia^2 + var_variancia,
    
    # Autocovariância lag 1
    media_autocov_lag1 = mean(autocov_lag1, na.rm = TRUE),
    var_autocov_lag1 = var(autocov_lag1, na.rm = TRUE),
    bias_autocov_lag1 = media_autocov_lag1 - acv1_true,
    eqm_autocov_lag1 = bias_autocov_lag1^2 + var_autocov_lag1,
    
    # Autocovariância lag 2
    media_autocov_lag2 = mean(autocov_lag2, na.rm = TRUE),
    var_autocov_lag2 = var(autocov_lag2, na.rm = TRUE),
    bias_autocov_lag2 = media_autocov_lag2 - acv2_true,
    eqm_autocov_lag2 = bias_autocov_lag2^2 + var_autocov_lag2,
    
    # Autocovariância lag 3
    media_autocov_lag3 = mean(autocov_lag3, na.rm = TRUE),
    var_autocov_lag3 = var(autocov_lag3, na.rm = TRUE),
    bias_autocov_lag3 = media_autocov_lag3 - acv3_true,
    eqm_autocov_lag3 = bias_autocov_lag3^2 + var_autocov_lag3
  )




# Pensando em visualização --- 



# ggplot(tabela_mc, aes(x = serie, y = media)) +
#   geom_boxplot(fill = "skyblue") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # valor real da média
#   labs(title = "Distribuição da Média Estimada nas Simulações",
#        y = "Média estimada", x = "Série")
# 
# 
# ggplot(teste_completo, aes(x = serie, y = bias_media)) +
#   geom_col(fill = "orange") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(title = "Viés das Médias Estimadas", y = "Viés", x = "Série")
# 
# eqm_resumo <- teste_completo %>%
#   select(serie, eqm_media, eqm_variancia, eqm_autocov_lag1, eqm_autocov_lag2, eqm_autocov_lag3) %>%
#   pivot_longer(-serie, names_to = "estatistica", values_to = "eqm")
# 
# ggplot(eqm_resumo, aes(x = serie, y = eqm, fill = estatistica)) +
#   geom_col(position = position_dodge()) +
#   labs(title = "Erro Quadrático Médio dos Estimadores", y = "EQM", x = "Série")
# 
# ggplot(tabela_mc, aes(x = sim, y = media, color = serie)) +
#   geom_line(alpha = 0.7) +
#   labs(title = "Médias Estimadas ao Longo das Simulações", x = "Número da Simulação", y = "Média Estimada")


teste_media_variancia <- 
  teste_completo %>%
  select(
    serie,
    media_media, var_media, bias_media, eqm_media,
    media_variancia, var_variancia, bias_variancia, eqm_variancia
  )
