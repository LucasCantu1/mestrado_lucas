# Pacotes ----- 
{
  library(ggplot2)
  library(tidyr)
  library(dplyr)
}

# Função Principal ---- 

fn_gera_serie_colada <- function(n, phi = 0.8, sigma = 1, tau = 10, m = 5, seed = NULL) {
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
  y_incompleta[tau:(tau + m - 1)] <- NA
  
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
  
  # Dataframe final
  df <- data.frame(
    tempo = 1:n,
    completa = y,
    incompleta = y_incompleta,
    imput_media = y_media,
    imput_hotdeck = y_hotdeck,
    serie_colada = y_colada
  )
  
  return(df)
}

# Plot para o trabalho ----

# Extração dos parâmetros da série
tau_sim <- 40
m_sim <- 30
n_sim <- 100 
phi_sim = 0.5
sigma_sim = 1
seed_sim = 123456

# Gerar dados
df <- fn_gera_serie_colada(n = n_sim, phi = phi_sim, sigma = sigma_sim, tau = tau_sim, m = m_sim, seed = seed_sim)

# Determinar os dois trechos possíveis
idx_antes <- 1:(tau_sim - 1)
idx_depois <- (tau_sim + m_sim):n_sim

# Seleciona o maior trecho
trecho_maior <- if (length(idx_antes) >= length(idx_depois)) {
  df[idx_antes, c("tempo", "completa")] %>% mutate(tipo = "Maior trecho observado")
} else {
  df[idx_depois, c("tempo", "completa")] %>% mutate(tipo = "Maior trecho observado")
}
colnames(trecho_maior)[2] <- "valor"

# Transformar o restante do df para formato longo
df_long <- df %>%
  pivot_longer(cols = -tempo, names_to = "tipo", values_to = "valor") %>%
  mutate(tipo = recode(tipo,
                       completa = "Série completa",
                       incompleta = "Série com falhas",
                       imput_media = "Imputação: média",
                       imput_hotdeck = "Imputação: hot-deck",
                       serie_colada = "Série acoplada (sem NA)"
  ))

# Adicionar o quinto cenário (maior trecho)
df_final <- bind_rows(df_long, trecho_maior)

# Definir ordem dos gráficos
df_final$tipo <- factor(df_final$tipo, levels = c(
  "Série completa",
  "Série com falhas",
  "Série acoplada (sem NA)",
  "Imputação: média",
  "Maior trecho observado",
  "Imputação: hot-deck"
))

# Plot final
# Plot final sem eixo y
ggplot(df_final, aes(x = tempo, y = valor)) +
  geom_line(color = "steelblue", size = 1, na.rm = FALSE) +
  facet_wrap(~ tipo, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    title = "",
    x = "Tempo (t)",
    y = ""
  )
