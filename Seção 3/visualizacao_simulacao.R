library(ggplot2)
library(tidyr)
library(purrr)

arma_with_imputations <- function(n, p, q, ar_coefs = NULL, ma_coefs = NULL, sigma = 1, tau, m) {
  
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
  
  # Hotdeck simples (preencher com valores conhecidos — ex: média local ou aleatório)
  # Exemplo simples: amostragem dos valores não-missing
  set.seed(123)  # reprodutível
  observed_values <- y_missing[!is.na(y_missing)]
  y_hotdeck <- y_missing
  y_hotdeck[is.na(y_hotdeck)] <- sample(observed_values, sum(is.na(y_hotdeck)), replace = TRUE)
  
  # Sua função para maior segmento (perfeita!)
  encontra_maior_segmento <- function(x) {
    is_na <- is.na(x)
    rle_na <- rle(!is_na)
    comprimentos <- rle_na$lengths
    valores <- rle_na$values
    
    if (all(!valores)) {
      return(numeric(0)) # Tudo é NA
    }
    
    idx_max <- which.max(comprimentos * valores) # Maior bloco sem NA
    posicoes <- cumsum(comprimentos)
    
    start <- ifelse(idx_max == 1, 1, posicoes[idx_max - 1] + 1)
    end <- posicoes[idx_max]
    
    return(x[start:end])
  }
  
  # Aplicação na série
  maior_bloco_valores <- encontra_maior_segmento(y_missing)
  
  # Construir vetor y_largest_block com mesmo tamanho n
  y_largest_block <- rep(NA, n)
  if (length(maior_bloco_valores) > 0) {
    # Onde colocar o maior bloco? No início!
    y_largest_block[1:length(maior_bloco_valores)] <- maior_bloco_valores
  }
  
  # Colagem dos trechos observados (mantendo n linhas, com NA)
  y_collapsed_full <- rep(NA, n)
  collapsed_values <- y_missing[!is.na(y_missing)]
  y_collapsed_full[seq_along(collapsed_values)] <- collapsed_values
  
  # Construir o data.frame final
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


# Gerar a série
set.seed(123)
df_arma <- arma_with_imputations(n = 100, p = 2, q = 1, tau = 30, m = 8)



# Transformar para formato long
df_long <- pivot_longer(df_arma, -t, names_to = "variable", values_to = "value")

# Criar um vetor de renomeação
renomear <- c(
  y_true = "Série Completa",
  y_missing = "Série com dados faltantes",
  y_collapsed = "Série acoplada",
  y_mean_imp = "Imputação pela média",
  y_hotdeck = "Imputação Hotdeck",
  y_largest_block = "Maior bloco sem NA"
)

# Renomear a coluna 'variable'
df_long <- df_long %>%
  mutate(variable = recode(variable, !!!renomear))

# Definir ordem dos fatores (ordenando os três primeiros, o resto na sequência)
ordem <- c("Série Completa", "Série com dados faltantes", "Série acoplada")
df_long$variable <- factor(df_long$variable, levels = c(ordem, setdiff(unique(df_long$variable), ordem)))

ggplot(df_long, aes(x = t, y = value)) +
  geom_line(na.rm = TRUE, color = "steelblue", size = 1) +  # linha mais grossa
  facet_wrap(~ variable, ncol = 1, scales = "free_y", strip.position = "top") +  # painéis em 1 coluna
  labs(title = "Séries ARMA com diferentes imputações",
       x = "Tempo (t)", y = "") +
  theme_minimal(base_size = 14) +  # base menor
  theme(
    strip.text = element_text(size = 11, face = "plain", margin = margin(b = 6)),  # título do painel
    axis.text = element_text(size = 11),  # texto dos eixos
    axis.text.y = element_blank(),         # remove os valores do eixo y
    axis.ticks.y = element_blank(),        # remove os ticks do eixo y
    axis.title = element_text(size = 12),  # título dos eixos
    panel.spacing = unit(0.8, "lines"),  # espaço entre painéis
    plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),  # título centralizado e sem negrito
    plot.margin = margin(12, 12, 12, 12),  # margem do plot
    panel.grid.minor = element_blank()  # tirar grid menor
  )

