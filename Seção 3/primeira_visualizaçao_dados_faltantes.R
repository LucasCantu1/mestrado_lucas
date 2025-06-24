# Carregar pacotes
library(ggplot2)

# Definir semente para reprodutibilidade
set.seed(123)

# Criar série temporal
n <- 100
time <- 1:n
value <- cumsum(rnorm(n))  # Série simulada com tendência

# Introduzir 30 NA a partir da observação 51
value[51:(51 + 30 - 1)] <- NA

# Colocar em um data frame
df <- data.frame(time = time, value = value)

# Fazer o gráfico (sem indicar os NAs)
ggplot(df, aes(x = time, y = value)) +
  geom_line(color = "steelblue", size = 1) +     # Apenas a linha, com quebras nos NAs
  geom_point(color = "steelblue") +              # Pontos para os valores não faltantes
  labs(title = "",
       x = "Tempo (t)",
       y = "") +
  theme_minimal(base_size = 14)



