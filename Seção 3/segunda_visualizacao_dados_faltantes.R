# Carregar pacotes
library(ggplot2)
library(dplyr)

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

# Remover os NA
df_sem_na <- df %>% filter(!is.na(value))

# Reindexar o tempo se quiser "colar" as partes da série
df_sem_na <- df_sem_na %>% 
  mutate(time_reindex = row_number())

# Fazer o gráfico com a série "colada"
ggplot(df_sem_na, aes(x = time_reindex, y = value)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue") +
  labs(title = "",
       x = "Tempo(t)",
       y = "") +
  theme_minimal(base_size = 14)
