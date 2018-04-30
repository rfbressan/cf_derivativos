## Baseado no método de McNeil2000
## VaR e ES para NATU3 e LAME4

# Inicio ------------------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(broom)
library(gridExtra)
library(kableExtra)
library(xtable)
library(WeightedPortTest)
library(xts)
library(PerformanceAnalytics)
library(tidyquant)
source("./derivativos_fun.R") # Carrega a funcao roll_fit para fazer o backtest

# AMOSTRA COM DADOS A PARTIR DE 01-01-2006
start <- as.Date("2008-01-01")
end <- as.Date("2009-01-01")
backstart <- end + 1
w1 <- 0.5 # Peso do primeiro ativo na carteira
w <- c(w1, 1 - w1) # Vetor de pesos da carteira
u_quant <- 0.92 # quantile for treshold u
options(xtable.booktabs = TRUE) # utilizar o pacote booktabs no Latex
options(xtable.tabular.environment = "longtable") # Utiliza o pacote longtable no Latex
options(xtable.floating = FALSE) # Nao pode ser utilizado em conjunto com longtable

# XTS de precos e de retornos
assets <- c("NATU3", "LAME4")
precos <- prices(assets, start) %>% 
  group_by(symbol)
retornos <- precos %>% 
  group_by(symbol) %>% 
  tq_transmute(select = preco,
               mutate_fun = dailyReturn) 

port_ret <- retornos %>% 
  tq_portfolio(assets_col = symbol,
               returns_col = daily.returns,
               weights = w,
               col_rename = "daily.returns") %>% 
  mutate(symbol = "PORT") %>% 
  spread(key = symbol, value = daily.returns)
  
# VaR da carteira ---------------------------------------------------------

var_port <- retornos %>% 
  spread(key = symbol, value = daily.returns) %>% 
  tq_mutate(select = LAME4:NATU3,
            mutate_fun = rollapply,
            width = 252,
            FUN = VaR_port,
            by.column = FALSE,
            align = "right",
            weights = w,
            col_rename = "VaRnorm") %>% 
  tq_mutate(select = LAME4:NATU3,
            mutate_fun = rollapply,
            width = 252,
            FUN = VaR_port,
            by.column = FALSE,
            align = "right",
            method = "historical",
            weights = w,
            col_rename = "VaRhis") %>% 
  mutate(VaRnorm = Lag(-VaRnorm),
         VaRhis = Lag(-VaRhis)) %>% 
  left_join(port_ret, by = "index") %>% 
  gather(key = symbol, value = value, -index)

# Basicamente o mesmo resultado historico usando quantile
# var_his <- port_ret %>% 
#   tq_mutate(select = PORT,
#             mutate_fun = rollapply,
#             width = 252,
#             FUN = quantile,
#             align = "right",
#             probs = 0.01,
#             col_rename = "VaRnorm")


# Gera o grafico dos ativos -----------------------------------------------

paleta <- c("CF-Azul" = "#114A5C",
            "CF-Verde" = "#159399",
            "CF-Amarelo" = "#E7CA3C",
            "CF-amarelo" = "#F9DE70")
plot_precos <- ggplot(precos, aes(x = index, y = preco, group = symbol)) + 
  geom_line(aes(color = symbol)) +
  labs(title = "Evolução dos preços ajustados.",
       x = "",
       y = "Preço (R$)",
       caption = "Fonte: Yahoo!Finance") +
  scale_colour_brewer(palette = "Dark2",
                      labels = c("LAME4", "NATU3"),
                      name = "Ação") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_economist_white() 
pdf(file = "./apresentacao/figs/precos.pdf")
plot_precos
dev.off()

# Gráfico com retornos do portfolio e níveis de VaR
ret_var <- var_port %>% 
  dplyr::filter(symbol %in% c("PORT", "VaRhis", "VaRnorm")) %>% 
  dplyr::filter(!is.na(value))
plot_ret <- ggplot(ret_var, aes(x = index, y = value, group = symbol)) + 
  geom_line(aes(color = symbol)) +
  labs(title = "Retornos da carteira e VaR",
       x = "",
       y = "Retornos",
       caption = "Fonte: Yahoo!Finance e cálculos dos autores.") +
  scale_colour_brewer(palette = "Dark2",
                      labels = c("Carteira", "VaRhis", "VaRnorm"),
                      name = "") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_economist_white() 
pdf(file = "./apresentacao/figs/ret_var.pdf")
plot_ret
dev.off()
