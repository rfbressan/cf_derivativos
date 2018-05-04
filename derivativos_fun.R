## Funcoes utilizadas em derivativos.R
## artigo_fun.R

# Carrega os pacotes necessarios se faltantes
library(timetk)

# prices ------------------------------------------------------------------
# Extrai a serie de precos ajustados dos ativos em analise
prices <- function(assets, start) {
  n <- length(assets)
  l_prices <- vector(mode = "list", length = n)
  for (i in seq_len(n)) {
    tb <- read_csv(paste0("./input/cf-", assets[i], ".csv"), 
                   col_types = cols_only(Date = col_date(), 
                                         `Adj Close` = col_double()))
    l_prices[[i]] <- xts(tb$`Adj Close`, order.by = tb$Date)[paste0(start, "/")]
    colnames(l_prices[[i]]) <- assets[i]
  }
  xts_prices <- na.locf(do.call(merge, l_prices))
  tbl_prices <- tk_tbl(xts_prices) %>% 
    gather(key = symbol, value = preco, -index)
  return(tbl_prices)
}

# VaR_port ----------------------------------------------------------------

VaR_port <- function(R, 
                     p = 0.99,
                     method = "gaussian",
                     portfolio_method = "component",
                     weights = NULL,
                     invert = TRUE){
  var_p <- PerformanceAnalytics::VaR(R = R,
                                     p = p,
                                     method = method,
                                     portfolio_method = portfolio_method,
                                     weights = weights,
                                     invert = invert)
  if (method == "gaussian")
    return(var_p$VaR)
  else if (method == "historical")
    return(var_p)
}

# ES_port ----------------------------------------------------------------

ES_port <- function(R, 
                     p = 0.99,
                     method = "gaussian",
                     portfolio_method = "component",
                     weights = NULL,
                     invert = TRUE){
  es_p <- PerformanceAnalytics::ETL(R = R,
                                     p = p,
                                     method = method,
                                     portfolio_method = portfolio_method,
                                     weights = weights,
                                     invert = invert)
  if (method == "gaussian")
    return(es_p$ES)
  else if (method == "historical")
    return(es_p)
}

