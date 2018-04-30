## Funcoes utilizadas em derivativos.R
## artigo_fun.R

# Carrega os pacotes necessarios se faltantes
library(timetk)

cores <- detectCores() # Quantos cores estao rodando

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

# roll_fit ----------------------------------------------------------------

roll_fit <- function(data, window.size, n.roll, spec, models) {
  cat("\nroll_fit do ativo inicio:", as.character(Sys.time()))
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll), is.null(spec), is.null(models)))
    stop("roll_fit: Devem ser passados todos os argumentos!")
  if(class(spec) != "uGARCHspec") stop("roll_fit: spec deve ser da classe uGARCspec")
  
  garch_fit.list <- mclapply(1:n.roll, function(i){
    garch.fit <- try(ugarchfit(spec,
                               data[i:(window.size+i)],
                               solver = "hybrid"))
    return(garch.fit)
  },
  mc.cores = cores) # Fim do mclapply
  # Indice para ordenar as respostas das funcoes roll_fit_X
  id.xts <- index(data[(window.size+2):(window.size+n.roll)])
  
  tmp.list <- map(models, 
                  ~switch (.x,
                           cevt = roll_fit_cevt(garch_fit.list, id.xts),
                           cnorm = roll_fit_cnorm(garch_fit.list, id.xts),
                           ct = roll_fit_ct(garch_fit.list, id.xts),
                           uevt = roll_fit_uevt(garch_fit.list, id.xts),
                           unorm = roll_fit_unorm(data, window.size, n.roll),
                           ut = roll_fit_ut(data, window.size, n.roll),
                           riskmetrics = roll_fit_riskmetrics(data, window.size, n.roll),
                           cat("Modelo não definido:", .x)
                           ) # fim switch
                  ) # fim map
  names(tmp.list) <- models
  tmp.list <- enframe(tmp.list)
  colnames(tmp.list) <- c("model_type", "roll")
  tmp.list <- subset(tmp.list, subset = !is.null(tmp.list$roll))
  cat("\nroll_fit do ativo fim:", as.character(Sys.time()), "\n")
  return(tmp.list)
}
# roll_fit_cevt -----------------------------------------------------------
roll_fit_cevt <- function(garch_list, id_xts) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!length(garch_list)) stop("roll_fit_cevt: argumento de tamanho invalido!")
  # Extrair residuos padronizados
  # Aplicar o gpfFit
  # Calcular zq e sq
  # Para cada interacao. n.roll vezes
  # Devolver tudo em um data.frame
  # data deve ser o xts com as perdas de todo o periodo
  tmp.list <- mclapply(seq_along(garch_list), function(i){
    # 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
    if(inherits(garch_list[[i]], 'try-error') || 
       convergence(garch_list[[i]])!=0 || 
       is.null(garch_list[[i]]@fit$cvar)){
      lapply_ans <- t(cbind(rep(NA, 9)))
      # Se algum erro, envia como resposta NA, mas nao interrompe a estimacao
      # Depois no xts retornado pela funcao roll_fit, preencher os NA com os dados
      # da estimacao anterior com na.locf
    } else{
      # Salvar o coeficiente beta1 apenas para mostrar a evolucao deste ao longo
      # do período
      beta1 <- coef(garch_list[[i]])["beta1"]
      # Retirar os residuos padronizados, e o ultimo mu_t e sigma_t
      resid_z <- coredata(residuals(garch_list[[i]], standardize = TRUE))
      mu_t <- last(fitted(garch_list[[i]]))
      sigma_t <- last(sigma(garch_list[[i]]))
      # Ajustar uma gpd aos dados de residuos padronizados
      gpd.fit <- gpdFit(resid_z, u = quantile(resid_z, u_quant))
      xi = gpd.fit@fit$par.ests[1]
      beta = gpd.fit@fit$par.ests[2]
      xi_se = gpd.fit@fit$par.ses[1]
      beta_se = gpd.fit@fit$par.ses[2]
      # Obter as medias de risco do residuo padronizado
      zq975 <- gpdRiskMeasures(gpd.fit, prob = 0.975)$quantile
      zq990 <- gpdRiskMeasures(gpd.fit, prob = 0.990)$quantile
      sq975 <- gpdRiskMeasures(gpd.fit, prob = 0.975)$shortfall
      sq990 <- gpdRiskMeasures(gpd.fit, prob = 0.990)$shortfall
      # Agora escalar e deslocar as medidas de risco de acordo com o modelo garch
      # Zq = mu_t+sqrt(sigma_t)*zq
      # Sq = mu_t+sqrt(sigma_t)*sq
      Zq975 <- mu_t+sigma_t*zq975
      Zq990 <- mu_t+sigma_t*zq990
      Sq975 <- mu_t+sigma_t*sq975
      Sq990 <- mu_t+sigma_t*sq990
      
      # Verbose mode
      #cat("Estimacao numero:", i, "as", as.character(Sys.time()))
      lapply_ans <- cbind(beta1, xi, beta, xi_se, beta_se, Zq975, Zq990, Sq975, Sq990)
    } # fim do else
    return(lapply_ans)
  },
  mc.cores = cores) # fim do mclapply
  #stopCluster(cluster)
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]
  ans <- xts(ans, order.by = id_xts)
  colnames(ans) <- c("beta1", "xi", "beta", "xi_se", "beta_se", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("beta1", "xi", "beta", "xi_se", "beta_se")]
  toc <- Sys.time()
  cat("\ncevt:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit_cevt

# roll_fit_cnorm -----------------------------------------------------------
roll_fit_cnorm <- function(garch_list, id_xts) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!length(garch_list)) stop("roll_fit_cnorm: argumento de tamanho invalido!")
  
  tmp.list <- mclapply(seq_along(garch_list), function(i){
    # 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
    if(inherits(garch_list[[i]], 'try-error') || 
       convergence(garch_list[[i]])!=0 || 
       is.null(garch_list[[i]]@fit$cvar)){
      lapply_ans <- t(cbind(rep(NA, 5)))
      # Se algum erro, envia como resposta NA, mas nao interrompe a estimacao
      # Depois no xts retornado pela funcao roll_fit, preencher os NA com os dados
      # da estimacao anterior com na.locf
    } else{
      # Retirar os residuos padronizados, e o ultimo mu_t e sigma_t
      resid_z <- coredata(residuals(garch_list[[i]], standardize = TRUE))
      mu_t <- last(fitted(garch_list[[i]]))
      sigma_t <- last(sigma(garch_list[[i]]))
      # Ajustar os dados de residuos padronizados
      fit <- fitdist(distribution = "norm", resid_z)
      # Obter os quantis do residuo padronizado
      zq975 <- qdist(distribution = "norm", mu = fit$pars["mu"], sigma = fit$pars["sigma"], p = 0.975)
      zq990 <- qdist(distribution = "norm", mu = fit$pars["mu"], sigma = fit$pars["sigma"], p = 0.990)
      # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
      # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
      sq975 <- integrate(function(x){
        qnorm(x, mean = fit$pars["mu"], sd = fit$pars["sigma"])},
        0.975,
        1)$value / (1-0.975)
      sq990 <- integrate(function(x){
        qnorm(x, mean = fit$pars["mu"], sd = fit$pars["sigma"])},
        0.990,
        1)$value / (1-0.990)
      
      # Agora escalar e deslocar as medidas de risco de acordo com o modelo garch
      # Zq = mu_t+sqrt(sigma_t)*zq
      # Sq = mu_t+sqrt(sigma_t)*sq
      Zq975 <- mu_t+sigma_t*zq975
      Zq990 <- mu_t+sigma_t*zq990
      Sq975 <- mu_t+sigma_t*sq975
      Sq990 <- mu_t+sigma_t*sq990

      lapply_ans <- cbind(fit$pars["sigma"], Zq975, Zq990, Sq975, Sq990)
    } # fim do else
    return(lapply_ans)
  },
  mc.cores = cores) # fim do mclapply
  #stopCluster(cluster)
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]
  ans <- xts(ans, order.by = id_xts)
  colnames(ans) <- c("sigma", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("sigma")]
  toc <- Sys.time()
  cat("\ncnorm:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
} # fim da roll_fit_cnorm

# roll_fit_ct -----------------------------------------------------------
roll_fit_ct <- function(garch_list, id_xts) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!length(garch_list)) stop("roll_fit_ct: argumento de tamanho invalido!")
  
    tmp.list <- mclapply(seq_along(garch_list), function(i){
    # 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
    if(inherits(garch_list[[i]], 'try-error') || 
       convergence(garch_list[[i]])!=0 || 
       is.null(garch_list[[i]]@fit$cvar)){
      lapply_ans <- t(cbind(rep(NA, 6)))
      # Se algum erro, envia como resposta NA, mas nao interrompe a estimacao
      # Depois no xts retornado pela funcao roll_fit, preencher os NA com os dados
      # da estimacao anterior com na.locf
    } else{
      # Retirar os residuos padronizados, e o ultimo mu_t e sigma_t
      resid_z <- coredata(residuals(garch_list[[i]], standardize = TRUE))
      mu_t <- last(fitted(garch_list[[i]]))
      sigma_t <- last(sigma(garch_list[[i]]))
      # Ajustar os dados de residuos padronizados
      fit <- fitdist(distribution = "std", resid_z)
      # Obter os quantis do residuo padronizado
      zq975 <- qdist(distribution = "std", 
                     mu = fit$pars["mu"], sigma = fit$pars["sigma"], shape = fit$pars["shape"],
                     p = 0.975)
      zq990 <- qdist(distribution = "norm", 
                     mu = fit$pars["mu"], sigma = fit$pars["sigma"], shape = fit$pars["shape"],
                     p = 0.990)
      # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
      # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
      # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
      # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
      sq975 <- integrate(function(x){
        qdist(distribution = "std", x, 
              mu = fit$pars["mu"], sigma = fit$pars["sigma"], shape = fit$pars["shape"])},
        0.975,
        1)$value / (1-0.975)
      sq990 <- integrate(function(x){
        qdist(distribution = "std", x, 
              mu = fit$pars["mu"], sigma = fit$pars["sigma"], shape = fit$pars["shape"])},
        0.990,
        1)$value / (1-0.990)
      
      # Agora escalar e deslocar as medidas de risco de acordo com o modelo garch
      # Zq = mu_t+sqrt(sigma_t)*zq
      # Sq = mu_t+sqrt(sigma_t)*sq
      Zq975 <- mu_t+sigma_t*zq975
      Zq990 <- mu_t+sigma_t*zq990
      Sq975 <- mu_t+sigma_t*sq975
      Sq990 <- mu_t+sigma_t*sq990
      
      lapply_ans <- cbind(fit$pars["sigma"], fit$pars["shape"], Zq975, Zq990, Sq975, Sq990)
    } # fim do else
    return(lapply_ans)
  },
  mc.cores = cores) # fim do mclapply
  #stopCluster(cluster)
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]
  ans <- xts(ans, order.by = id_xts)
  colnames(ans) <- c("sigma", "shape", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("sigma", "shape")]
  toc <- Sys.time()
  cat("\nct:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit_ct

# roll_fit_unorm ----------------------------------------------------------
# Ajusta os dados para um modelo Normal incondicional
roll_fit_unorm <- function(data, window.size, n.roll) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_unorm: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll)))
    stop("roll_fit_unorm: Devem ser passados todos os argumentos!")
  
  tmp.list <- mclapply(1:n.roll, function(i){
    xts <- data[i:(window.size+i)]
    mean <- mean(xts)
    sd <- sd(xts)
    Zq975 <- qnorm(0.975, mean, sd)
    Zq990 <- qnorm(0.990, mean, sd)
    # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
    # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
    Sq975 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.975,
      1)$value / (1-0.975)
    Sq990 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.990,
      1)$value / (1-0.990)
    return(cbind(mean, sd, Zq975, Zq990, Sq975, Sq990))
  },
  mc.cores = cores) # Fim do lapply
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("mu", "sigma", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("mu", "sigma")]
  toc <- Sys.time()
  cat("\nunorm:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
}

# roll_fit_ut ----------------------------------------------------------
# Ajusta os dados para um modelo t-Student incondicional
roll_fit_ut <- function(data, window.size, n.roll) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_ut: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll)))
    stop("roll_fit_ut: Devem ser passados todos os argumentos!")
  
  tmp.list <- mclapply(1:n.roll, function(i){
    xts <- data[i:(window.size+i)]
    t_fit <- fitdist(distribution = "std", xts)
    t_mu <- t_fit$pars["mu"]
    t_sigma <- t_fit$pars["sigma"]
    t_shape <- t_fit$pars["shape"]
    Zq975 <- qdist(distribution = "std", p = 0.975,
                   mu = t_mu,
                   sigma = t_sigma,
                   shape = t_shape)
    Zq990 <- qdist(distribution = "std", p = 0.990,
                   mu = t_mu,
                   sigma = t_sigma,
                   shape = t_shape)
    # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
    # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
    Sq975 <- integrate(function(x){
      qdist(distribution = "std", x, mu = t_mu, sigma = t_sigma, shape = t_shape)},
      0.975,
      1)$value / (1-0.975)
    Sq990 <- integrate(function(x){
      qdist(distribution = "std", x, mu = t_mu, sigma = t_sigma, shape = t_shape)},
      0.990,
      1)$value / (1-0.990)
    return(cbind(t_mu, t_sigma, t_shape, Zq975, Zq990, Sq975, Sq990))
  },
  mc.cores = cores) # Fim do mclapply
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("mu", "sigma", "shape", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("mu", "sigma", "shape")]
  toc <- Sys.time()
  cat("\nut:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
}

# roll_fit_uevt -----------------------------------------------------------
roll_fit_uevt <- function(garch_list, id_xts) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!length(garch_list)) stop("roll_fit_uevt: argumento de tamanho invalido!")
  # Extrair residuos padronizados
  # Aplicar o gpfFit
  # Calcular zq e sq com base na media e variancia INCONDICIONAIS
  # Para cada interacao. n.roll vezes
  # Devolver tudo em um data.frame (ou lista)
  # data deve ser o xts com as perdas de todo o periodo
  tmp.list <- mclapply(seq_along(garch_list), function(i){
    # 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
    if(inherits(garch_list[[i]], 'try-error') || 
       convergence(garch_list[[i]])!=0 || 
       is.null(garch_list[[i]]@fit$cvar)){
      lapply_ans <- t(cbind(rep(NA, 6)))
      # Se algum erro, envia como resposta NA, mas nao interrompe a estimacao
      # Depois no xts retornado pela funcao roll_fit, preencher os NA com os dados 
      # da estimacao anterior com na.locf
    } else{
      # Retirar os residuos padronizados, e o ultimo mu_t e sigma_t
      resid_z <- coredata(residuals(garch_list[[i]], standardize = TRUE))
      mu_t <- uncmean(garch_list[[i]])
      sigma_t <- sqrt(uncvariance(garch_list[[i]])) 
      # Ajustar uma gpd aos dados de residuos padronizados
      gpd.fit <- gpdFit(resid_z, u = quantile(resid_z, u_quant))
      xi = gpd.fit@fit$par.ests[1]
      beta = gpd.fit@fit$par.ests[2]
      # xi_se = gpd.fit@fit$par.ses[1]
      # beta_se = gpd.fit@fit$par.ses[2]
      # Obter as medias de risco do residuo padronizado
      zq975 <- gpdRiskMeasures(gpd.fit, prob = 0.975)$quantile
      zq990 <- gpdRiskMeasures(gpd.fit, prob = 0.990)$quantile
      sq975 <- gpdRiskMeasures(gpd.fit, prob = 0.975)$shortfall
      sq990 <- gpdRiskMeasures(gpd.fit, prob = 0.990)$shortfall
      # Agora escalar e deslocar as medidas de risco de acordo com o modelo garch
      # Zq = mu_t+sigma_t*zq
      # Sq = mu_t+sigma_t*sq
      Zq975 <- mu_t+sigma_t*zq975
      Zq990 <- mu_t+sigma_t*zq990
      Sq975 <- mu_t+sigma_t*sq975
      Sq990 <- mu_t+sigma_t*sq990
      
      # Verbose mode
      #cat("Estimacao numero:", i, "as", as.character(Sys.time()))
      lapply_ans <- cbind(xi, beta, Zq975, Zq990, Sq975, Sq990)
    } # fim do else
    return(lapply_ans)
  },
  mc.cores = cores) # fim do mclapply
  #stopCluster(cluster)
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = id_xts)
  colnames(ans) <- c("xi", "beta", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("xi", "beta")]
  toc <- Sys.time()
  cat("\nuevt:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit_uevt

# roll_fit_riskmetrics -----------------------------------------------------------
roll_fit_riskmetrics <- function(data, window.size, n.roll){
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_riskmetrics: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll)))
    stop("roll_fit_riskmetrics: Devem ser passados todos os argumentos!")
  # Eh uma especificacao de Garch(1,1) com mu = 0, omega = 0, lambda = beta1 e 1-lambda = alpha1
  # Como os parametros do Garch sao fixos, pode-se utilizar o metodo ugarchfilter
  lambda <- 0.94 # Valor apontado como ideal para dados diarios
  ruspec <- ugarchspec(mean.model = list(armaOrder = c(0,0),
                                         include.mean = FALSE),
                       variance.model = list(model = "iGARCH",
                                             garchOrder = c(1, 1)),
                       distribution.model = "norm",
                       fixed.pars = list(omega = 0,
                                         alpha1 = (1-lambda))) # Beta eh calculado no modelo iGarch
  filter <- ugarchfilter(ruspec, data[(window.size+1):(window.size+n.roll)])
  resid <- residuals(filter)
  sigma <- sigma(filter) # variancia!! eh necessario tirar a raiz para obter o desv. padrao
  Zq975 <- sigma*qnorm(0.975)
  Zq990 <- sigma*qnorm(0.990)
  Sq975 <- (coredata(sigma)*dnorm(qnorm(0.975)))/0.025 # Eq 4.7 p. 37 de Pfaff2013
  Sq990 <- (coredata(sigma)*dnorm(qnorm(0.990)))/0.01
  
  ans <- cbind(lambda, Zq975, Zq990, Sq975, Sq990) # 
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("lambda", "Zq975", "Zq990", "Sq975", "Sq990") # 
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("lambda")]
  toc <- Sys.time()
  cat("\nriskmetrics:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
}

# es_test -----------------------------------------------------------------
# Teste nao parametrico para os residuos das violacoes ao VaR, conforme 
# teste de ES de MacNeil2000
# ATENCAO!! Este teste so faz sentido ser realizados apos os teste de VaR
# aprovarem o modelo
es_test <- function(alpha = 0.0275, losses, ES, VaR, conf.level = 0.95, n.boot = 1000) {
  # Check for univariate series
  if(!all(dim(losses)[2] == 1,
          dim(ES)[2] == 1,
          dim(VaR)[2] == 1))
    stop("\nSeries are not univariate!")
  # Check for the same length in all series
  N <-  length(losses)
  if(!all(length(VaR) == N,
          length(ES) == N))
    stop("\nLength of series are not equal!")
  idx <-  which(losses > VaR)
  s <-  ES[idx] # ES values when VaR violation occurred
  x <- losses[idx] # losses that violated VaR
  
  # One-sided test. H0: mean of x is less than or equal to mean of s
  # We do not want to reject H0
  boot <- boot.t.test(x, s, reps = n.boot, alternative = "greater")
  
  ans <-  tibble(expected.exceed = floor(alpha*N),
                 actual.exceed = length(idx)) %>% 
    bind_cols(boot)
  # conditional expected shortfall is systematically underestimated
  ans$H0 <- "Mean of Excess Violations of VaR is less than or equal to zero"
  ans$Decision <- ifelse(ans$p.value < (1 - conf.level), "Reject H0", "Failed to reject H0")
  return(ans)
} # end of es_test

# boot.t.test -------------------------------------------------------------
# Codigo importado de https://github.com/tpepler/nonpar
# Adaptado apenas para as checagens dos argumentos
boot.t.test <- function(x, y, reps = 1000, mu = 0, alternative = c("two.sided", "less", "greater")){
  # Bootstrap t-test as described in Efron and Tibshirani (1993), (Algorithm 16.2, p224)
  if(is.null(x) | is.null(y)) 
    stop("\nArguments to boot.t.test cannot be NULL!")
  if((length(x) <= 10) | (length(y) <= 10)) 
    stop("\nboot.t.test: Sample lengths cannot be less than 10!")
  
  nx <- length(x)
  ny <- length(y)
  t.obs <- (mean(x) - mean(y) - mu) / sqrt(var(x) / nx + var(y) / ny)
  comb.mean <- mean(c(x, y))
  x.c <- x - mean(x) + comb.mean
  y.c <- y - mean(y) + comb.mean
  t.boot <- rep(NA, times = reps)
  
  bootFunc <- function(){
    bootx <- x.c[sample(1:nx, size = nx, replace = TRUE)]
    booty <- y.c[sample(1:ny, size = ny, replace = TRUE)]
    return((mean(bootx) - mean(booty) - mu) / sqrt(var(bootx) / nx + var(booty) / ny))
  }
  
  t.boot <- replicate(reps, expr = bootFunc())
  
  if(alternative[1] == "two.sided"){
    pval <- length(t.boot[abs(t.boot) >= abs(t.obs)]) / reps
    h1phrase <- "not equal to"
  }
  
  if(alternative[1] == "less"){
    pval <- length(t.boot[t.boot <= t.obs]) / reps
    h1phrase <- "less than"
  }
  
  if(alternative[1] == "greater"){
    pval <- length(t.boot[t.boot >= t.obs]) / reps
    h1phrase <- "greater than"
  }
  
  cat("\nBootstrap Two Sample t-test\n")
  cat(paste("\nt = ", round(t.obs, 3), ", p-value = ", round(pval, 4), "\n", sep=""))
  cat(paste("Alternative hypothesis: true difference in means is ", h1phrase, " ", mu, "\n\n", sep=""))
  
  return(tibble(mu0 = mu,
                statistic = t.obs,
                alternative = alternative[1],
                p.value = pval))
}

# var_test ----------------------------------------------------------------
# Retirado e modificado de rugarch
# https://bitbucket.org/alexiosg/rugarch
# De source/R/rugarch-tests.R
var_test <- function(cover = 0.025, loss, var, conf.level = 0.95) {
  N <- length(loss)
  VaRn <- floor(N * cover)
  if(N != length(var)) stop("\nlength of realized losses not equal to length of VaR!")
  tmp <- LR.cc.test(p = cover, lr_loss = loss, lr_var = var, conf_level = conf.level)
  ans  <- tibble(
    expected.exceed = floor(cover*tmp$TN),
    actual.exceed = tmp$N,
    uc.H0 = "Correct Exceedances",
    uc.LRstat = tmp$stat.uc,
    uc.critical = tmp$crit.val.uc,
    uc.LRp = tmp$p.value.uc,
    uc.Decision = ifelse(uc.LRp<(1-conf.level), "Reject H0", "Fail to Reject H0"),
    cc.H0 = "Correct Exceedances & Independent",
    cc.LRstat = tmp$stat.cc,
    cc.critical = tmp$crit.val.cc,
    cc.LRp = tmp$p.value.cc,
    cc.Decision = ifelse(cc.LRp<(1-conf.level), "Reject H0", "Fail to Reject H0"))
  return(ans)  
}

# LR.cc.test --------------------------------------------------------------
# Retirado e modificado de rugarch
# https://bitbucket.org/alexiosg/rugarch
# De source/R/rugarch-tests.R
# Modificado para fazer bootstrap de stat.ind que pode apresentar por varias vezes
# o valor NaN
########################################################################
# Available in a number of locations/textbooks.
# This code originally presented in a webcast by Guy Yollin Feb 2006 (from Insightful)
# Functions to perform Hypothesis test
# on VaR models based on # of exceedances
# calc LR.uc statistic
LR.cc.test <- function (p, lr_loss, lr_var, conf_level = 0.95) 
{
  result <- .LR.cc(p = p, actual = lr_loss, VaR = lr_var, reps = 1000)
  crit.val.uc <- qchisq(conf_level, df = 1)
  crit.val.cc <- qchisq(conf_level, df = 2)
  p.value.cc <- 1 - pchisq(result$stat.cc, df = 2)
  p.value.uc <- 1 - pchisq(result$stat.uc, df = 1)
  reject <- ifelse(p.value.cc < 1 - conf_level, TRUE, FALSE)
  return(list(stat.cc = result$stat.cc, 
              stat.uc = result$stat.uc, 
              p.value.cc = p.value.cc, 
              p.value.uc = p.value.uc, 
              conf.level = conf_level, 
              reject = reject, 
              N = result$N, 
              TN = result$TN, 
              crit.val.uc = crit.val.uc, 
              crit.val.cc = crit.val.cc))
} # Fim da LR.cc.test

.LR.cc <- function (p, actual, VaR, reps = 1000) 
{
  VaR.ind <- ifelse(actual > VaR, 1, 0)
  N <- sum(VaR.ind)
  TN <- length(VaR.ind)
  # Bootstrap routine
  ind.boot <- rep(0, times = reps)
  bootFunc <- function(){
    bootvar <- VaR.ind[sample(1:TN, size = TN, replace = TRUE)]
    T00 <- sum(c(0, ifelse(bootvar[2:TN] == 0 & bootvar[1:(TN - 1)] == 0, 1, 0)))
    T11 <- sum(c(0, ifelse(bootvar[2:TN] == 1 & bootvar[1:(TN - 1)] == 1, 1, 0)))
    T01 <- sum(c(0, ifelse(bootvar[2:TN] == 1 & bootvar[1:(TN - 1)] == 0, 1, 0)))
    T10 <- sum(c(0, ifelse(bootvar[2:TN] == 0 & bootvar[1:(TN - 1)] == 1, 1, 0)))
    T0 <- T00 + T01
    T1 <- T10 + T11
    pi0 <- T01/T0
    pi1 <- T11/T1
    pe <- (T01 + T11)/(T0 + T1)
    stat.ind <- -2 *( (T00 + T10)*log(1 - pe) + (T01 + T11)*log(pe)) + 2 * (T00*log(1 - pi0)+T01*log(pi0)+T10*log(1 - pi1)+T11*log(pi1))
    return(stat.ind)
  }
  
  ind.boot <- replicate(reps, expr = bootFunc()) # ind.boot may contain NaN
  stat.ind <- mean(ind.boot, na.rm = TRUE)
  # stat.ind = -2 * log((1 - pe)^(T00 + T10) * pe^(T01 + T11)) + 2 * log((1 - pi0)^T00 * pi0^T01 * (1 - pi1)^T10 * pi1^T11)
  stat.uc <- .LR.uc(p = p, TN = TN, N = N)
  stat.cc <- stat.uc + stat.ind
  return(list(stat.cc = stat.cc, stat.uc = stat.uc, N = N, 
              TN = TN))
}

.LR.uc <- function (p, TN, N) 
{
  stat.uc <- -2 *( (TN - N)*log(1 - p)+ N*log(p) ) + 2 * ( (TN - N)*log(1 - N/TN)+N*log(N/TN) )
  return(stat.uc)
}

.Log <- function(x){
  ans <- log(x)
  #if(!is.finite(ans)) ans = sign(ans) * 1e10
  ans
}

# vartest -----------------------------------------------------------------
# Funcao para englobar o teste de Kupiec1995 e de Christoffersen2004
# Chama as funcoes VaRTest e VaRDurTest e agrega apenas as informacoes de interesse
vartest <- function(alpha = 0.025, actual_ret, VaR, conf.level = 0.95) {
  if(length(actual_ret) != length(VaR))
    stop("\nvartest: Length of returns and VaR series are different!")
  
  ans_var <- flatten_df(VaRTest(alpha = alpha, 
                                actual = actual_ret, 
                                VaR = VaR, 
                                conf.level = conf.level))
  ans_dur <- flatten_df(VaRDurTest(alpha = alpha,
                                   actual = actual_ret,
                                   VaR = VaR,
                                   conf.level = conf.level))
  ans <- ans_var %>% 
    select(uc.LRstat, uc.LRp) %>% 
    bind_cols(select(ans_dur, uLL, rLL, LRp))
  
  return(ans)
}



