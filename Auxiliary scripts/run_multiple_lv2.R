run_multiple_lv2 <- function(lv_vector, correction = TRUE, family = "negative.binomial",
                             method = "VA", n.init = 10, seed = 1:10,...){

  require(parallel)

  n_cores <- detectCores() - 2 #Define o número de cores a serem usados
  message(paste("Using ", n_cores, " cores")) #Mensagem do número de cores usados
  #registerDoParallel(cores=n_cores)
  cl <- makeCluster(n_cores) #Fazendo os clusters
  #Incluindo os pacotes usados e objetos que será usado em cada cluster
  clusterEvalQ(cl, c("lv_vector", "gllvm2"))
  clusterEvalQ(cl, {library(gllvm)})


  models <- parLapply(cl, lv_vector, gllvm, ...)
  stopCluster(cl)


  AICc_tab <- gllvm_AICc_tab(models, names = as.character(num.lv), order = FALSE, correction = correction)

  return(list(AICc_tab = AICc_tab, models = models))
}



gllvm(formula = ~ treatments * isolation,
       num.lv = 0,
       y = com_herb_det_SS1, X = SS1_predictors,
       family = "negative.binomial",
       method = "VA",
       n.init = 10, seed = 11:20, starting.val = 0)



dobro <- function(x,y){x*2/y}


list <- c(1,2,3,4)
lapply(list, dobro, y = 3)
