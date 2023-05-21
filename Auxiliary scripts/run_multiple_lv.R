run_multiple_lv <- function(num.lv, correction = TRUE,...){

  #begin <- Sys.time()

  models <- list()

  for(i in 1:length(num.lv)){
    models[[i]] <- gllvm(num.lv = num.lv[i], ...)
  }

  AICc_tab <- gllvm_AICc_tab(models, names = as.character(num.lv), order = FALSE, correction = correction)

  #message(paste("Elapsed time =",round(Sys.time()-begin,2)))

  return(list(AICc_tab = AICc_tab, models = models))
}





#run_multiple_lv2 <- function(num.lv, correction = TRUE, ...){
#  require(foreach)
#  require(doParallel)
  #require(gllvm)

#  cl <- makeCluster(detectCores() - 2)

#  registerDoParallel(cl)

#  begin <- Sys.time()

#  models <- list()

#  models <- foreach(i = 1:length(num.lv))%dopar%{
#    gllvm(num.lv = num.lv[i], ...)
#  }

#  AICc_tab <- gllvm_AICc_tab(models, names = as.character(num.lv), order = FALSE, correction = correction)

#  message(paste("Elapsed time =",round(Sys.time()-begin,2)))
#  stopCluster(cl)

#  return(list(AICc_tab = AICc_tab, models = models))
#}
