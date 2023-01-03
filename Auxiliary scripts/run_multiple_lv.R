run_multiple_lv <- function(num.lv, ...){

  models <- list()

  for(i in 1:length(num.lv)){
    models[[i]] <- gllvm(num.lv = num.lv[i], ...)
  }

  AICc_tab <- gllvm_AICc_tab(models, names = as.character(num.lv), order = FALSE)
  return(AICc_tab)
}
