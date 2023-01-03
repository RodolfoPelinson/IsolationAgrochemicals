run_multiple_gllvm <- function(formulas, num.lv = 0, names, no_effect = FALSE, X = NULL, X_row = NULL, ...){
  models <- list()

  if(isTRUE(no_effect)){
    first <- 2
    last <- length(formulas)+1
    names <- c("No Effect", names)
    rev_i <- 1
  }else{
    first <- 1
    last <- length(formulas)
    rev_i <- 0
  }

  for(i in first:last){
    models[[i]] <- gllvm(formula = formulas[[i-rev_i]], num.lv = num.lv, X, ...)
    #names(models)[[i]] <- formulas[[i-rev_i]]
  }
  if(isTRUE(no_effect)){
    models[[1]] <- gllvm(num.lv = num.lv, X = X_row, formula = NULL, ...)
  }

  AICc_tab <- gllvm_AICc_tab(models, names = names, order = FALSE)
  return(AICc_tab)
}
