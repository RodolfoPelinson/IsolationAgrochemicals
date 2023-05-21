gllvm_AICc_tab <- function(x, names, order = FALSE, correction = FALSE){
  aiccs <- rep(NA, length(x))
  df <- rep(NA, length(x))
  nobs <- rep(NA, length(x))
  for(i in 1:length(x)){

    if(isTRUE(correction)){
      aiccs[i] <- gllvm::AICc(x[[i]])
    }else{
      aiccs[i] <- AIC(x[[i]])
    }

    df[i] <- attributes(logLik(x[[i]]))$df
    nobs[i] <- dim(x[[i]]$y)[1]*dim(x[[i]]$y)[2]
  }

  if(order == TRUE){
    aiccs_ord <- order(aiccs)
    aiccs <- aiccs[aiccs_ord]
    df <- df[aiccs_ord]
    nobs <- nobs[aiccs_ord]
    names <- names[aiccs_ord]
  }

  delta_aicc <- aiccs - min(aiccs)

  if(isTRUE(correction)){
    result <- data.frame(model = names,AICc = aiccs, delta_AICc = delta_aicc, df = df, nobs = nobs)
  }else{
    result <- data.frame(model = names,AIC = aiccs, delta_AIC = delta_aicc, df = df, nobs = nobs)
  }


  return(result)
}


