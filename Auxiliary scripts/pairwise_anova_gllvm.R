pairwise_anova_gllvm <- function(model, X_row = NULL, row.eff = FALSE, fact , by = NULL, n.init = 10, seed = 11:20, methop_p = "fdr", ...){

  ####
  #model = model_herb_det_selection_SS3$models$`Isolation * Treatments`
  #X_row = data.frame(ID = ID_SS2_3_4)
  #fact = "treatments"
  #by = "isolation"
  #row.eff = ~ (1|ID)
  #####


  X  <- model$X
  y  <- model$y
  family <- model$family
  method <- model$method
  num.lv <- model$num.lv


  if(is.null(by)){
    fact_pairwise <- X[colnames(X)==fact][,1]
    fact_levels <- levels(as.factor(fact_pairwise))

    combinations <- combn(fact_levels, 2)

    anova_table <- matrix(nrow = ncol(combinations), ncol = 4)

    for(i in 1:ncol(combinations)){
      new_x <- data.frame(as.factor(as.character((X[fact_pairwise == combinations[1,i] | fact_pairwise == combinations[2,i], fact]))))
      colnames(new_x) <- fact

      if(is.null(X_row) == FALSE){
        new_X_row <- data.frame(as.factor(as.character((X_row[fact_pairwise == combinations[1,i] | fact_pairwise == combinations[2,i],]))))
        colnames(new_X_row) <- colnames(X_row)
      }else{
        new_X_row <- NULL
      }

      new_y <- y[fact_pairwise == combinations[1,i] | fact_pairwise == combinations[2,i],]

      mod0 <- gllvm(y = new_y,
                    num.lv = num.lv,
                    X = new_X_row,
                    row.eff = row.eff,
                    family = family,
                    method = method,
                    n.init = n.init,
                    seed = seed, ...)

      if(is.null(X_row) == FALSE){
        new_x <- data.frame(new_x,new_X_row)
      }

      mod1 <- gllvm(formula = as.formula(paste("~", fact)),
                    y = new_y,
                    num.lv = num.lv,
                    X = new_x,
                    row.eff = row.eff,
                    family = family,
                    method = method,
                    n.init = n.init,
                    seed = seed, ...)

      anova <- anova(mod0,mod1)

      anova_table[i,] <- c(paste(combinations[1,i],"-",combinations[2,i]),anova$D[2],anova$Df.diff[2],anova$P.value[2])
    }
    colnames(anova_table) <- c("model","D","Diff.Df","p")
    anova_table <- as.data.frame(anova_table)
    anova_table$p_adj <- p.adjust(as.numeric(anova_table$p), method = methop_p)
    anova_table$D <- round(as.numeric(anova_table$D),3)
    anova_table$p <- round(as.numeric(anova_table$p),3)
    return(anova_table)

  }else{




    by_pairwise <- X[colnames(X)==by][,1]
    by_levels <- levels(as.factor(by_pairwise))

    comparisons <- list()

    for(i in 1:length(by_levels)){


      new_x <- X[ X[which(colnames(X) == by)] == by_levels[i]  ,]
      new_X_row <- data.frame(X_row[ X[which(colnames(X) == by)] == by_levels[i]  ,])
      colnames(new_X_row) <- colnames(X_row)
      new_y <- y[ X[which(colnames(X) == by)] == by_levels[i]  ,]


      fact_pairwise <- new_x[colnames(new_x)==fact][,1]
      fact_levels <- levels(as.factor(fact_pairwise))
      combinations <- combn(fact_levels, 2)

      anova_table <- matrix(nrow = ncol(combinations), ncol = 4)


      for(j in 1:ncol(combinations)){
        new_x2 <- data.frame(as.factor(as.character((new_x[fact_pairwise == combinations[1,j] | fact_pairwise == combinations[2,j], fact]))))
        colnames(new_x2) <- fact

        if(is.null(X_row) == FALSE){
          new_X_row2 <- data.frame(as.factor(as.character((new_X_row[fact_pairwise == combinations[1,j] | fact_pairwise == combinations[2,j],]))))
          colnames(new_X_row2) <- colnames(new_X_row)
        }else{
          new_X_row2 <- NULL
        }

        new_y2 <- new_y[fact_pairwise == combinations[1,j] | fact_pairwise == combinations[2,j],]
        #new_y2 <- new_y2[,colSums(new_y2)>0]

        mod0 <- gllvm(formula = as.formula(paste("~ 1")),
                      y = new_y2,
                      num.lv = num.lv,
                      X = new_X_row2,
                      row.eff = row.eff,
                      family = family,
                      method = method,
                      n.init = n.init,
                      seed = seed, ...)

        if(is.null(X_row) == FALSE){
          new_x2 <- data.frame(new_x2,new_X_row2)
        }

        mod1 <- gllvm(formula = as.formula(paste("~", fact)),
                      y = new_y2,
                      num.lv = num.lv,
                      X = new_x2,
                      row.eff = row.eff,
                      family = family,
                      method = method,
                      n.init = n.init,
                      seed = seed, ...)

        anova <- anova(mod0,mod1)

        anova_table[j,] <- c(paste(combinations[1,j],"-",combinations[2,j]),anova$D[2],anova$Df.diff[2],anova$P.value[2])
      }
      colnames(anova_table) <- c("model","D","Diff.Df","p")
      anova_table <- as.data.frame(anova_table)
      anova_table$p_adj <- round(p.adjust(as.numeric(anova_table$p), method = methop_p),3)
      anova_table$D <- round(as.numeric(anova_table$D),3)
      anova_table$p <- round(as.numeric(anova_table$p),3)
      comparisons[[i]] <- anova_table
    }
    names(comparisons) <- by_levels
    return(comparisons)

  }


}
