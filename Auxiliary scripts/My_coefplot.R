My_coefplot <- function (mles, lower, upper, order = NULL,species_labels = NULL, xlab = NULL, cex.axis = 1,y_spa = 0, rect = F,rect_lim = 5, font = 1, main = NULL, cex.main = 1, at.xaxis = NULL, break.axis = NULL, ...)
{

  if(is.null(order) == FALSE){
    mles <- mles[order]
    lower <- lower[order]
    upper <- upper[order]
    if(is.null(species_labels) == FALSE){
      species_labels <- species_labels[order]
    }
  }

  col.seq <- rep("grey", length(mles))
  col.seq[which(lower < 0 & upper < 0)] <- "red"
  col.seq[which(lower > 0 & upper > 0)] <- "blue"

#  col.seq[which(lower < 0 & upper < 0 & mles < -10)] <- "indianred1"
#  col.seq[which(lower > 0 & upper > 0 & mles > 10)] <- "steelblue1"


  lwd.seq <- rep(1, length(mles))
  lwd.seq[which(lower < 0 & upper < 0)] <- 2
  lwd.seq[which(lower > 0 & upper > 0)] <- 2

  #lwd.seq[which(mles < -10 | mles > 10)] <- 1

  cex.seq <- rep(1, length(mles))
  cex.seq[which(lower < 0 & upper < 0)] <- 1.4
  cex.seq[which(lower > 0 & upper > 0)] <- 1.4

  #cex.seq[which(mles < -10 | mles > 10)] <- 1


  At.y <- rev(1:length(mles))
  ylim <- c(min(At.y-y_spa), max(At.y+y_spa))


  if(is.null(break.axis)){
    plot(x = NULL, y = NULL, yaxt = "n", ylim = ylim,
         ylab = "", xlab = xlab, xlim = c(min(lower), max(upper)), xaxt = "n")#,...)
    title(main = main, cex.main = cex.main, adj = 1, line= 0.25)

    if(is.null(at.xaxis)){
      axis(1, gap.axis = -10)
    }else{
      axis(1, gap.axis = -10, at = at.xaxis)
    }
  }else{

    if(is.null(at.xaxis)){
      stop("Must inform at.xaxis")
    }else{
      new_mles <- mles
      new_upper <- upper
      new_lower <- lower


      for(i in 1:length(mles)){
        if(mles[i] < break.axis[1]){
          mles[i] <- mles[i] - (break.axis[2] - break.axis[1])
          upper[i] <- upper[i] - (break.axis[2] - break.axis[1])
          lower[i] <- lower[i] - (break.axis[2] - break.axis[1])
        }
      }


      plot(x = NULL, y = NULL, yaxt = "n", ylim = ylim,
           ylab = "", xlab = xlab, xlim = c(min(lower), max(upper)), xaxt = "n")#,...)
      title(main = main, cex.main = cex.main, adj = 1, line= 0.25)

      #at.xaxis <- seq(from = 200, to = -200, by = -1)

      new_at.xaxis <- at.xaxis
      new_at.xaxis_up_zero <- new_at.xaxis[new_at.xaxis >= 0]
      new_at.xaxis_below_zero <- new_at.xaxis[new_at.xaxis < 0]


      if(break.axis[1] > 0){
        new_at.xaxis_up_zero <-  new_at.xaxis_up_zero[abs(new_at.xaxis_up_zero) < abs(break.axis[1]) | abs(new_at.xaxis_up_zero) > abs(break.axis[2])]
        new_at.xaxis_below_zero_labels <- as.character(new_at.xaxis_below_zero)
        new_at.xaxis_up_zero_labels <- as.character(new_at.xaxis_up_zero)
        new_at.xaxis_up_zero[abs(new_at.xaxis_up_zero) > abs(break.axis[2])] <- new_at.xaxis_up_zero[abs(new_at.xaxis_up_zero) > abs(break.axis[2])] - (break.axis[2] - break.axis[1])
      }

      if(break.axis[1] < 0){
        new_at.xaxis_below_zero <-  new_at.xaxis_below_zero[abs(new_at.xaxis_below_zero) < abs(break.axis[1]) | abs(new_at.xaxis_below_zero) > abs(break.axis[2])]
        new_at.xaxis_below_zero_labels <- as.character(new_at.xaxis_below_zero)
        new_at.xaxis_up_zero_labels <- as.character(new_at.xaxis_up_zero)
        new_at.xaxis_below_zero[abs(new_at.xaxis_below_zero) > abs(break.axis[2])] <- new_at.xaxis_below_zero[abs(new_at.xaxis_below_zero) > abs(break.axis[2])] - (break.axis[2] - break.axis[1])
      }

      new_at.xaxis <- c(new_at.xaxis_up_zero, new_at.xaxis_below_zero)
      new_at.xaxis_labels <- c(new_at.xaxis_up_zero_labels, new_at.xaxis_below_zero_labels)



      axis(1, gap.axis = -10, at = new_at.xaxis, labels = new_at.xaxis_labels)
      plotrix::axis.break(1, break.axis[1], breakcol="black", style="slash")

    }


  }



  if(isTRUE(rect)){
  rect(xleft = min(lower)-1, ybottom = 0, xright =  max(upper)+1, ytop = max(At.y)-rect_lim+0.5, density = NULL, border = "transparent", col = rgb(col2rgb("cornflowerblue", alpha = FALSE)[1],
                                                                                                          col2rgb("cornflowerblue", alpha = FALSE)[2],
                                                                                                          col2rgb("cornflowerblue", alpha = FALSE)[3],
                                                                                                          alpha = 40, maxColorValue = 255))
  rect(xleft = min(lower)-1, ybottom = max(At.y)-rect_lim+0.5, xright =  max(upper)+1, ytop = max(At.y)+1, density = NULL, border = "transparent", col = rgb(col2rgb("coral3", alpha = FALSE)[1],
                                                                                                                                col2rgb("coral3", alpha = FALSE)[2],
                                                                                                                                col2rgb("coral3", alpha = FALSE)[3],
                                                                                                                                alpha = 40, maxColorValue = 255))
  }

  points(x = mles, y = At.y, col = col.seq, pch = 16, cex = cex.seq)
  #segments(x0 = lower,
  #         x1 = upper,
  #         y1 = At.y, y0 = At.y, col = col.seq, lwd =lwd.seq)


  arrows(y1 = At.y, y0 = At.y, x1 = upper, x0 = lower,
         code = 3, angle = 90, length = 0.025,col = col.seq, lwd =lwd.seq)


  abline(v = 0, lty = 3)


  if(length(font) > 1){
    for(i in 1:length(species_labels)){
      axis(2, at = At.y[i], labels = species_labels[i], las = 1, cex.axis = cex.axis, font = font[i])
    }
  }else{
    axis(2, at = At.y, labels = species_labels, las = 1, cex.axis = cex.axis, font = font)
  }



}



