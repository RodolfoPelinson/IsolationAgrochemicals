library(colorspace)
col_control_30 <- "#56B4E9"
col_control_120 <- darken(col_control_30, amount = 0.25, space = "HCL")
col_control_480 <- darken(col_control_30, amount = 0.5, space = "HCL")

col_pasture_30 <- "#F0E442"
col_pasture_120 <- darken(col_pasture_30, amount = 0.25, space = "HCL")
col_pasture_480 <- darken(col_pasture_30, amount = 0.5, space = "HCL")

col_sugarcane_30 <- "#D55E00"
col_sugarcane_120 <- darken(col_sugarcane_30, amount = 0.25, space = "HCL")
col_sugarcane_480 <- darken(col_sugarcane_30, amount = 0.5, space = "HCL")


#SS1

plot_com_SS1_herb_treat_30 <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS1_predictors$isolation == "30",],
              groups = SS1_predictors$treatments[SS1_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)





  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",])

  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)

  lines(x = c(centroid_control_30[1], centroid_pasture_30[1]),
        y = c(centroid_control_30[2], centroid_pasture_30[2]))

  lines(x = c(centroid_control_30[1], centroid_sugarcane_30[1]),
        y = c(centroid_control_30[2], centroid_sugarcane_30[2]))

  lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
        y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))

  box()

  title(main = "30 m", adj = 1, line = 0.25)
}

plot_com_SS1_herb_treat_120 <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS1_predictors$isolation == "120",],
              groups = SS1_predictors$treatments[SS1_predictors$isolation == "120"],
              draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_120[1], centroid_pasture_120[1]),
        y = c(centroid_control_120[2], centroid_pasture_120[2]))

  lines(x = c(centroid_control_120[1], centroid_sugarcane_120[1]),
        y = c(centroid_control_120[2], centroid_sugarcane_120[2]))

  lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
        y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))


  box()

  title(main = "120 m", adj = 1, line = 0.25)
}

plot_com_SS1_herb_treat_480 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS1_predictors$isolation == "480",],
              groups = SS1_predictors$treatments[SS1_predictors$isolation == "480"],
              draw = "polygon", border = FALSE , lty = 1, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)
  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",])



  #lines(x = c(centroid_sugarcane_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_480[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_480[1], centroid_pasture_480[1]),
        y = c(centroid_control_480[2], centroid_pasture_480[2]))

  lines(x = c(centroid_control_480[1], centroid_sugarcane_480[1]),
        y = c(centroid_control_480[2], centroid_sugarcane_480[2]))

  lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
        y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))


  box()

  title(main = "480 m", adj = 1, line = 0.25)}

plot_com_SS1_herb_iso_control <- function(){

  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS1_predictors$treatments == "control",],
              groups = SS1_predictors$isolation[SS1_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_480[1]),
        y = c(centroid_sugarcane_120[2], centroid_sugarcane_480[2]))



  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  box()

  title(main = "Control", adj = 1, line = 0.25)

  }

plot_com_SS1_herb_iso_pasture <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS1_predictors$treatments == "pasture",],
              groups = SS1_predictors$isolation[SS1_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_pasture_120[1], centroid_pasture_480[1]),
        y = c(centroid_pasture_120[2], centroid_pasture_480[2]))


  box()

  title(main = "Pasture", adj = 1, line = 0.25)
}

plot_com_SS1_herb_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane",],
              groups = SS1_predictors$isolation[SS1_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_480[1]),
        y = c(centroid_sugarcane_120[2], centroid_sugarcane_480[2]))


  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}


#SS2
plot_com_SS2_herb_treat_30 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS2_predictors$isolation == "30",],
              groups = SS2_predictors$treatments[SS2_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  #ordiellipse(scaled_lvs$sites,
  #            groups = SS2_predictors$treatments,
  #            draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",])

  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control",])
  #centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture",])
  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane",])


  lines(x = c(centroid_control_30[1], centroid_pasture_30[1]),
        y = c(centroid_control_30[2], centroid_pasture_30[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)

  box()

  title(main = "30 m", adj = 1, line = 0.25)
}

plot_com_SS2_herb_treat_120 <- function(){par(mar = c(4,4,1.25,.1))
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS2_predictors$isolation == "120",],
              groups = SS2_predictors$treatments[SS2_predictors$isolation == "120"],
              draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)



  lines(x = c(centroid_control_120[1], centroid_pasture_120[1]),
        y = c(centroid_control_120[2], centroid_pasture_120[2]))


  box()

  title(main = "120 m", adj = 1, line = 0.25)}

plot_com_SS2_herb_treat_480 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS2_predictors$isolation == "480",],
              groups = SS2_predictors$treatments[SS2_predictors$isolation == "480"],
              draw = "polygon", border = FALSE , lty = 1, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)
  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",])



  #lines(x = c(centroid_sugarcane_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_480[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_480[1], centroid_pasture_480[1]),
        y = c(centroid_control_480[2], centroid_pasture_480[2]))


  box()

  title(main = "480 m", adj = 1, line = 0.25)
}

plot_com_SS2_herb_iso_control <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS2_predictors$treatments == "control",],
              groups = SS2_predictors$isolation[SS2_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_control_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",])
  centroid_control_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",])
  centroid_control_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }



  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  box()

  title(main = "Control", adj = 1, line = 0.25)
}

plot_com_SS2_herb_iso_pasture <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS2_predictors$treatments == "pasture",],
              groups = SS2_predictors$isolation[SS2_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  #lines(x = c(centroid_pasture_30[1], centroid_pasture_120[1]),
  #      y = c(centroid_pasture_30[2], centroid_pasture_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_pasture_480[1]),
  #      y = c(centroid_pasture_120[2], centroid_pasture_480[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_pasture_480[1]),
  #      y = c(centroid_pasture_30[2], centroid_pasture_480[2]))

  box()

  title(main = "Pasture", adj = 1, line = 0.25)
}

plot_com_SS2_herb_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane",],
              groups = SS2_predictors$isolation[SS2_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)



  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}


#SS3
plot_com_SS3_herb_treat_30 <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS3_predictors$isolation == "30",],
              groups = SS3_predictors$treatments[SS3_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  #ordiellipse(scaled_lvs$sites,
  #            groups = SS3_predictors$treatments,
  #            draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",])

  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control",])
  #centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture",])
  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane",])


  #lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_sugarcane_30[2], centroid_sugarcane_30[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)

  lines(x = c(centroid_control_30[1], centroid_pasture_30[1]),
        y = c(centroid_control_30[2], centroid_pasture_30[2]))

  box()

  title(main = "30 m", adj = 1, line = 0.25)
}

plot_com_SS3_herb_treat_120 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS3_predictors$isolation == "120",],
              groups = SS3_predictors$treatments[SS3_predictors$isolation == "120"],
              draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)



  lines(x = c(centroid_control_120[1], centroid_pasture_120[1]),
        y = c(centroid_control_120[2], centroid_pasture_120[2]))


  box()

  title(main = "120 m", adj = 1, line = 0.25)
}

plot_com_SS3_herb_treat_480 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS3_predictors$isolation == "480",],
              groups = SS3_predictors$treatments[SS3_predictors$isolation == "480"],
              draw = "polygon", border = FALSE , lty = 1, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)
  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",])



  #lines(x = c(centroid_sugarcane_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_480[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_480[1], centroid_pasture_480[1]),
        y = c(centroid_control_480[2], centroid_pasture_480[2]))

  box()

  title(main = "480 m", adj = 1, line = 0.25)
}

plot_com_SS3_herb_iso_control <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS3_predictors$treatments == "control",],
              groups = SS3_predictors$isolation[SS3_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_control_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",])
  centroid_control_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",])
  centroid_control_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }



  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  lines(x = c(centroid_control_30[1], centroid_control_120[1]),
        y = c(centroid_control_30[2], centroid_control_120[2]))

  box()

  title(main = "Control", adj = 1, line = 0.25)
}

plot_com_SS3_herb_iso_pasture <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS3_predictors$treatments == "pasture",],
              groups = SS3_predictors$isolation[SS3_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  #lines(x = c(centroid_pasture_30[1], centroid_pasture_120[1]),
  #      y = c(centroid_pasture_30[2], centroid_pasture_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_pasture_480[1]),
  #      y = c(centroid_pasture_120[2], centroid_pasture_480[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_pasture_480[1]),
  #      y = c(centroid_pasture_30[2], centroid_pasture_480[2]))

  box()

  title(main = "Pasture", adj = 1, line = 0.25)
}

plot_com_SS3_herb_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane",],
              groups = SS3_predictors$isolation[SS3_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_120[1]),
        y = c(centroid_sugarcane_30[2], centroid_sugarcane_120[2]))

  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}


#SS4
plot_com_SS4_herb_treat_30 <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "30",],
              groups = SS4_predictors$treatments[SS4_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  ordiellipse(scaled_lvs$sites[SS1_predictors$isolation == "30",],
              groups = SS1_predictors$treatments[SS1_predictors$isolation == "30"],
              draw = "polygon", border = c(col_control_120, col_pasture_120, col_sugarcane_120),  col = NULL, kind = "sd", alpha = 50, lwd = 2)

  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",])

  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control",])
  #centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture",])
  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane",])


  #lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_sugarcane_30[2], centroid_sugarcane_30[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  #points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  #points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  #points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)

  box()

  title(main = "30 m", adj = 1, line = 0.25)
}

plot_com_SS4_herb_treat_120 <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  #ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "120",],
  #            groups = SS4_predictors$treatments[SS4_predictors$isolation == "120"],
  #            draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)


  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "120",],
              groups = SS4_predictors$treatments[SS4_predictors$isolation == "120"],
              draw = "polygon", border = c(col_control_120, col_pasture_120, col_sugarcane_120) , lty = 2, col = NULL, kind = "sd", alpha = 50)

  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "120",],
              groups = SS4_predictors$treatments_contpast_sug[SS4_predictors$isolation == "120"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_sugarcane_120), kind = "sd", alpha = 50)

  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "120",],
              groups = SS4_predictors$treatments_contpast_sug[SS4_predictors$isolation == "120"],
              draw = "polygon", border = c(col_pasture_120, col_sugarcane_120),  col = NULL, kind = "sd", alpha = 50, lwd = 2)


  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  #points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  #points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  #points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)



  #lines(x = c(centroid_control_120[1], centroid_pasture_120[1]),
  #      y = c(centroid_control_120[2], centroid_pasture_120[2]))


  box()

  title(main = "120 m", adj = 1, line = 0.25)
}

plot_com_SS4_herb_treat_480 <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "480",],
              groups = SS4_predictors$treatments[SS4_predictors$isolation == "480"],
              draw = "polygon", border = c(col_control_120, col_pasture_120, col_sugarcane_120) , lty = 2, col = NULL, kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",])



  #lines(x = c(centroid_sugarcane_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_480[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  #points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  #points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  #points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  #lines(x = c(centroid_control_480[1], centroid_pasture_480[1]),
  #      y = c(centroid_control_480[2], centroid_pasture_480[2]))

  #lines(x = c(centroid_control_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_control_480[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))

  box()

  title(main = "480 m", adj = 1, line = 0.25)
}

plot_com_SS4_herb_iso_control <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS4_predictors$treatments == "control",],
              groups = SS4_predictors$isolation[SS4_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_control_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",])
  centroid_control_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",])
  centroid_control_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }



  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  lines(x = c(centroid_control_480[1], centroid_control_120[1]),
        y = c(centroid_control_480[2], centroid_control_120[2]))

  box()

  title(main = "Control", adj = 1, line = 0.25)}

plot_com_SS4_herb_iso_pasture <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS4_predictors$treatments == "pasture",],
              groups = SS4_predictors$isolation[SS4_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  #lines(x = c(centroid_pasture_30[1], centroid_pasture_120[1]),
  #      y = c(centroid_pasture_30[2], centroid_pasture_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_pasture_480[1]),
  #      y = c(centroid_pasture_120[2], centroid_pasture_480[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_pasture_480[1]),
  #      y = c(centroid_pasture_30[2], centroid_pasture_480[2]))

  box()

  title(main = "Pasture", adj = 1, line = 0.25)
}

plot_com_SS4_herb_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane",],
              groups = SS4_predictors$isolation[SS4_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_sugarcane_480[1], centroid_sugarcane_120[1]),
        y = c(centroid_sugarcane_480[2], centroid_sugarcane_120[2]))

  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}



#Predators

#SS1
plot_com_SS1_pred_treat_30 <- function(){
  plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS1_predictors$isolation == "30",],
              groups = SS1_predictors$treatments[SS1_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  #ordiellipse(scaled_lvs$sites,
  #            groups = SS1_predictors$treatments,
  #            draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",])

  #centroid_control_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control",])
  #centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture",])
  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane",])


  #lines(x = c(centroid_control_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_control_30[2], centroid_sugarcane_30[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)

  lines(x = c(centroid_control_30[1], centroid_pasture_30[1]),
        y = c(centroid_control_30[2], centroid_pasture_30[2]))

  lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
        y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))

  lines(x = c(centroid_control_30[1], centroid_sugarcane_30[1]),
        y = c(centroid_control_30[2], centroid_sugarcane_30[2]))



  title(main = "30 m", adj = 1, line = 0.25)

  box()

}

plot_com_SS1_pred_treat_120 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS1_predictors$isolation == "120",],
              groups = SS1_predictors$treatments[SS1_predictors$isolation == "120"],
              draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_120[1], centroid_pasture_120[1]),
        y = c(centroid_control_120[2], centroid_pasture_120[2]))

  lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
        y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))

  lines(x = c(centroid_control_120[1], centroid_sugarcane_120[1]),
        y = c(centroid_control_120[2], centroid_sugarcane_120[2]))




  box()

  title(main = "120 m", adj = 1, line = 0.25)
}

plot_com_SS1_pred_treat_480 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS1_predictors$isolation == "480",],
              groups = SS1_predictors$treatments[SS1_predictors$isolation == "480"],
              draw = "polygon", border = FALSE , lty = 1, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)
  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",])



  #lines(x = c(centroid_sugarcane_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_480[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_480[1], centroid_pasture_480[1]),
        y = c(centroid_control_480[2], centroid_pasture_480[2]))

  lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
        y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))


  lines(x = c(centroid_control_480[1], centroid_sugarcane_480[1]),
        y = c(centroid_control_480[2], centroid_sugarcane_480[2]))

  box()

  title(main = "480 m", adj = 1, line = 0.25)}

plot_com_SS1_pred_iso_control <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS1_predictors$treatments == "control",],
              groups = SS1_predictors$isolation[SS1_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_control_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "30",])
  centroid_control_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "120",])
  centroid_control_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "control" & SS1_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_control_30[1], centroid_control_120[1]),
        y = c(centroid_control_30[2], centroid_control_120[2]))

  lines(x = c(centroid_control_120[1], centroid_control_480[1]),
        y = c(centroid_control_120[2], centroid_control_480[2]))

  lines(x = c(centroid_control_30[1], centroid_control_480[1]),
        y = c(centroid_control_30[2], centroid_control_480[2]))







  box()

  title(main = "Control", adj = 1, line = 0.25)
}

plot_com_SS1_pred_iso_pasture <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS1_predictors$treatments == "pasture",],
              groups = SS1_predictors$isolation[SS1_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "pasture" & SS1_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_pasture_30[1], centroid_pasture_120[1]),
        y = c(centroid_pasture_30[2], centroid_pasture_120[2]))

  lines(x = c(centroid_pasture_120[1], centroid_pasture_480[1]),
        y = c(centroid_pasture_120[2], centroid_pasture_480[2]))

  lines(x = c(centroid_pasture_30[1], centroid_pasture_480[1]),
        y = c(centroid_pasture_30[2], centroid_pasture_480[2]))




  box()

  title(main = "Pasture", adj = 1, line = 0.25)
}

plot_com_SS1_pred_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane",],
              groups = SS1_predictors$isolation[SS1_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS1_predictors$treatments == "sugar_cane" & SS1_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS1[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_120[1]),
        y = c(centroid_sugarcane_30[2], centroid_sugarcane_120[2]))

  lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_480[1]),
        y = c(centroid_sugarcane_120[2], centroid_sugarcane_480[2]))

  lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_480[1]),
        y = c(centroid_sugarcane_30[2], centroid_sugarcane_480[2]))

  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}



#SS2
plot_com_SS2_pred_treat_30 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS2_predictors$isolation == "30",],
              groups = SS2_predictors$treatments[SS2_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  #ordiellipse(scaled_lvs$sites,
  #            groups = SS2_predictors$treatments,
  #            draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",])

  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control",])
  #centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture",])
  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane",])


  #lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_sugarcane_30[2], centroid_sugarcane_30[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)



  lines(x = c(centroid_control_30[1], centroid_pasture_30[1]),
        y = c(centroid_control_30[2], centroid_pasture_30[2]))



  title(main = "30 m", adj = 1, line = 0.25)

  box()
}

plot_com_SS2_pred_treat_120 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS2_predictors$isolation == "120",],
              groups = SS2_predictors$treatments[SS2_predictors$isolation == "120"],
              draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_120[1], centroid_pasture_120[1]),
        y = c(centroid_control_120[2], centroid_pasture_120[2]))


  box()

  title(main = "120 m", adj = 1, line = 0.25)

}

plot_com_SS2_pred_treat_480 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS2_predictors$isolation == "480",],
              groups = SS2_predictors$treatments[SS2_predictors$isolation == "480"],
              draw = "polygon", border = FALSE , lty = 1, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)
  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",])




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  lines(x = c(centroid_control_480[1], centroid_pasture_480[1]),
        y = c(centroid_control_480[2], centroid_pasture_480[2]))



  box()

  title(main = "480 m", adj = 1, line = 0.25)}

plot_com_SS2_pred_iso_control <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS2_predictors$treatments == "control",],
              groups = SS2_predictors$isolation[SS2_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_control_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "30",])
  centroid_control_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "120",])
  centroid_control_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "control" & SS2_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_control_30[1], centroid_control_120[1]),
        y = c(centroid_control_30[2], centroid_control_120[2]))

  lines(x = c(centroid_control_120[1], centroid_control_480[1]),
        y = c(centroid_control_120[2], centroid_control_480[2]))

  lines(x = c(centroid_control_30[1], centroid_control_480[1]),
        y = c(centroid_control_30[2], centroid_control_480[2]))







  box()

  title(main = "Control", adj = 1, line = 0.25)
}

plot_com_SS2_pred_iso_pasture <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS2_predictors$treatments == "pasture",],
              groups = SS2_predictors$isolation[SS2_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "pasture" & SS2_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_pasture_30[1], centroid_pasture_120[1]),
        y = c(centroid_pasture_30[2], centroid_pasture_120[2]))

  lines(x = c(centroid_pasture_120[1], centroid_pasture_480[1]),
        y = c(centroid_pasture_120[2], centroid_pasture_480[2]))

  lines(x = c(centroid_pasture_30[1], centroid_pasture_480[1]),
        y = c(centroid_pasture_30[2], centroid_pasture_480[2]))




  box()

  title(main = "Pasture", adj = 1, line = 0.25)

}

plot_com_SS2_pred_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane",],
              groups = SS2_predictors$isolation[SS2_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS2_predictors$treatments == "sugar_cane" & SS2_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS2[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_120[1]),
        y = c(centroid_sugarcane_30[2], centroid_sugarcane_120[2]))

  lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_480[1]),
        y = c(centroid_sugarcane_120[2], centroid_sugarcane_480[2]))

  lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_480[1]),
        y = c(centroid_sugarcane_30[2], centroid_sugarcane_480[2]))

  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}




#SS3
plot_com_SS3_pred_treat_30 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS3_predictors$isolation == "30",],
              groups = SS3_predictors$treatments[SS3_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  #ordiellipse(scaled_lvs$sites,
  #            groups = SS3_predictors$treatments,
  #            draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",])

  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control",])
  #centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture",])
  #centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane",])


  #lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_sugarcane_30[2], centroid_sugarcane_30[2]))

  #lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
  #      y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)





  title(main = "30 m", adj = 1, line = 0.25)

  box()}

plot_com_SS3_pred_treat_120 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS3_predictors$isolation == "120",],
              groups = SS3_predictors$treatments[SS3_predictors$isolation == "120"],
              draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  box()

  title(main = "120 m", adj = 1, line = 0.25)

}

plot_com_SS3_pred_treat_480 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS3_predictors$isolation == "480",],
              groups = SS3_predictors$treatments[SS3_predictors$isolation == "480"],
              draw = "polygon", border = FALSE , lty = 1, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)
  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",])



  #lines(x = c(centroid_sugarcane_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_480[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)


  box()

  title(main = "480 m", adj = 1, line = 0.25)
}

plot_com_SS3_pred_iso_control <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS3_predictors$treatments == "control",],
              groups = SS3_predictors$isolation[SS3_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_control_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "30",])
  centroid_control_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "120",])
  centroid_control_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "control" & SS3_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_control_30[1], centroid_control_120[1]),
        y = c(centroid_control_30[2], centroid_control_120[2]))




  box()

  title(main = "Control", adj = 1, line = 0.25)}

plot_com_SS3_pred_iso_pasture <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS3_predictors$treatments == "pasture",],
              groups = SS3_predictors$isolation[SS3_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "pasture" & SS3_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  lines(x = c(centroid_pasture_30[1], centroid_pasture_120[1]),
        y = c(centroid_pasture_30[2], centroid_pasture_120[2]))



  box()

  title(main = "Pasture", adj = 1, line = 0.25)
}

plot_com_SS3_pred_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane",],
              groups = SS3_predictors$isolation[SS3_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS3_predictors$treatments == "sugar_cane" & SS3_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS3[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_120[1]),
        y = c(centroid_sugarcane_30[2], centroid_sugarcane_120[2]))

  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}





#SS4
plot_com_SS4_pred_treat_30 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "30",],
              groups = SS4_predictors$treatments[SS4_predictors$isolation == "30"],
              draw = "polygon", border = FALSE,  col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",], bg = col_sugarcane_30, pch = 21)


  centroid_control_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",])
  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",])
  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)



  lines(x = c(centroid_control_30[1], centroid_pasture_30[1]),
        y = c(centroid_control_30[2], centroid_pasture_30[2]))

  lines(x = c(centroid_pasture_30[1], centroid_sugarcane_30[1]),
        y = c(centroid_pasture_30[2], centroid_sugarcane_30[2]))




  title(main = "30 m", adj = 1, line = 0.25)

  box()
}

plot_com_SS4_pred_treat_120 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "120",],
              groups = SS4_predictors$treatments[SS4_predictors$isolation == "120"],
              draw = "polygon", border = FALSE , lty = 2, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)

  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",], bg = col_sugarcane_30, pch = 21)

  centroid_control_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",])

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_120[2]))

  #lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
  #      y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))




  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)



  lines(x = c(centroid_control_120[1], centroid_pasture_120[1]),
        y = c(centroid_control_120[2], centroid_pasture_120[2]))

  lines(x = c(centroid_pasture_120[1], centroid_sugarcane_120[1]),
        y = c(centroid_pasture_120[2], centroid_sugarcane_120[2]))



  box()

  title(main = "120 m", adj = 1, line = 0.25)
}

plot_com_SS4_pred_treat_480 <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)

  ordiellipse(scaled_lvs$sites[SS4_predictors$isolation == "480",],
              groups = SS4_predictors$treatments[SS4_predictors$isolation == "480"],
              draw = "polygon", border = FALSE , lty = 1, col = c(col_control_120, col_pasture_120, col_sugarcane_120), kind = "sd", alpha = 50)
  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",], bg = col_control_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",], bg = col_pasture_30, pch = 21)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",], bg = col_sugarcane_30, pch = 21)



  centroid_control_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_30, pch = 23, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_pasture_30, pch = 23, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_sugarcane_30, pch = 23, cex = 1.5)



  lines(x = c(centroid_control_480[1], centroid_pasture_480[1]),
        y = c(centroid_control_480[2], centroid_pasture_480[2]))

  lines(x = c(centroid_pasture_480[1], centroid_sugarcane_480[1]),
        y = c(centroid_pasture_480[2], centroid_sugarcane_480[2]))

  box()

  title(main = "480 m", adj = 1, line = 0.25)}

plot_com_SS4_pred_iso_control <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS4_predictors$treatments == "control",],
              groups = SS4_predictors$isolation[SS4_predictors$treatments == "control"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_control_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "30",])
  centroid_control_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "120",])
  centroid_control_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "control" & SS4_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_control_30[1], centroid_control_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_control_120[1], centroid_control_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_control_480[1], centroid_control_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  box()

  title(main = "Control", adj = 1, line = 0.25)
}

plot_com_SS4_pred_iso_pasture <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS4_predictors$treatments == "pasture",],
              groups = SS4_predictors$isolation[SS4_predictors$treatments == "pasture"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_pasture_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "30",])
  centroid_pasture_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "120",])
  centroid_pasture_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "pasture" & SS4_predictors$isolation == "480",])



  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }

  points(centroid_pasture_30[1], centroid_pasture_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_pasture_120[1], centroid_pasture_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_pasture_480[1], centroid_pasture_480[2], bg = col_control_480, pch = 24, cex = 1.5)

  #lines(x = c(centroid_pasture_30[1], centroid_pasture_120[1]),
  #      y = c(centroid_pasture_30[2], centroid_pasture_120[2]))

  #lines(x = c(centroid_sugarcane_120[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_120[2], centroid_sugarcane_480[2]))

  #lines(x = c(centroid_sugarcane_30[1], centroid_sugarcane_480[1]),
  #      y = c(centroid_sugarcane_30[2], centroid_sugarcane_480[2]))

  box()

  title(main = "Pasture", adj = 1, line = 0.25)}

plot_com_SS4_pred_iso_sugarcane <- function(){plot(NA,xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "LV2", xlab = "LV1", axes = F)
  axis(1 , gap.axis = -10)
  axis(2 , gap.axis = -10)
  abline(h=0,v=0, lty =2)
  ordiellipse(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane",],
              groups = SS4_predictors$isolation[SS4_predictors$treatments == "sugar_cane"],
              draw = "polygon", border = "transparent", lty = c(3,2,1), col = c(col_control_30, col_control_120, col_control_480), kind = "sd", alpha = 50)


  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",], col = col_control_30, pch = 15)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",], col = col_control_120, pch = 16)
  points(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",], col = col_control_480, pch = 17)

  centroid_sugarcane_30 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "30",])
  centroid_sugarcane_120 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "120",])
  centroid_sugarcane_480 <- colMeans(scaled_lvs$sites[SS4_predictors$treatments == "sugar_cane" & SS4_predictors$isolation == "480",])


  for(i in 1:nrow(scaled_lvs$new_species)){
    points(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], cex = 4, pch = 21, bg = col_SS4[i])
    text(x = scaled_lvs$new_species[i,1], y = scaled_lvs$new_species[i,2], labels = substr(rownames(scaled_lvs$species)[i], 1, 3), cex = 0.9)
  }


  points(centroid_sugarcane_30[1], centroid_sugarcane_30[2], bg = col_control_30, pch = 22, cex = 1.5)
  points(centroid_sugarcane_120[1], centroid_sugarcane_120[2], bg = col_control_120, pch = 21, cex = 1.5)
  points(centroid_sugarcane_480[1], centroid_sugarcane_480[2], bg = col_control_480, pch = 24, cex = 1.5)


  box()

  title(main = "Sugarcane", adj = 1, line = 0.25)
}


