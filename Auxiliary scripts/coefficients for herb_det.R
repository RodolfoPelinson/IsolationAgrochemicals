com_herb_det <- com[,Trait$trophic == "consumer"]
com_herb_det_SS1 <- com_SS1[,Trait_SS1$trophic == "consumer"]
com_herb_det_SS2 <- com_SS2[,Trait_SS2$trophic == "consumer"]
com_herb_det_SS3 <- com_SS3[,Trait_SS3$trophic == "consumer"]
com_herb_det_SS4 <- com_SS4[,Trait_SS4$trophic == "consumer"]

z<-qnorm(0.975, mean = 0, sd = 1)






############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################ TREATMENTS
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################






############################################################################################################################
###################################################### SS2 #################################################################
############################################################################################################################
fit_herb_det_SS2_treatments<- gllvm(com_herb_det_SS2, X = data.frame(ID = ID_SS2_3_4,
                                                                     treatments = relevel(treatments_SS2_3_4_contpast_sug, ref = "control_pasture"),
                                                                     isolation = isolation_SS2_3_4),
                                    formula = ~ treatments + isolation,
                                    family = "negative.binomial",
                                    method = "VA",
                                    row.eff = ~ (1|ID),
                                    n.init = 20, num.lv = 0, seed = 1:20)



#effect_SS2 of TREATMENTS
effect_SS2_control_sugar_cane <-fit_herb_det_SS2_treatments$params$Xcoef[,1]
effect_SS2_control_sugar_cane_se <-fit_herb_det_SS2_treatments$sd$Xcoef[,1]
#effect_SS2_control_sugar_cane_se[which(effect_SS2_control_sugar_cane_se>100)] <- fit_herb_det_SS2_treatments_control$sd$Xcoef[,2][which(effect_SS2_control_sugar_cane_se>100)]

#effect_SS2_control_sugar_cane[which(effect_SS2_control_sugar_cane_se>100)] <- 0
effect_SS2_control_sugar_cane_se[which(effect_SS2_control_sugar_cane_se>100)] <- 0

effect_SS2_control_sugar_cane_upper <- effect_SS2_control_sugar_cane + (z*effect_SS2_control_sugar_cane_se)
effect_SS2_control_sugar_cane_lower <- effect_SS2_control_sugar_cane - (z*effect_SS2_control_sugar_cane_se)


ab_SS2_herb <- order(Trait_SS2$total_ab[which(Trait_SS2$trophic == "consumer")], decreasing = T)
effect_SS2_control_sugar_cane <- effect_SS2_control_sugar_cane[ab_SS2_herb]
effect_SS2_control_sugar_cane_lower <- effect_SS2_control_sugar_cane_lower[ab_SS2_herb]
effect_SS2_control_sugar_cane_upper <- effect_SS2_control_sugar_cane_upper[ab_SS2_herb]




############################################################################################################################
###################################################### SS3 #################################################################
############################################################################################################################

#effect_SS3 of TREATMENTS
#30

fit_herb_det_SS3_treatments_control <- gllvm(com_herb_det_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                              treatments = relevel(treatments_SS2_3_4_contpast_sug, ref = "control_pasture"),
                                                                              isolation = relevel(isolation_SS2_3_4, ref = "30")),
                                             formula = ~ treatments * isolation,
                                             family = "negative.binomial",
                                             method = "VA",
                                             row.eff = ~ (1|ID),
                                             n.init = 20, num.lv = 0, seed = 1:20, control = list(optimizer = "nlminb"))



effect_SS3_30_control_sugar_cane <-fit_herb_det_SS3_treatments_control$params$Xcoef[,1]
effect_SS3_30_control_sugar_cane_se <-fit_herb_det_SS3_treatments_control$sd$Xcoef[,1]
#effect_SS3_30_control_sugar_cane_se[which(effect_SS3_30_control_sugar_cane_se>100)] <- fit_herb_det_SS3_treatments_control$sd$Xcoef[,2][which(effect_SS3_30_control_sugar_cane_se>100)]

#effect_SS3_30_control_sugar_cane[which(effect_SS3_30_control_sugar_cane_se>100)] <- 0
effect_SS3_30_control_sugar_cane_se[which(effect_SS3_30_control_sugar_cane_se>100)] <- 0


effect_SS3_30_control_sugar_cane_lower <- effect_SS3_30_control_sugar_cane - (z*effect_SS3_30_control_sugar_cane_se)
effect_SS3_30_control_sugar_cane_upper <- effect_SS3_30_control_sugar_cane + (z*effect_SS3_30_control_sugar_cane_se)


ab_SS3_herb <- order(Trait_SS3$total_ab[which(Trait_SS3$trophic == "consumer")], decreasing = T)

effect_SS3_30_control_sugar_cane <- effect_SS3_30_control_sugar_cane[ab_SS3_herb]
effect_SS3_30_control_sugar_cane_lower <- effect_SS3_30_control_sugar_cane_lower[ab_SS3_herb]
effect_SS3_30_control_sugar_cane_upper <- effect_SS3_30_control_sugar_cane_upper[ab_SS3_herb]



#120

fit_herb_det_SS3_treatments_control <- gllvm(com_herb_det_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                              treatments = relevel(treatments_SS2_3_4_contpast_sug, ref = "control_pasture"),
                                                                              isolation = relevel(isolation_SS2_3_4, ref = "120")),
                                             formula = ~ treatments * isolation,
                                             family = "negative.binomial",
                                             method = "VA",
                                             row.eff = ~ (1|ID),
                                             n.init = 20, num.lv = 0, seed = 1:20, control = list(optimizer = "nlminb"))



effect_SS3_120_control_sugar_cane <-fit_herb_det_SS3_treatments_control$params$Xcoef[,1]
effect_SS3_120_control_sugar_cane_se <-fit_herb_det_SS3_treatments_control$sd$Xcoef[,1]
#effect_SS3_120_control_sugar_cane_se[which(effect_SS3_120_control_sugar_cane_se>100)] <- fit_herb_det_SS3_treatments_control$sd$Xcoef[,2][which(effect_SS3_120_control_sugar_cane_se>100)]

#effect_SS3_120_control_sugar_cane[which(effect_SS3_120_control_sugar_cane_se>100)] <- 0
effect_SS3_120_control_sugar_cane_se[which(effect_SS3_120_control_sugar_cane_se>100)] <- 0


effect_SS3_120_control_sugar_cane_lower <- effect_SS3_120_control_sugar_cane - (z*effect_SS3_120_control_sugar_cane_se)
effect_SS3_120_control_sugar_cane_upper <- effect_SS3_120_control_sugar_cane + (z*effect_SS3_120_control_sugar_cane_se)


ab_SS3_herb <- order(Trait_SS3$total_ab[which(Trait_SS3$trophic == "consumer")], decreasing = T)

effect_SS3_120_control_sugar_cane <- effect_SS3_120_control_sugar_cane[ab_SS3_herb]
effect_SS3_120_control_sugar_cane_lower <- effect_SS3_120_control_sugar_cane_lower[ab_SS3_herb]
effect_SS3_120_control_sugar_cane_upper <- effect_SS3_120_control_sugar_cane_upper[ab_SS3_herb]




#480

fit_herb_det_SS3_treatments_control <- gllvm(com_herb_det_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                              treatments = relevel(treatments_SS2_3_4_contpast_sug, ref = "control_pasture"),
                                                                              isolation = relevel(isolation_SS2_3_4, ref = "480")),
                                             formula = ~ treatments * isolation,
                                             family = "negative.binomial",
                                             method = "VA",
                                             row.eff = ~ (1|ID),
                                             n.init = 20, num.lv = 0, seed = 1:20, control = list(optimizer = "nlminb"))


effect_SS3_480_control_sugar_cane <-fit_herb_det_SS3_treatments_control$params$Xcoef[,1]
effect_SS3_480_control_sugar_cane_se <-fit_herb_det_SS3_treatments_control$sd$Xcoef[,1]
#effect_SS3_480_control_sugar_cane_se[which(effect_SS3_480_control_sugar_cane_se>100)] <- fit_herb_det_SS3_treatments_control$sd$Xcoef[,2][which(effect_SS3_480_control_sugar_cane_se>100)]

#effect_SS3_480_control_sugar_cane[which(effect_SS3_480_control_sugar_cane_se>100)] <- 0
effect_SS3_480_control_sugar_cane_se[which(effect_SS3_480_control_sugar_cane_se>100)] <- 0


effect_SS3_480_control_sugar_cane_lower <- effect_SS3_480_control_sugar_cane - (z*effect_SS3_480_control_sugar_cane_se)
effect_SS3_480_control_sugar_cane_upper <- effect_SS3_480_control_sugar_cane + (z*effect_SS3_480_control_sugar_cane_se)


ab_SS3_herb <- order(Trait_SS3$total_ab[which(Trait_SS3$trophic == "consumer")], decreasing = T)

effect_SS3_480_control_sugar_cane <- effect_SS3_480_control_sugar_cane[ab_SS3_herb]
effect_SS3_480_control_sugar_cane_lower <- effect_SS3_480_control_sugar_cane_lower[ab_SS3_herb]
effect_SS3_480_control_sugar_cane_upper <- effect_SS3_480_control_sugar_cane_upper[ab_SS3_herb]




############################################################################################################################
###################################################### SS4 #################################################################
############################################################################################################################



#effect_SS4 of TREATMENTS
#30

fit_herb_det_SS4_treatments_control <- gllvm(com_herb_det_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                              treatments = relevel(treatments_SS2_3_4, ref = "control"),
                                                                              isolation = relevel(isolation_SS2_3_4, ref = "30")),
                                             formula = ~ treatments * isolation,
                                             family = "negative.binomial",
                                             method = "VA",
                                             row.eff = ~ (1|ID),
                                             n.init = 20, num.lv = 0, seed = 1:20)



effect_SS4_30_control_pasture <-fit_herb_det_SS4_treatments_control$params$Xcoef[,1]
effect_SS4_30_control_pasture_se <-fit_herb_det_SS4_treatments_control$sd$Xcoef[,1]
#effect_SS4_30_control_pasture_se[which(effect_SS4_30_control_pasture_se>100)] <- fit_herb_det_SS4_treatments_control$sd$Xcoef[,2][which(effect_SS4_30_control_pasture_se>100)]

#effect_SS4_30_control_pasture[which(effect_SS4_30_control_pasture_se>100)] <- 0
effect_SS4_30_control_pasture_se[which(effect_SS4_30_control_pasture_se>100)] <- 0


effect_SS4_30_control_pasture_lower <- effect_SS4_30_control_pasture - (z*effect_SS4_30_control_pasture_se)
effect_SS4_30_control_pasture_upper <- effect_SS4_30_control_pasture + (z*effect_SS4_30_control_pasture_se)



effect_SS4_30_control_sugar_cane <-fit_herb_det_SS4_treatments_control$params$Xcoef[,2]
effect_SS4_30_control_sugar_cane_se <-fit_herb_det_SS4_treatments_control$sd$Xcoef[,2]
#effect_SS4_30_control_sugar_cane_se[which(effect_SS4_30_control_sugar_cane_se>100)] <- fit_herb_det_SS4_treatments_control$sd$Xcoef[,2][which(effect_SS4_30_control_sugar_cane_se>100)]

#effect_SS4_30_control_sugar_cane[which(effect_SS4_30_control_sugar_cane_se>100)] <- 0
effect_SS4_30_control_sugar_cane_se[which(effect_SS4_30_control_sugar_cane_se>100)] <- 0


effect_SS4_30_control_sugar_cane_lower <- effect_SS4_30_control_sugar_cane - (z*effect_SS4_30_control_sugar_cane_se)
effect_SS4_30_control_sugar_cane_upper <- effect_SS4_30_control_sugar_cane + (z*effect_SS4_30_control_sugar_cane_se)




fit_herb_det_SS4_treatments_pasture <- gllvm(com_herb_det_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                              treatments = relevel(treatments_SS2_3_4, ref = "pasture"),
                                                                              isolation = relevel(isolation_SS2_3_4, ref = "30")),
                                             formula = ~ treatments * isolation,
                                             family = "negative.binomial",
                                             method = "VA",
                                             row.eff = ~ (1|ID),
                                             n.init = 20, num.lv = 0, seed = 1:20)


effect_SS4_30_pasture_sugar_cane <-fit_herb_det_SS4_treatments_pasture$params$Xcoef[,2]
effect_SS4_30_pasture_sugar_cane_se <-fit_herb_det_SS4_treatments_pasture$sd$Xcoef[,2]
#effect_SS4_30_pasture_sugar_cane_se[which(effect_SS4_30_pasture_sugar_cane_se>100)] <- fit_herb_det_SS4_treatments_pasture$sd$Xcoef[,2][which(effect_SS4_30_pasture_sugar_cane_se>100)]

#effect_SS4_30_pasture_sugar_cane[which(effect_SS4_30_pasture_sugar_cane_se>100)] <- 0
effect_SS4_30_pasture_sugar_cane_se[which(effect_SS4_30_pasture_sugar_cane_se>100)] <- 0


effect_SS4_30_pasture_sugar_cane_lower <- effect_SS4_30_pasture_sugar_cane - (z*effect_SS4_30_pasture_sugar_cane_se)
effect_SS4_30_pasture_sugar_cane_upper <- effect_SS4_30_pasture_sugar_cane + (z*effect_SS4_30_pasture_sugar_cane_se)








ab_SS4_herb <- order(Trait_SS4$total_ab[which(Trait_SS4$trophic == "consumer")], decreasing = T)

effect_SS4_30_control_pasture <- effect_SS4_30_control_pasture[ab_SS4_herb]
effect_SS4_30_control_pasture_lower <- effect_SS4_30_control_pasture_lower[ab_SS4_herb]
effect_SS4_30_control_pasture_upper <- effect_SS4_30_control_pasture_upper[ab_SS4_herb]

effect_SS4_30_control_sugar_cane <- effect_SS4_30_control_sugar_cane[ab_SS4_herb]
effect_SS4_30_control_sugar_cane_lower <- effect_SS4_30_control_sugar_cane_lower[ab_SS4_herb]
effect_SS4_30_control_sugar_cane_upper <- effect_SS4_30_control_sugar_cane_upper[ab_SS4_herb]


effect_SS4_30_pasture_sugar_cane <- effect_SS4_30_pasture_sugar_cane[ab_SS4_herb]
effect_SS4_30_pasture_sugar_cane_lower <- effect_SS4_30_pasture_sugar_cane_lower[ab_SS4_herb]
effect_SS4_30_pasture_sugar_cane_upper <- effect_SS4_30_pasture_sugar_cane_upper[ab_SS4_herb]







#120

fit_herb_det_SS4_treatments_control <- gllvm(com_herb_det_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                              treatments = relevel(treatments_SS2_3_4_contpast_sug, ref = "control_pasture"),
                                                                              isolation = relevel(isolation_SS2_3_4, ref = "120")),
                                             formula = ~ treatments * isolation,
                                             family = "negative.binomial",
                                             method = "VA",
                                             row.eff = ~ (1|ID),
                                             n.init = 20, num.lv = 0, seed = 1:20)



effect_SS4_120_control_sugar_cane <-fit_herb_det_SS4_treatments_control$params$Xcoef[,1]
effect_SS4_120_control_sugar_cane_se <-fit_herb_det_SS4_treatments_control$sd$Xcoef[,1]
#effect_SS4_120_control_sugar_cane_se[which(effect_SS4_120_control_sugar_cane_se>100)] <- fit_herb_det_SS4_treatments_control$sd$Xcoef[,2][which(effect_SS4_120_control_sugar_cane_se>100)]

#effect_SS4_120_control_sugar_cane[which(effect_SS4_120_control_sugar_cane_se>100)] <- 0
effect_SS4_120_control_sugar_cane_se[which(effect_SS4_120_control_sugar_cane_se>100)] <- 0


effect_SS4_120_control_sugar_cane_lower <- effect_SS4_120_control_sugar_cane - (z*effect_SS4_120_control_sugar_cane_se)
effect_SS4_120_control_sugar_cane_upper <- effect_SS4_120_control_sugar_cane + (z*effect_SS4_120_control_sugar_cane_se)


ab_SS4_herb <- order(Trait_SS4$total_ab[which(Trait_SS4$trophic == "consumer")], decreasing = T)

effect_SS4_120_control_sugar_cane <- effect_SS4_120_control_sugar_cane[ab_SS4_herb]
effect_SS4_120_control_sugar_cane_lower <- effect_SS4_120_control_sugar_cane_lower[ab_SS4_herb]
effect_SS4_120_control_sugar_cane_upper <- effect_SS4_120_control_sugar_cane_upper[ab_SS4_herb]







































############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################ ISOLATION
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################




############################################################################################################################
###################################################### SS1 #################################################################
############################################################################################################################

#effect_SS1 of  isolation




fit_herb_det_SS1_isolation_30480 <- gllvm(com_herb_det_SS1, X = data.frame(ID = ID_SS1,
                                                                            isolation = relevel(isolation30_120480, ref = "30")),
                                           formula = ~ isolation,
                                           family = "negative.binomial",
                                           method = "VA",
                                           #row.eff = ~ (1|ID),
                                           n.init = 20, num.lv = 0, seed = 1:20)





effect_SS1_30480 <-fit_herb_det_SS1_isolation_30480$params$Xcoef[,1]
effect_SS1_30480_se <-fit_herb_det_SS1_isolation_30480$sd$Xcoef[,1]
#effect_SS1_30480_se[which(effect_SS1_30480_se>100)] <- fit_herb_det_SS1_isolation_30480$sd$Xcoef[,2][which(effect_SS1_30480_se>100)]

#effect_SS1_30480[which(effect_SS1_30480_se>100)] <- 0
effect_SS1_30480_se[which(effect_SS1_30480_se>100)] <- 0




effect_SS1_30480_lower <- effect_SS1_30480 - (z*effect_SS1_30480_se)
effect_SS1_30480_upper <- effect_SS1_30480 + (z*effect_SS1_30480_se)



ab_SS1_herb <- order(Trait_SS1$total_ab[which(Trait_SS1$trophic == "consumer")], decreasing = T)

effect_SS1_30480 <- effect_SS1_30480[ab_SS1_herb]
effect_SS1_30480_lower <- effect_SS1_30480_lower[ab_SS1_herb]
effect_SS1_30480_upper <- effect_SS1_30480_upper[ab_SS1_herb]









############################################################################################################################
###################################################### SS2 #################################################################
############################################################################################################################

#effect_SS2 of  isolation


fit_herb_det_SS2_isolation_30_120 <- gllvm(com_herb_det_SS2, X = data.frame(ID = ID_SS2_3_4,
                                                                            isolation = relevel(isolation_SS2_3_4, ref = "30"),
                                                                            treatments = treatments_SS2_3_4_contpast_sug),
                                           formula = ~ isolation + treatments,
                                           family = "negative.binomial",
                                           method = "VA",
                                           row.eff = ~ (1|ID),
                                           n.init = 20, num.lv = 0, seed = 1:20)


effect_SS2_30_120 <-fit_herb_det_SS2_isolation_30_120$params$Xcoef[,1]
effect_SS2_30_120_se <-fit_herb_det_SS2_isolation_30_120$sd$Xcoef[,1]
#effect_SS2_30_120_480_se[which(effect_SS2_30_120_480_se>100)] <- fit_herb_det_SS2_isolation_30_120$sd$Xcoef[,2][which(effect_SS2_30_120_480_se>100)]

#effect_SS2_30_120_480[which(effect_SS2_30_120_480_se>100)] <- 0
effect_SS2_30_120_se[which(effect_SS2_30_120_se>100)] <- 0




effect_SS2_30_120_lower <- effect_SS2_30_120 - (z*effect_SS2_30_120_se)
effect_SS2_30_120_upper <- effect_SS2_30_120 + (z*effect_SS2_30_120_se)



ab_SS2_herb <- order(Trait_SS2$total_ab[which(Trait_SS2$trophic == "consumer")], decreasing = T)

effect_SS2_30_120 <- effect_SS2_30_120[ab_SS2_herb]
effect_SS2_30_120_lower <- effect_SS2_30_120_lower[ab_SS2_herb]
effect_SS2_30_120_upper <- effect_SS2_30_120_upper[ab_SS2_herb]






effect_SS2_30_480 <-fit_herb_det_SS2_isolation_30_120$params$Xcoef[,2]
effect_SS2_30_480_se <-fit_herb_det_SS2_isolation_30_120$sd$Xcoef[,2]
#effect_SS2_30_480_480_se[which(effect_SS2_30_480_480_se>100)] <- fit_herb_det_SS2_isolation_30_480$sd$Xcoef[,2][which(effect_SS2_30_480_480_se>100)]

#effect_SS2_30_480_480[which(effect_SS2_30_480_480_se>100)] <- 0
effect_SS2_30_480_se[which(effect_SS2_30_480_se>100)] <- 0




effect_SS2_30_480_lower <- effect_SS2_30_480 - (z*effect_SS2_30_480_se)
effect_SS2_30_480_upper <- effect_SS2_30_480 + (z*effect_SS2_30_480_se)



ab_SS2_herb <- order(Trait_SS2$total_ab[which(Trait_SS2$trophic == "consumer")], decreasing = T)

effect_SS2_30_480 <- effect_SS2_30_480[ab_SS2_herb]
effect_SS2_30_480_lower <- effect_SS2_30_480_lower[ab_SS2_herb]
effect_SS2_30_480_upper <- effect_SS2_30_480_upper[ab_SS2_herb]






fit_herb_det_SS2_isolation_120_480 <- gllvm(com_herb_det_SS2, X = data.frame(ID = ID_SS2_3_4,
                                                                            isolation = relevel(isolation_SS2_3_4, ref = "120"),
                                                                            treatments = treatments_SS2_3_4_contpast_sug),
                                           formula = ~ isolation + treatments,
                                           family = "negative.binomial",
                                           method = "VA",
                                           row.eff = ~ (1|ID),
                                           n.init = 20, num.lv = 0, seed = 1:20)


effect_SS2_120_480 <-fit_herb_det_SS2_isolation_120_480$params$Xcoef[,2]
effect_SS2_120_480_se <-fit_herb_det_SS2_isolation_120_480$sd$Xcoef[,2]
#effect_SS2_120_480_480_se[which(effect_SS2_120_480_480_se>100)] <- fit_herb_det_SS2_isolation_120_480$sd$Xcoef[,2][which(effect_SS2_120_480_480_se>100)]

#effect_SS2_120_480_480[which(effect_SS2_120_480_480_se>100)] <- 0
effect_SS2_120_480_se[which(effect_SS2_120_480_se>100)] <- 0




effect_SS2_120_480_lower <- effect_SS2_120_480 - (z*effect_SS2_120_480_se)
effect_SS2_120_480_upper <- effect_SS2_120_480 + (z*effect_SS2_120_480_se)



ab_SS2_herb <- order(Trait_SS2$total_ab[which(Trait_SS2$trophic == "consumer")], decreasing = T)

effect_SS2_120_480 <- effect_SS2_120_480[ab_SS2_herb]
effect_SS2_120_480_lower <- effect_SS2_120_480_lower[ab_SS2_herb]
effect_SS2_120_480_upper <- effect_SS2_120_480_upper[ab_SS2_herb]




############################################################################################################################
###################################################### SS3 #################################################################
############################################################################################################################


#effect_SS3 of  isolation - CONTROL


fit_herb_det_SS3_isolation_30_120_control <- gllvm(com_herb_det_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                            isolation = relevel(isolation30120_480_SS2_3_4, ref = "30_120"),
                                                                            treatments = relevel(treatments_SS2_3_4, ref = "control")),
                                           formula = ~ isolation * treatments,
                                           family = "negative.binomial",
                                           method = "VA",
                                           row.eff = ~ (1|ID),
                                           n.init = 20, num.lv = 0, seed = 1:20, control = list(optimizer = "nlminb"))



com_herb_det_SS3$Heterelmis[treatments_SS2_3_4 == "control" & isolation30120_480_SS2_3_4 == "30_120" ]
com_herb_det_SS3$Heterelmis[treatments_SS2_3_4 == "control" & isolation30120_480_SS2_3_4 == "480" ]



effect_SS3_30_120_480_control <-fit_herb_det_SS3_isolation_30_120_control$params$Xcoef[,1]
effect_SS3_30_120_480_control_se <-fit_herb_det_SS3_isolation_30_120_control$sd$Xcoef[,1]
#effect_SS3_30_120_480_control_se[which(effect_SS3_30_120_480_control_se>100)] <- fit_herb_det_SS3_isolation_30_120$sd$Xcoef[,2][which(effect_SS3_30_120_480_control_se>100)]

#effect_SS3_30_120_480_control[which(effect_SS3_30_120_480_control_se>100)] <- 0
effect_SS3_30_120_480_control_se[which(effect_SS3_30_120_480_control_se>100)] <- 0




effect_SS3_30_120_480_control_lower <- effect_SS3_30_120_480_control - (z*effect_SS3_30_120_480_control_se)
effect_SS3_30_120_480_control_upper <- effect_SS3_30_120_480_control + (z*effect_SS3_30_120_480_control_se)



ab_SS3_herb <- order(Trait_SS3$total_ab[which(Trait_SS3$trophic == "consumer")], decreasing = T)

effect_SS3_30_120_480_control <- effect_SS3_30_120_480_control[ab_SS3_herb]
effect_SS3_30_120_480_control_lower <- effect_SS3_30_120_480_control_lower[ab_SS3_herb]
effect_SS3_30_120_480_control_upper <- effect_SS3_30_120_480_control_upper[ab_SS3_herb]




##########################################################################

#effect_SS3 of  isolation - pasture


fit_herb_det_SS3_isolation_30_pasture <- gllvm(com_herb_det_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                                    isolation = relevel(isolation_SS2_3_4, ref = "30"),
                                                                                    treatments = relevel(treatments_SS2_3_4, ref = "pasture")),
                                                   formula = ~ isolation * treatments,
                                                   family = "negative.binomial",
                                                   method = "VA",
                                                   row.eff = ~ (1|ID),
                                                   n.init = 20, num.lv = 0, seed = 1:20, control = list(optimizer = "nlminb"))




effect_SS3_pasture_30_120 <-fit_herb_det_SS3_isolation_30_pasture$params$Xcoef[,1]
effect_SS3_pasture_30_120_se <-fit_herb_det_SS3_isolation_30_pasture$sd$Xcoef[,1]
#effect_SS3_pasture_30_120_se[which(effect_SS3_pasture_30_120_se>100)] <- fit_herb_det_SS3_isolation_30$sd$Xcoef[,1][which(effect_SS3_pasture_30_120_se>100)]

#effect_SS3_pasture_30_120[which(effect_SS3_pasture_30_120_se>100)] <- 0
effect_SS3_pasture_30_120_se[which(effect_SS3_pasture_30_120_se>100)] <- 0


effect_SS3_pasture_30_480 <-fit_herb_det_SS3_isolation_30_pasture$params$Xcoef[,2]
effect_SS3_pasture_30_480_se <-fit_herb_det_SS3_isolation_30_pasture$sd$Xcoef[,2]
#effect_SS3_pasture_30_480_se[which(effect_SS3_pasture_30_480_se>100)] <- fit_herb_det_SS3_isolation_30$sd$Xcoef[,2][which(effect_SS3_pasture_30_480_se>100)]

#effect_SS3_pasture_30_480[which(effect_SS3_pasture_30_480_se>100)] <- 0
effect_SS3_pasture_30_480_se[which(effect_SS3_pasture_30_480_se>100)] <- 0

fit_herb_det_SS3_isolation_120_pasture <- gllvm(com_herb_det_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                                    isolation = relevel(isolation_SS2_3_4, ref = "120"),
                                                                                    treatments = relevel(treatments_SS2_3_4, ref = "pasture")),
                                                   formula = ~ isolation * treatments,
                                                   family = "negative.binomial",
                                                   method = "VA",
                                                   row.eff = ~ (1|ID),
                                                   n.init = 20, num.lv = 0, seed = 1:20, control = list(optimizer = "nlminb"))


effect_SS3_pasture_120_480 <-fit_herb_det_SS3_isolation_120_pasture$params$Xcoef[,2]
effect_SS3_pasture_120_480_se <-fit_herb_det_SS3_isolation_120_pasture$sd$Xcoef[,2]
#effect_SS3_pasture_120_480_se[which(effect_SS3_pasture_120_480_se>100)] <- fit_herb_det_SS3_isolation_120$sd$Xcoef[,1][which(effect_SS3_pasture_120_480_se>100)]

#effect_SS3_pasture_120_480[which(effect_SS3_pasture_120_480_se>100)] <- 0
effect_SS3_pasture_120_480_se[which(effect_SS3_pasture_120_480_se>30)] <- 0

effect_SS3_pasture_30_120_lower <- effect_SS3_pasture_30_120 - (z*effect_SS3_pasture_30_120_se)
effect_SS3_pasture_30_120_upper <- effect_SS3_pasture_30_120 + (z*effect_SS3_pasture_30_120_se)

effect_SS3_pasture_30_480_lower <- effect_SS3_pasture_30_480 - (z*effect_SS3_pasture_30_480_se)
effect_SS3_pasture_30_480_upper <- effect_SS3_pasture_30_480 + (z*effect_SS3_pasture_30_480_se)

effect_SS3_pasture_120_480_lower <- effect_SS3_pasture_120_480 - (z*effect_SS3_pasture_120_480_se)
effect_SS3_pasture_120_480_upper <- effect_SS3_pasture_120_480 + (z*effect_SS3_pasture_120_480_se)



ab_SS3_herb <- order(Trait_SS3$total_ab[which(Trait_SS3$trophic == "consumer")], decreasing = T)

effect_SS3_pasture_30_120 <- effect_SS3_pasture_30_120[ab_SS3_herb]
effect_SS3_pasture_30_120_lower <- effect_SS3_pasture_30_120_lower[ab_SS3_herb]
effect_SS3_pasture_30_120_upper <- effect_SS3_pasture_30_120_upper[ab_SS3_herb]

effect_SS3_pasture_30_480 <- effect_SS3_pasture_30_480[ab_SS3_herb]
effect_SS3_pasture_30_480_lower <- effect_SS3_pasture_30_480_lower[ab_SS3_herb]
effect_SS3_pasture_30_480_upper <- effect_SS3_pasture_30_480_upper[ab_SS3_herb]

effect_SS3_pasture_120_480 <- effect_SS3_pasture_120_480[ab_SS3_herb]

effect_SS3_pasture_120_480_lower <- effect_SS3_pasture_120_480_lower[ab_SS3_herb]
effect_SS3_pasture_120_480_upper <- effect_SS3_pasture_120_480_upper[ab_SS3_herb]

effect_SS3_pasture_120_480[10] <- 0
effect_SS3_pasture_120_480_lower[10] <- 0
effect_SS3_pasture_120_480_upper[10] <- 0








#########################################################
#effect_SS3 of  isolation - sugar_cane


fit_herb_det_SS3_isolation_30_120_sugar_cane <- gllvm(com_herb_det_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                                    isolation = relevel(isolation30120_480_SS2_3_4, ref = "30_120"),
                                                                                    treatments = relevel(treatments_SS2_3_4, ref = "sugar_cane")),
                                                   formula = ~ isolation * treatments,
                                                   family = "negative.binomial",
                                                   method = "VA",
                                                   row.eff = ~ (1|ID),
                                                   n.init = 20, num.lv = 0, seed = 1:20, control = list(optimizer = "nlminb"))





effect_SS3_30_120_480_sugar_cane <-fit_herb_det_SS3_isolation_30_120_sugar_cane$params$Xcoef[,1]
effect_SS3_30_120_480_sugar_cane_se <-fit_herb_det_SS3_isolation_30_120_sugar_cane$sd$Xcoef[,1]
#effect_SS3_30_120_480_sugar_cane_se[which(effect_SS3_30_120_480_sugar_cane_se>100)] <- fit_herb_det_SS3_isolation_30_120$sd$Xcoef[,2][which(effect_SS3_30_120_480_sugar_cane_se>100)]

#effect_SS3_30_120_480_sugar_cane[which(effect_SS3_30_120_480_sugar_cane_se>100)] <- 0
effect_SS3_30_120_480_sugar_cane_se[which(effect_SS3_30_120_480_sugar_cane_se>100)] <- 0




effect_SS3_30_120_480_sugar_cane_lower <- effect_SS3_30_120_480_sugar_cane - (z*effect_SS3_30_120_480_sugar_cane_se)
effect_SS3_30_120_480_sugar_cane_upper <- effect_SS3_30_120_480_sugar_cane + (z*effect_SS3_30_120_480_sugar_cane_se)



ab_SS3_herb <- order(Trait_SS3$total_ab[which(Trait_SS3$trophic == "consumer")], decreasing = T)

effect_SS3_30_120_480_sugar_cane <- effect_SS3_30_120_480_sugar_cane[ab_SS3_herb]
effect_SS3_30_120_480_sugar_cane_lower <- effect_SS3_30_120_480_sugar_cane_lower[ab_SS3_herb]
effect_SS3_30_120_480_sugar_cane_upper <- effect_SS3_30_120_480_sugar_cane_upper[ab_SS3_herb]










############################################################################################################################
###################################################### SS4 #################################################################
############################################################################################################################




#effect_SS4 of  isolation - CONTROL


fit_herb_det_SS4_isolation_30_120_control <- gllvm(com_herb_det_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                                    isolation = relevel(isolation30_120480_SS2_3_4, ref = "30"),
                                                                                    treatments = relevel(treatments_SS2_3_4, ref = "control")),
                                                   formula = ~ isolation * treatments,
                                                   family = "negative.binomial",
                                                   method = "VA",
                                                   row.eff = ~ (1|ID),
                                                   n.init = 20, num.lv = 0, seed = 1:20)





effect_SS4_30_120_480_control <-fit_herb_det_SS4_isolation_30_120_control$params$Xcoef[,1]
effect_SS4_30_120_480_control_se <-fit_herb_det_SS4_isolation_30_120_control$sd$Xcoef[,1]
#effect_SS4_30_120_480_control_se[which(effect_SS4_30_120_480_control_se>100)] <- fit_herb_det_SS4_isolation_30_120$sd$Xcoef[,2][which(effect_SS4_30_120_480_control_se>100)]

#effect_SS4_30_120_480_control[which(effect_SS4_30_120_480_control_se>100)] <- 0
effect_SS4_30_120_480_control_se[which(effect_SS4_30_120_480_control_se>100)] <- 0




effect_SS4_30_120_480_control_lower <- effect_SS4_30_120_480_control - (z*effect_SS4_30_120_480_control_se)
effect_SS4_30_120_480_control_upper <- effect_SS4_30_120_480_control + (z*effect_SS4_30_120_480_control_se)



ab_SS4_herb <- order(Trait_SS4$total_ab[which(Trait_SS4$trophic == "consumer")], decreasing = T)

effect_SS4_30_120_480_control <- effect_SS4_30_120_480_control[ab_SS4_herb]
effect_SS4_30_120_480_control_lower <- effect_SS4_30_120_480_control_lower[ab_SS4_herb]
effect_SS4_30_120_480_control_upper <- effect_SS4_30_120_480_control_upper[ab_SS4_herb]



##########################################################################

#effect_SS4 of  isolation - pasture


fit_herb_det_SS4_isolation_30_pasture <- gllvm(com_herb_det_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                                isolation = relevel(isolation_SS2_3_4, ref = "30"),
                                                                                treatments = relevel(treatments_SS2_3_4, ref = "pasture")),
                                               formula = ~ isolation * treatments,
                                               family = "negative.binomial",
                                               method = "VA",
                                               row.eff = ~ (1|ID),
                                               n.init = 20, num.lv = 0, seed = 1:20)




effect_SS4_pasture_30_120 <-fit_herb_det_SS4_isolation_30_pasture$params$Xcoef[,1]
effect_SS4_pasture_30_120_se <-fit_herb_det_SS4_isolation_30_pasture$sd$Xcoef[,1]
#effect_SS4_pasture_30_120_se[which(effect_SS4_pasture_30_120_se>100)] <- fit_herb_det_SS4_isolation_30$sd$Xcoef[,1][which(effect_SS4_pasture_30_120_se>100)]

#effect_SS4_pasture_30_120[which(effect_SS4_pasture_30_120_se>100)] <- 0
effect_SS4_pasture_30_120_se[which(effect_SS4_pasture_30_120_se>100)] <- 0


effect_SS4_pasture_30_480 <-fit_herb_det_SS4_isolation_30_pasture$params$Xcoef[,2]
effect_SS4_pasture_30_480_se <-fit_herb_det_SS4_isolation_30_pasture$sd$Xcoef[,2]
#effect_SS4_pasture_30_480_se[which(effect_SS4_pasture_30_480_se>100)] <- fit_herb_det_SS4_isolation_30$sd$Xcoef[,2][which(effect_SS4_pasture_30_480_se>100)]

#effect_SS4_pasture_30_480[which(effect_SS4_pasture_30_480_se>100)] <- 0
effect_SS4_pasture_30_480_se[which(effect_SS4_pasture_30_480_se>100)] <- 0

fit_herb_det_SS4_isolation_120_pasture <- gllvm(com_herb_det_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                                 isolation = relevel(isolation_SS2_3_4, ref = "120"),
                                                                                 treatments = relevel(treatments_SS2_3_4, ref = "pasture")),
                                                formula = ~ isolation * treatments,
                                                family = "negative.binomial",
                                                method = "VA",
                                                row.eff = ~ (1|ID),
                                                n.init = 20, num.lv = 0, seed = 1:20)


effect_SS4_pasture_120_480 <-fit_herb_det_SS4_isolation_120_pasture$params$Xcoef[,2]
effect_SS4_pasture_120_480_se <-fit_herb_det_SS4_isolation_120_pasture$sd$Xcoef[,2]
#effect_SS4_pasture_120_480_se[which(effect_SS4_pasture_120_480_se>100)] <- fit_herb_det_SS4_isolation_120$sd$Xcoef[,1][which(effect_SS4_pasture_120_480_se>100)]

#effect_SS4_pasture_120_480[which(effect_SS4_pasture_120_480_se>100)] <- 0
effect_SS4_pasture_120_480_se[which(effect_SS4_pasture_120_480_se>30)] <- 0

effect_SS4_pasture_30_120_lower <- effect_SS4_pasture_30_120 - (z*effect_SS4_pasture_30_120_se)
effect_SS4_pasture_30_120_upper <- effect_SS4_pasture_30_120 + (z*effect_SS4_pasture_30_120_se)

effect_SS4_pasture_30_480_lower <- effect_SS4_pasture_30_480 - (z*effect_SS4_pasture_30_480_se)
effect_SS4_pasture_30_480_upper <- effect_SS4_pasture_30_480 + (z*effect_SS4_pasture_30_480_se)

effect_SS4_pasture_120_480_lower <- effect_SS4_pasture_120_480 - (z*effect_SS4_pasture_120_480_se)
effect_SS4_pasture_120_480_upper <- effect_SS4_pasture_120_480 + (z*effect_SS4_pasture_120_480_se)



ab_SS4_herb <- order(Trait_SS4$total_ab[which(Trait_SS4$trophic == "consumer")], decreasing = T)

effect_SS4_pasture_30_120 <- effect_SS4_pasture_30_120[ab_SS4_herb]
effect_SS4_pasture_30_120_lower <- effect_SS4_pasture_30_120_lower[ab_SS4_herb]
effect_SS4_pasture_30_120_upper <- effect_SS4_pasture_30_120_upper[ab_SS4_herb]

effect_SS4_pasture_30_480 <- effect_SS4_pasture_30_480[ab_SS4_herb]
effect_SS4_pasture_30_480_lower <- effect_SS4_pasture_30_480_lower[ab_SS4_herb]
effect_SS4_pasture_30_480_upper <- effect_SS4_pasture_30_480_upper[ab_SS4_herb]

effect_SS4_pasture_120_480 <- effect_SS4_pasture_120_480[ab_SS4_herb]
effect_SS4_pasture_120_480_lower <- effect_SS4_pasture_120_480_lower[ab_SS4_herb]
effect_SS4_pasture_120_480_upper <- effect_SS4_pasture_120_480_upper[ab_SS4_herb]








#########################################################
#effect_SS4 of  isolation - sugar_cane


fit_herb_det_SS4_isolation_30_120_sugar_cane <- gllvm(com_herb_det_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                                       isolation = relevel(isolation30_120480_SS2_3_4, ref = "30"),
                                                                                       treatments = relevel(treatments_SS2_3_4, ref = "sugar_cane")),
                                                      formula = ~ isolation * treatments,
                                                      family = "negative.binomial",
                                                      method = "VA",
                                                      row.eff = ~ (1|ID),
                                                      n.init = 20, num.lv = 0, seed = 1:20)





effect_SS4_30_120_480_sugar_cane <-fit_herb_det_SS4_isolation_30_120_sugar_cane$params$Xcoef[,1]
effect_SS4_30_120_480_sugar_cane_se <-fit_herb_det_SS4_isolation_30_120_sugar_cane$sd$Xcoef[,1]
#effect_SS4_30_120_480_sugar_cane_se[which(effect_SS4_30_120_480_sugar_cane_se>100)] <- fit_herb_det_SS4_isolation_30_120$sd$Xcoef[,2][which(effect_SS4_30_120_480_sugar_cane_se>100)]

#effect_SS4_30_120_480_sugar_cane[which(effect_SS4_30_120_480_sugar_cane_se>100)] <- 0
effect_SS4_30_120_480_sugar_cane_se[which(effect_SS4_30_120_480_sugar_cane_se>100)] <- 0




effect_SS4_30_120_480_sugar_cane_lower <- effect_SS4_30_120_480_sugar_cane - (z*effect_SS4_30_120_480_sugar_cane_se)
effect_SS4_30_120_480_sugar_cane_upper <- effect_SS4_30_120_480_sugar_cane + (z*effect_SS4_30_120_480_sugar_cane_se)



ab_SS4_herb <- order(Trait_SS4$total_ab[which(Trait_SS4$trophic == "consumer")], decreasing = T)

effect_SS4_30_120_480_sugar_cane <- effect_SS4_30_120_480_sugar_cane[ab_SS4_herb]
effect_SS4_30_120_480_sugar_cane_lower <- effect_SS4_30_120_480_sugar_cane_lower[ab_SS4_herb]
effect_SS4_30_120_480_sugar_cane_upper <- effect_SS4_30_120_480_sugar_cane_upper[ab_SS4_herb]






































