

com_predators <- com[,Trait$trophic == "predator"]
com_predators_SS1 <- com_SS1[,Trait_SS1$trophic == "predator"]
com_predators_SS2 <- com_SS2[,Trait_SS2$trophic == "predator"]
com_predators_SS3 <- com_SS3[,Trait_SS3$trophic == "predator"]
com_predators_SS4 <- com_SS4[,Trait_SS4$trophic == "predator"]

z<-qnorm(0.975, mean = 0, sd = 1)


###################################################### All surveys #################################################################


fit_predators_survey_1 <- gllvm(com_predators, X = data.frame(ID = ID, survey = relevel(SS, ref = "1")),
                                              formula = ~ survey,
                                              family = "negative.binomial",
                                              method = "VA",
                                              row.eff = ~ (1|ID),
                                              n.init = 10, num.lv = 0, seed = 1:10)


effect_1_2 <-fit_predators_survey_1$params$Xcoef[,1]
effect_1_2_se <-fit_predators_survey_1$sd$Xcoef[,1]

effect_1_2_se[which(effect_1_2_se>50)] <- 0



fit_predators_survey_2 <- gllvm(com_predators, X = data.frame(ID = ID, survey = relevel(SS, ref = "2")),
                                formula = ~ survey,
                                family = "negative.binomial",
                                method = "VA",
                                row.eff = ~ (1|ID),
                                n.init = 10, num.lv = 0, seed = 1:10)


effect_2_3 <-fit_predators_survey_2$params$Xcoef[,2]
effect_2_3_se <-fit_predators_survey_2$sd$Xcoef[,2]

effect_2_3_se[which(effect_2_3_se>50)] <- 0



fit_predators_survey_3 <- gllvm(com_predators, X = data.frame(ID = ID, survey = relevel(SS, ref = "3")),
                                formula = ~ survey,
                                family = "negative.binomial",
                                method = "VA",
                                row.eff = ~ (1|ID),
                                n.init = 10, num.lv = 0, seed = 1:10)


effect_3_4 <-fit_predators_survey_3$params$Xcoef[,3]
effect_3_4_se <-fit_predators_survey_3$sd$Xcoef[,3]

effect_3_4_se[which(effect_3_4_se>50)] <- 0



effect_1_2_lower <- effect_1_2 - (z*effect_1_2_se)
effect_1_2_upper <- effect_1_2 + (z*effect_1_2_se)

effect_2_3_lower <- effect_2_3 - (z*effect_2_3_se)
effect_2_3_upper <- effect_2_3 + (z*effect_2_3_se)

effect_3_4_lower <- effect_3_4 - (z*effect_3_4_se)
effect_3_4_upper <- effect_3_4 + (z*effect_3_4_se)


ab_pred <- order(Trait$total_ab[which(Trait$trait == "insect_predator")], decreasing = T)

effect_1_2 <- effect_1_2[ab_pred]
effect_1_2_lower <- effect_1_2_lower[ab_pred]
effect_1_2_upper <- effect_1_2_upper[ab_pred]

effect_2_3 <- effect_2_3[ab_pred]
effect_2_3_lower <- effect_2_3_lower[ab_pred]
effect_2_3_upper <- effect_2_3_upper[ab_pred]

effect_3_4 <- effect_3_4[ab_pred]
effect_3_4_lower <- effect_3_4_lower[ab_pred]
effect_3_4_upper <- effect_3_4_upper[ab_pred]



###################################################### SS2 #################################################################
fit_predators_SS2_treatments<- gllvm(com_predators_SS2, X = data.frame(ID = ID_SS2_3_4,
                                                                                treatments = relevel(treatments_SS2_3_4_contpast_sug, ref = "control_pasture")),
                                              formula = ~ treatments,
                                              family = "negative.binomial",
                                              method = "VA",
                                              row.eff = ~ (1|ID),
                                              n.init = 10, num.lv = 0, seed = 11:20)



#effect_SS2 of TREATMENTS
effect_SS2_control_sugar_cane <-fit_predators_SS2_treatments$params$Xcoef[,1]
effect_SS2_control_sugar_cane_se <-fit_predators_SS2_treatments$sd$Xcoef[,1]
#effect_SS2_control_sugar_cane_se[which(effect_SS2_control_sugar_cane_se>100)] <- fit_predators_SS2_treatments_control$sd$Xcoef[,2][which(effect_SS2_control_sugar_cane_se>100)]

#effect_SS2_control_sugar_cane[which(effect_SS2_control_sugar_cane_se>100)] <- 0
effect_SS2_control_sugar_cane_se[which(effect_SS2_control_sugar_cane_se>100)] <- 0

effect_SS2_control_sugar_cane_upper <- effect_SS2_control_sugar_cane + (z*effect_SS2_control_sugar_cane_se)
effect_SS2_control_sugar_cane_lower <- effect_SS2_control_sugar_cane - (z*effect_SS2_control_sugar_cane_se)


ab_SS2_pred <- order(Trait_SS2$total_ab[which(Trait_SS2$trait == "insect_predator")], decreasing = T)
effect_SS2_control_sugar_cane <- effect_SS2_control_sugar_cane[ab_SS2_pred]
effect_SS2_control_sugar_cane_lower <- effect_SS2_control_sugar_cane_lower[ab_SS2_pred]
effect_SS2_control_sugar_cane_upper <- effect_SS2_control_sugar_cane_upper[ab_SS2_pred]





###################################################### SS3 - Treatments #################################################################

#effect_SS3 of TREATMENTS


fit_predators_SS3_treatments_control <- gllvm(com_predators_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                                treatments = relevel(treatments_SS2_3_4, ref = "control"),
                                                                                isolation = isolation_SS2_3_4),
                                              formula = ~ treatments + isolation,
                                              family = "negative.binomial",
                                              method = "VA",
                                              row.eff = ~ (1|ID),
                                              n.init = 10, num.lv = 0, seed = 11:20)




effect_SS3_control_pasture <-fit_predators_SS3_treatments_control$params$Xcoef[,1]
effect_SS3_control_pasture_se <-fit_predators_SS3_treatments_control$sd$Xcoef[,1]
#effect_SS3_control_pasture_se[which(effect_SS3_control_pasture_se>100)] <- fit_predators_SS3_treatments_control$sd$Xcoef[,1][which(effect_SS3_control_pasture_se>100)]

#effect_SS3_control_pasture[which(effect_SS3_control_pasture_se>100)] <- 0
effect_SS3_control_pasture_se[which(effect_SS3_control_pasture_se>100)] <- 0


effect_SS3_control_sugar_cane <-fit_predators_SS3_treatments_control$params$Xcoef[,2]
effect_SS3_control_sugar_cane_se <-fit_predators_SS3_treatments_control$sd$Xcoef[,2]
#effect_SS3_control_sugar_cane_se[which(effect_SS3_control_sugar_cane_se>100)] <- fit_predators_SS3_treatments_control$sd$Xcoef[,2][which(effect_SS3_control_sugar_cane_se>100)]

#effect_SS3_control_sugar_cane[which(effect_SS3_control_sugar_cane_se>100)] <- 0
effect_SS3_control_sugar_cane_se[which(effect_SS3_control_sugar_cane_se>100)] <- 0

fit_predators_SS3_treatments_pasture <- gllvm(com_predators_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                                treatments = relevel(treatments_SS2_3_4, ref = "pasture"),
                                                                                isolation = isolation_SS2_3_4),
                                              formula = ~ treatments + isolation,
                                              family = "negative.binomial",
                                              method = "VA",
                                              row.eff = ~ (1|ID),
                                              n.init = 10, num.lv = 0, seed = 11:20)


effect_SS3_pasture_sugar_cane <-fit_predators_SS3_treatments_pasture$params$Xcoef[,2]
effect_SS3_pasture_sugar_cane_se <-fit_predators_SS3_treatments_pasture$sd$Xcoef[,2]
#effect_SS3_pasture_sugar_cane_se[which(effect_SS3_pasture_sugar_cane_se>100)] <- fit_predators_SS3_treatments_pasture$sd$Xcoef[,1][which(effect_SS3_pasture_sugar_cane_se>100)]

#effect_SS3_pasture_sugar_cane[which(effect_SS3_pasture_sugar_cane_se>100)] <- 0
effect_SS3_pasture_sugar_cane_se[which(effect_SS3_pasture_sugar_cane_se>100)] <- 0

effect_SS3_control_pasture_lower <- effect_SS3_control_pasture - (z*effect_SS3_control_pasture_se)
effect_SS3_control_pasture_upper <- effect_SS3_control_pasture + (z*effect_SS3_control_pasture_se)

effect_SS3_control_sugar_cane_lower <- effect_SS3_control_sugar_cane - (z*effect_SS3_control_sugar_cane_se)
effect_SS3_control_sugar_cane_upper <- effect_SS3_control_sugar_cane + (z*effect_SS3_control_sugar_cane_se)

effect_SS3_pasture_sugar_cane_lower <- effect_SS3_pasture_sugar_cane - (z*effect_SS3_pasture_sugar_cane_se)
effect_SS3_pasture_sugar_cane_upper <- effect_SS3_pasture_sugar_cane + (z*effect_SS3_pasture_sugar_cane_se)



ab_SS3_pred <- order(Trait_SS3$total_ab[which(Trait_SS3$trait == "insect_predator")], decreasing = T)

effect_SS3_control_pasture <- effect_SS3_control_pasture[ab_SS3_pred]
effect_SS3_control_pasture_lower <- effect_SS3_control_pasture_lower[ab_SS3_pred]
effect_SS3_control_pasture_upper <- effect_SS3_control_pasture_upper[ab_SS3_pred]

effect_SS3_control_sugar_cane <- effect_SS3_control_sugar_cane[ab_SS3_pred]
effect_SS3_control_sugar_cane_lower <- effect_SS3_control_sugar_cane_lower[ab_SS3_pred]
effect_SS3_control_sugar_cane_upper <- effect_SS3_control_sugar_cane_upper[ab_SS3_pred]

effect_SS3_pasture_sugar_cane <- effect_SS3_pasture_sugar_cane[ab_SS3_pred]
effect_SS3_pasture_sugar_cane_lower <- effect_SS3_pasture_sugar_cane_lower[ab_SS3_pred]
effect_SS3_pasture_sugar_cane_upper <- effect_SS3_pasture_sugar_cane_upper[ab_SS3_pred]




###################################################### SS3 - Isolation #################################################################

#effect_SS3 of  isolation


fit_predators_SS3_isolation_30_120 <- gllvm(com_predators_SS3, X = data.frame(ID = ID_SS2_3_4,
                                                                                 isolation = relevel(isolation30120_480, ref = "30_120"),
                                                                              treatments = treatments_SS2_3_4),
                                              formula = ~ isolation + treatments,
                                              family = "negative.binomial",
                                              method = "VA",
                                              row.eff = ~ (1|ID),
                                              n.init = 10, num.lv = 0, seed = 11:20)





effect_SS3_30_120_480 <-fit_predators_SS3_isolation_30_120$params$Xcoef[,1]
effect_SS3_30_120_480_se <-fit_predators_SS3_isolation_30_120$sd$Xcoef[,1]
#effect_SS3_30_120_480_se[which(effect_SS3_30_120_480_se>100)] <- fit_predators_SS3_isolation_30_120$sd$Xcoef[,2][which(effect_SS3_30_120_480_se>100)]

#effect_SS3_30_120_480[which(effect_SS3_30_120_480_se>100)] <- 0
effect_SS3_30_120_480_se[which(effect_SS3_30_120_480_se>100)] <- 0




effect_SS3_30_120_480_lower <- effect_SS3_30_120_480 - (z*effect_SS3_30_120_480_se)
effect_SS3_30_120_480_upper <- effect_SS3_30_120_480 + (z*effect_SS3_30_120_480_se)



ab_SS3_pred <- order(Trait_SS3$total_ab[which(Trait_SS3$trait == "insect_predator")], decreasing = T)

effect_SS3_30_120_480 <- effect_SS3_30_120_480[ab_SS3_pred]
effect_SS3_30_120_480_lower <- effect_SS3_30_120_480_lower[ab_SS3_pred]
effect_SS3_30_120_480_upper <- effect_SS3_30_120_480_upper[ab_SS3_pred]











###################################################### SS4 #################################################################


fit_predators_SS4_treatments<- gllvm(com_predators_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                       treatments = relevel(treatments_SS2_3_4, ref = "control"),
                                                                       isolation = isolation_SS2_3_4),
                                     formula = ~ treatments + isolation,
                                     family = "negative.binomial",
                                     method = "VA",
                                     row.eff = ~ (1|ID),
                                     n.init = 10, num.lv = 0, seed = 11:20)



#effect_SS4 of TREATMENTS
effect_SS4_control_sugar_cane <-fit_predators_SS4_treatments$params$Xcoef[,2]
effect_SS4_control_sugar_cane_se <-fit_predators_SS4_treatments$sd$Xcoef[,2]
#effect_SS4_control_sugar_cane_se[which(effect_SS4_control_sugar_cane_se>100)] <- fit_predators_SS4_treatments_control$sd$Xcoef[,2][which(effect_SS4_control_sugar_cane_se>100)]

#effect_SS4_control_sugar_cane[which(effect_SS4_control_sugar_cane_se>100)] <- 0
effect_SS4_control_sugar_cane_se[which(effect_SS4_control_sugar_cane_se>100)] <- 0

effect_SS4_control_sugar_cane_upper <- effect_SS4_control_sugar_cane + (z*effect_SS4_control_sugar_cane_se)
effect_SS4_control_sugar_cane_lower <- effect_SS4_control_sugar_cane - (z*effect_SS4_control_sugar_cane_se)


ab_SS4_pred <- order(Trait_SS4$total_ab[which(Trait_SS4$trait == "insect_predator")], decreasing = T)
effect_SS4_control_sugar_cane <- effect_SS4_control_sugar_cane[ab_SS4_pred]
effect_SS4_control_sugar_cane_lower <- effect_SS4_control_sugar_cane_lower[ab_SS4_pred]
effect_SS4_control_sugar_cane_upper <- effect_SS4_control_sugar_cane_upper[ab_SS4_pred]




#effect_SS4 of  isolation


fit_predators_SS4_isolation_30 <- gllvm(com_predators_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                          isolation = relevel(isolation_SS2_3_4, ref = "30"),
                                                                          treatments = treatments_SS2_3_4),
                                        formula = ~ isolation + treatments,
                                        family = "negative.binomial",
                                        method = "VA",
                                        row.eff = ~ (1|ID),
                                        n.init = 10, num.lv = 0, seed = 11:20)




effect_SS4_30_120 <-fit_predators_SS4_isolation_30$params$Xcoef[,1]
effect_SS4_30_120_se <-fit_predators_SS4_isolation_30$sd$Xcoef[,1]
#effect_SS4_30_120_se[which(effect_SS4_30_120_se>100)] <- fit_predators_SS4_isolation_30$sd$Xcoef[,1][which(effect_SS4_30_120_se>100)]

#effect_SS4_30_120[which(effect_SS4_30_120_se>100)] <- 0
effect_SS4_30_120_se[which(effect_SS4_30_120_se>100)] <- 0


effect_SS4_30_480 <-fit_predators_SS4_isolation_30$params$Xcoef[,2]
effect_SS4_30_480_se <-fit_predators_SS4_isolation_30$sd$Xcoef[,2]
#effect_SS4_30_480_se[which(effect_SS4_30_480_se>100)] <- fit_predators_SS4_isolation_30$sd$Xcoef[,2][which(effect_SS4_30_480_se>100)]

#effect_SS4_30_480[which(effect_SS4_30_480_se>100)] <- 0
effect_SS4_30_480_se[which(effect_SS4_30_480_se>100)] <- 0

fit_predators_SS4_isolation_120 <- gllvm(com_predators_SS4, X = data.frame(ID = ID_SS2_3_4,
                                                                           isolation = relevel( isolation_SS2_3_4, ref = "120"),
                                                                           treatments = treatments_SS2_3_4),
                                         formula = ~  isolation + treatments,
                                         family = "negative.binomial",
                                         method = "VA",
                                         row.eff = ~ (1|ID),
                                         n.init = 10, num.lv = 0, seed = 11:20)


effect_SS4_120_480 <-fit_predators_SS4_isolation_120$params$Xcoef[,2]
effect_SS4_120_480_se <-fit_predators_SS4_isolation_120$sd$Xcoef[,2]
#effect_SS4_120_480_se[which(effect_SS4_120_480_se>100)] <- fit_predators_SS4_isolation_120$sd$Xcoef[,1][which(effect_SS4_120_480_se>100)]

#effect_SS4_120_480[which(effect_SS4_120_480_se>100)] <- 0
effect_SS4_120_480_se[which(effect_SS4_120_480_se>100)] <- 0

effect_SS4_30_120_lower <- effect_SS4_30_120 - (z*effect_SS4_30_120_se)
effect_SS4_30_120_upper <- effect_SS4_30_120 + (z*effect_SS4_30_120_se)

effect_SS4_30_480_lower <- effect_SS4_30_480 - (z*effect_SS4_30_480_se)
effect_SS4_30_480_upper <- effect_SS4_30_480 + (z*effect_SS4_30_480_se)

effect_SS4_120_480_lower <- effect_SS4_120_480 - (z*effect_SS4_120_480_se)
effect_SS4_120_480_upper <- effect_SS4_120_480 + (z*effect_SS4_120_480_se)



ab_SS4_pred <- order(Trait_SS4$total_ab[which(Trait_SS4$trait == "insect_predator")], decreasing = T)

effect_SS4_30_120 <- effect_SS4_30_120[ab_SS4_pred]
effect_SS4_30_120_lower <- effect_SS4_30_120_lower[ab_SS4_pred]
effect_SS4_30_120_upper <- effect_SS4_30_120_upper[ab_SS4_pred]

effect_SS4_30_480 <- effect_SS4_30_480[ab_SS4_pred]
effect_SS4_30_480_lower <- effect_SS4_30_480_lower[ab_SS4_pred]
effect_SS4_30_480_upper <- effect_SS4_30_480_upper[ab_SS4_pred]

effect_SS4_120_480 <- effect_SS4_120_480[ab_SS4_pred]
effect_SS4_120_480_lower <- effect_SS4_120_480_lower[ab_SS4_pred]
effect_SS4_120_480_upper <- effect_SS4_120_480_upper[ab_SS4_pred]









