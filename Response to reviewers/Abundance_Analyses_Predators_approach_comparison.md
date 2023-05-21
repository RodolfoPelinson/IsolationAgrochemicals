Abundance Analyses - Predators
================
Rodolfo Pelinson
26/12/2022

``` r
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(AICcmodavg)
```

Before anything, lets load and prepare some data sheets and vectors for
the analysis.

# Predators

First, lets make vectors of the abundances (summing up the abundance of
all predators taxa in each sample):

``` r
sum_com_predators_ab_SS1 <- rowSums(sum_com_orig_SS1[,which(Trait_SS1_sum_orig$trait == "insect_predator")])
sum_com_predators_ab_SS2 <- rowSums(sum_com_orig_SS2[,which(Trait_SS2_sum_orig$trait == "insect_predator")])
sum_com_predators_ab_SS3 <- rowSums(sum_com_orig_SS3[,which(Trait_SS3_sum_orig$trait == "insect_predator")])
sum_com_predators_ab_SS4 <- rowSums(sum_com_orig_SS4[,which(Trait_SS4_sum_orig$trait == "insect_predator")])

com_predators_ab_SS1 <- rowSums(com_SS1_orig[,which(Trait_SS1_orig$trait == "insect_predator")])
com_predators_ab_SS2 <- rowSums(com_SS2_orig[,which(Trait_SS2_orig$trait == "insect_predator")])
com_predators_ab_SS3 <- rowSums(com_SS3_orig[,which(Trait_SS3_orig$trait == "insect_predator")])
com_predators_ab_SS4 <- rowSums(com_SS4_orig[,which(Trait_SS4_orig$trait == "insect_predator")])
```

# First Survey (20 day)

Model Selection for mixed models

``` r
data_SS1 <- data.frame(isolation = isolation_SS1, treatments = treatments_SS1)

mod_pred_SS1_no_effect <- glmmTMB(com_predators_ab_SS1 ~ 1, family = nbinom2(link = "log"), data = data_SS1)
mod_pred_SS1_treatments <- glmmTMB(com_predators_ab_SS1 ~ treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_pred_SS1_isolation <- glmmTMB(com_predators_ab_SS1 ~ isolation, family = nbinom2(link = "log"), data = data_SS1)
mod_pred_SS1_isolation_treatments <- glmmTMB(com_predators_ab_SS1 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_pred_SS1_isolation_i_treatments <- glmmTMB(com_predators_ab_SS1 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS1)

model_selection_predators_SS1 <- aictab(list(mod_pred_SS1_no_effect,
                                             mod_pred_SS1_treatments,
                                             mod_pred_SS1_isolation,
                                             mod_pred_SS1_isolation_treatments,
                                             mod_pred_SS1_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_predators_SS1
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 231.59       0.00   0.58 -113.65
    ## Treatments              4 232.71       1.13   0.33 -111.86
    ## Isolation               4 235.99       4.40   0.06 -113.50
    ## Treatments + Isolation  6 237.50       5.91   0.03 -111.64
    ## Treatments * Isolation 10 248.90      17.32   0.00 -111.22

No effects.

``` r
mod_predators_SS1_treatments <- glmmTMB(com_predators_ab_SS1 ~ treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_predators_SS1_isolation <- glmmTMB(com_predators_ab_SS1 ~ isolation, family = nbinom2(link = "log"), data = data_SS1)
mod_predators_SS1_isolation_treatments <- glmmTMB(com_predators_ab_SS1 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_predators_SS1_isolation_i_treatments <- glmmTMB(com_predators_ab_SS1 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS1)

#Isolation
anova(mod_predators_SS1_treatments, mod_predators_SS1_isolation_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_predators_SS1_treatments: com_predators_ab_SS1 ~ treatments, zi=~0, disp=~1
    ## mod_predators_SS1_isolation_treatments: com_predators_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS1_treatments            4 231.72 238.94 -111.86   223.72       
    ## mod_predators_SS1_isolation_treatments  6 235.29 246.13 -111.64   223.29 0.4261
    ##                                        Chi Df Pr(>Chisq)
    ## mod_predators_SS1_treatments                            
    ## mod_predators_SS1_isolation_treatments      2     0.8081

``` r
#Treatments
anova(mod_predators_SS1_isolation, mod_predators_SS1_isolation_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_predators_SS1_isolation: com_predators_ab_SS1 ~ isolation, zi=~0, disp=~1
    ## mod_predators_SS1_isolation_treatments: com_predators_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS1_isolation             4 234.99 242.22 -113.50   226.99       
    ## mod_predators_SS1_isolation_treatments  6 235.29 246.13 -111.64   223.29 3.7015
    ##                                        Chi Df Pr(>Chisq)
    ## mod_predators_SS1_isolation                             
    ## mod_predators_SS1_isolation_treatments      2     0.1571

``` r
#Interaction
anova(mod_predators_SS1_isolation_treatments, mod_predators_SS1_isolation_i_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_predators_SS1_isolation_treatments: com_predators_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_predators_SS1_isolation_i_treatments: com_predators_ab_SS1 ~ isolation * treatments, zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_predators_SS1_isolation_treatments    6 235.29 246.13 -111.64   223.29
    ## mod_predators_SS1_isolation_i_treatments 10 242.43 260.50 -111.22   222.43
    ##                                           Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS1_isolation_treatments                           
    ## mod_predators_SS1_isolation_i_treatments 0.8564      4     0.9307

No effects.

# Second Survey (40 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS2 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_pred_SS2_no_effect <- glmmTMB(com_predators_ab_SS2 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_pred_SS2_treatments <- glmmTMB(com_predators_ab_SS2 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_pred_SS2_isolation <- glmmTMB(com_predators_ab_SS2 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_pred_SS2_isolation_treatments <- glmmTMB(com_predators_ab_SS2 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_pred_SS2_isolation_i_treatments <- glmmTMB(com_predators_ab_SS2 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)

model_selection_predators_SS2 <- aictab(list(mod_pred_SS2_no_effect,
                                             mod_pred_SS2_treatments,
                                             mod_pred_SS2_isolation,
                                             mod_pred_SS2_isolation_treatments,
                                             mod_pred_SS2_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_predators_SS2
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               3 854.57      47.93   0.00 -424.22
    ## Treatments              5 806.63       0.00   0.79 -398.14
    ## Isolation               5 858.30      51.67   0.00 -423.98
    ## Treatments + Isolation  7 809.42       2.78   0.20 -397.38
    ## Treatments * Isolation 11 815.16       8.53   0.01 -395.80

Effect of treatments.

LRT for mixed models

``` r
mod_predators_SS2_treatments <- glmmTMB(com_predators_ab_SS2 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_predators_SS2_isolation <- glmmTMB(com_predators_ab_SS2 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_predators_SS2_isolation_treatments <- glmmTMB(com_predators_ab_SS2 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_predators_SS2_isolation_i_treatments <- glmmTMB(com_predators_ab_SS2 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)

#Isolation
anova(mod_predators_SS2_treatments, mod_predators_SS2_isolation_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_predators_SS2_treatments: com_predators_ab_SS2 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS2_isolation_treatments: com_predators_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS2_treatments            5 806.29 822.25 -398.14   796.29       
    ## mod_predators_SS2_isolation_treatments  7 808.77 831.12 -397.38   794.77 1.5235
    ##                                        Chi Df Pr(>Chisq)
    ## mod_predators_SS2_treatments                            
    ## mod_predators_SS2_isolation_treatments      2     0.4669

``` r
#Treatments
anova(mod_predators_SS2_isolation, mod_predators_SS2_isolation_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_predators_SS2_isolation: com_predators_ab_SS2 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS2_isolation_treatments: com_predators_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS2_isolation             5 857.96 873.92 -423.98   847.96       
    ## mod_predators_SS2_isolation_treatments  7 808.77 831.12 -397.38   794.77 53.191
    ##                                        Chi Df Pr(>Chisq)    
    ## mod_predators_SS2_isolation                                 
    ## mod_predators_SS2_isolation_treatments      2  2.816e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_predators_SS2_isolation_treatments, mod_predators_SS2_isolation_i_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_predators_SS2_isolation_treatments: com_predators_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS2_isolation_i_treatments: com_predators_ab_SS2 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_predators_SS2_isolation_treatments    7 808.77 831.12 -397.38   794.77
    ## mod_predators_SS2_isolation_i_treatments 11 813.59 848.71 -395.80   791.59
    ##                                           Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS2_isolation_treatments                           
    ## mod_predators_SS2_isolation_i_treatments 3.1748      4      0.529

Effect of treatments.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS2_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_pred_SS2_no_effect_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ 1 , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_pred_SS2_treatments_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ treatments , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_pred_SS2_isolation_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ isolation , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_pred_SS2_isolation_treatments_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_pred_SS2_isolation_i_treatments_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS2_sum)

model_selection_predators_SS2_sum <- aictab(list(mod_pred_SS2_no_effect_sum,
                                             mod_pred_SS2_treatments_sum,
                                             mod_pred_SS2_isolation_sum,
                                             mod_pred_SS2_isolation_treatments_sum,
                                             mod_pred_SS2_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_predators_SS2_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 351.92      36.53   0.00 -173.82
    ## Treatments              4 315.40       0.00   0.79 -153.20
    ## Isolation               4 354.83      39.43   0.00 -172.91
    ## Treatments + Isolation  6 318.07       2.67   0.21 -151.93
    ## Treatments * Isolation 10 326.40      11.00   0.00 -149.97

Effect of treatments.

LRT for mixed models

``` r
mod_predators_SS2_treatments_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ treatments, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_predators_SS2_isolation_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ isolation, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_predators_SS2_isolation_treatments_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_predators_SS2_isolation_i_treatments_sum <- glmmTMB(sum_com_predators_ab_SS2 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS2_sum)

#Isolation
anova(mod_predators_SS2_treatments_sum, mod_predators_SS2_isolation_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_predators_SS2_treatments_sum: sum_com_predators_ab_SS2 ~ treatments, zi=~0, disp=~1
    ## mod_predators_SS2_isolation_treatments_sum: sum_com_predators_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_predators_SS2_treatments_sum            4 314.40 321.62 -153.20   306.40
    ## mod_predators_SS2_isolation_treatments_sum  6 315.86 326.69 -151.93   303.86
    ##                                             Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS2_treatments_sum                                   
    ## mod_predators_SS2_isolation_treatments_sum 2.5416      2     0.2806

``` r
#Treatments
anova(mod_predators_SS2_isolation_sum, mod_predators_SS2_isolation_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_predators_SS2_isolation_sum: sum_com_predators_ab_SS2 ~ isolation, zi=~0, disp=~1
    ## mod_predators_SS2_isolation_treatments_sum: sum_com_predators_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_predators_SS2_isolation_sum             4 353.83 361.05 -172.91   345.83
    ## mod_predators_SS2_isolation_treatments_sum  6 315.85 326.69 -151.93   303.85
    ##                                             Chisq Chi Df Pr(>Chisq)    
    ## mod_predators_SS2_isolation_sum                                        
    ## mod_predators_SS2_isolation_treatments_sum 41.971      2  7.693e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_predators_SS2_isolation_treatments_sum, mod_predators_SS2_isolation_i_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_predators_SS2_isolation_treatments_sum: sum_com_predators_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_predators_SS2_isolation_i_treatments_sum: sum_com_predators_ab_SS2 ~ isolation * treatments, zi=~0, disp=~1
    ##                                              Df    AIC    BIC  logLik deviance
    ## mod_predators_SS2_isolation_treatments_sum    6 315.85 326.69 -151.93   303.86
    ## mod_predators_SS2_isolation_i_treatments_sum 10 319.93 338.00 -149.97   299.93
    ##                                              Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS2_isolation_treatments_sum                          
    ## mod_predators_SS2_isolation_i_treatments_sum 3.924      4     0.4164

Effect of treatments.

# Third Survey (80 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS3 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_pred_SS3_no_effect <- glmmTMB(com_predators_ab_SS3 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_pred_SS3_treatments <- glmmTMB(com_predators_ab_SS3 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_pred_SS3_isolation <- glmmTMB(com_predators_ab_SS3 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_pred_SS3_isolation_treatments <- glmmTMB(com_predators_ab_SS3 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_pred_SS3_isolation_i_treatments <- glmmTMB(com_predators_ab_SS3 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)

model_selection_predators_SS3 <- aictab(list(mod_pred_SS3_no_effect,
                                             mod_pred_SS3_treatments,
                                             mod_pred_SS3_isolation,
                                             mod_pred_SS3_isolation_treatments,
                                             mod_pred_SS3_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_predators_SS3
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K    AICc Delta_AICc AICcWt      LL
    ## No Effect               3 1277.13       5.22   0.05 -635.50
    ## Treatments              5 1271.91       0.00   0.73 -630.78
    ## Isolation               5 1280.23       8.32   0.01 -634.94
    ## Treatments + Isolation  7 1274.90       2.99   0.16 -630.13
    ## Treatments * Isolation 11 1277.52       5.61   0.04 -626.98

Effect of agrochemical treatments.

LRT for mixed models

``` r
mod_predators_SS3_treatments <- glmmTMB(com_predators_ab_SS3 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_predators_SS3_isolation <- glmmTMB(com_predators_ab_SS3 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_predators_SS3_isolation_treatments <- glmmTMB(com_predators_ab_SS3 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_predators_SS3_isolation_i_treatments <- glmmTMB(com_predators_ab_SS3 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)

#Isolation
anova(mod_predators_SS3_treatments, mod_predators_SS3_isolation_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_predators_SS3_treatments: com_predators_ab_SS3 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS3_isolation_treatments: com_predators_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS3_treatments            5 1271.6 1287.5 -630.78   1261.6       
    ## mod_predators_SS3_isolation_treatments  7 1274.2 1296.6 -630.13   1260.2 1.3131
    ##                                        Chi Df Pr(>Chisq)
    ## mod_predators_SS3_treatments                            
    ## mod_predators_SS3_isolation_treatments      2     0.5186

``` r
#Treatments
anova(mod_predators_SS3_isolation, mod_predators_SS3_isolation_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_predators_SS3_isolation: com_predators_ab_SS3 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS3_isolation_treatments: com_predators_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS3_isolation             5 1279.9 1295.8 -634.94   1269.9       
    ## mod_predators_SS3_isolation_treatments  7 1274.2 1296.6 -630.13   1260.2 9.6354
    ##                                        Chi Df Pr(>Chisq)   
    ## mod_predators_SS3_isolation                                
    ## mod_predators_SS3_isolation_treatments      2   0.008086 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_predators_SS3_isolation_treatments, mod_predators_SS3_isolation_i_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_predators_SS3_isolation_treatments: com_predators_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS3_isolation_i_treatments: com_predators_ab_SS3 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_predators_SS3_isolation_treatments    7 1274.2 1296.6 -630.13   1260.2
    ## mod_predators_SS3_isolation_i_treatments 11 1276.0 1311.1 -626.98   1254.0
    ##                                           Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS3_isolation_treatments                           
    ## mod_predators_SS3_isolation_i_treatments 6.3013      4     0.1777

Effect of agrochemical treatments.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS3_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_pred_SS3_no_effect_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ 1 , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_pred_SS3_treatments_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ treatments , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_pred_SS3_isolation_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ isolation , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_pred_SS3_isolation_treatments_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_pred_SS3_isolation_i_treatments_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS3_sum)

model_selection_predators_SS3_sum <- aictab(list(mod_pred_SS3_no_effect_sum,
                                             mod_pred_SS3_treatments_sum,
                                             mod_pred_SS3_isolation_sum,
                                             mod_pred_SS3_isolation_treatments_sum,
                                             mod_pred_SS3_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_predators_SS3_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 433.95       3.41   0.13 -214.83
    ## Treatments              4 430.54       0.00   0.74 -210.77
    ## Isolation               4 437.32       6.78   0.02 -214.16
    ## Treatments + Isolation  6 434.57       4.04   0.10 -210.18
    ## Treatments * Isolation 10 439.93       9.39   0.01 -206.73

Effect of agrochemical treatments.

LRT for mixed models

``` r
mod_predators_SS3_treatments_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ treatments, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_predators_SS3_isolation_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ isolation, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_predators_SS3_isolation_treatments_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_predators_SS3_isolation_i_treatments_sum <- glmmTMB(sum_com_predators_ab_SS3 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS3_sum)

#Isolation
anova(mod_predators_SS3_treatments_sum, mod_predators_SS3_isolation_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_predators_SS3_treatments_sum: sum_com_predators_ab_SS3 ~ treatments, zi=~0, disp=~1
    ## mod_predators_SS3_isolation_treatments_sum: sum_com_predators_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_predators_SS3_treatments_sum            4 429.54 436.76 -210.77   421.54
    ## mod_predators_SS3_isolation_treatments_sum  6 432.36 443.20 -210.18   420.36
    ##                                             Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS3_treatments_sum                                   
    ## mod_predators_SS3_isolation_treatments_sum 1.1739      2      0.556

``` r
#Treatments
anova(mod_predators_SS3_isolation_sum, mod_predators_SS3_isolation_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_predators_SS3_isolation_sum: sum_com_predators_ab_SS3 ~ isolation, zi=~0, disp=~1
    ## mod_predators_SS3_isolation_treatments_sum: sum_com_predators_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_predators_SS3_isolation_sum             4 436.32 443.55 -214.16   428.32
    ## mod_predators_SS3_isolation_treatments_sum  6 432.36 443.20 -210.18   420.36
    ##                                             Chisq Chi Df Pr(>Chisq)  
    ## mod_predators_SS3_isolation_sum                                      
    ## mod_predators_SS3_isolation_treatments_sum 7.9564      2    0.01872 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_predators_SS3_isolation_treatments_sum, mod_predators_SS3_isolation_i_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_predators_SS3_isolation_treatments_sum: sum_com_predators_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_predators_SS3_isolation_i_treatments_sum: sum_com_predators_ab_SS3 ~ isolation * treatments, zi=~0, disp=~1
    ##                                              Df    AIC    BIC  logLik deviance
    ## mod_predators_SS3_isolation_treatments_sum    6 432.36 443.20 -210.18   420.36
    ## mod_predators_SS3_isolation_i_treatments_sum 10 433.46 451.52 -206.73   413.46
    ##                                               Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS3_isolation_treatments_sum                           
    ## mod_predators_SS3_isolation_i_treatments_sum 6.9055      4      0.141

Effect of agrochemical treatments.

# Fourth Survey (160 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS4 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_pred_SS4_no_effect <- glmmTMB(com_predators_ab_SS4 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_pred_SS4_treatments <- glmmTMB(com_predators_ab_SS4 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_pred_SS4_isolation <- glmmTMB(com_predators_ab_SS4 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_pred_SS4_isolation_treatments <- glmmTMB(com_predators_ab_SS4 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_pred_SS4_isolation_i_treatments <- glmmTMB(com_predators_ab_SS4 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)

model_selection_predators_SS4 <- aictab(list(mod_pred_SS4_no_effect,
                                             mod_pred_SS4_treatments,
                                             mod_pred_SS4_isolation,
                                             mod_pred_SS4_isolation_treatments,
                                             mod_pred_SS4_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_predators_SS4
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K    AICc Delta_AICc AICcWt      LL
    ## No Effect               3 1144.11       0.00   0.72 -568.98
    ## Treatments              5 1147.29       3.19   0.15 -568.47
    ## Isolation               5 1148.06       3.96   0.10 -568.86
    ## Treatments + Isolation  7 1151.34       7.23   0.02 -568.34
    ## Treatments * Isolation 11 1151.51       7.41   0.02 -563.97

No effects.

LRT for mixed models

``` r
mod_predators_SS4_treatments <- glmmTMB(com_predators_ab_SS4 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_predators_SS4_isolation <- glmmTMB(com_predators_ab_SS4 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_predators_SS4_isolation_treatments <- glmmTMB(com_predators_ab_SS4 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_predators_SS4_isolation_i_treatments <- glmmTMB(com_predators_ab_SS4 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)

#Isolation
anova(mod_predators_SS4_treatments, mod_predators_SS4_isolation_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_predators_SS4_treatments: com_predators_ab_SS4 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS4_isolation_treatments: com_predators_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS4_treatments            5 1147.0 1162.9 -568.47   1137.0       
    ## mod_predators_SS4_isolation_treatments  7 1150.7 1173.0 -568.34   1136.7 0.2626
    ##                                        Chi Df Pr(>Chisq)
    ## mod_predators_SS4_treatments                            
    ## mod_predators_SS4_isolation_treatments      2      0.877

``` r
#Treatments
anova(mod_predators_SS4_isolation, mod_predators_SS4_isolation_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_predators_SS4_isolation: com_predators_ab_SS4 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS4_isolation_treatments: com_predators_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_predators_SS4_isolation             5 1147.7 1163.7 -568.86   1137.7       
    ## mod_predators_SS4_isolation_treatments  7 1150.7 1173.0 -568.34   1136.7 1.0316
    ##                                        Chi Df Pr(>Chisq)
    ## mod_predators_SS4_isolation                             
    ## mod_predators_SS4_isolation_treatments      2      0.597

``` r
#Interaction
anova(mod_predators_SS4_isolation_treatments, mod_predators_SS4_isolation_i_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_predators_SS4_isolation_treatments: com_predators_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_predators_SS4_isolation_i_treatments: com_predators_ab_SS4 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_predators_SS4_isolation_treatments    7 1150.7 1173.0 -568.34   1136.7
    ## mod_predators_SS4_isolation_i_treatments 11 1149.9 1185.1 -563.97   1127.9
    ##                                           Chisq Chi Df Pr(>Chisq)  
    ## mod_predators_SS4_isolation_treatments                             
    ## mod_predators_SS4_isolation_i_treatments 8.7462      4    0.06777 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No effects. Actually, marginal effects of interaction.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS4_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_pred_SS4_no_effect_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ 1 , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_pred_SS4_treatments_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ treatments , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_pred_SS4_isolation_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ isolation , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_pred_SS4_isolation_treatments_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_pred_SS4_isolation_i_treatments_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS4_sum)

model_selection_predators_SS4_sum <- aictab(list(mod_pred_SS4_no_effect_sum,
                                             mod_pred_SS4_treatments_sum,
                                             mod_pred_SS4_isolation_sum,
                                             mod_pred_SS4_isolation_treatments_sum,
                                             mod_pred_SS4_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_predators_SS4_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 374.06       0.00   0.78 -184.89
    ## Treatments              4 377.78       3.73   0.12 -184.39
    ## Isolation               4 378.35       4.29   0.09 -184.68
    ## Treatments + Isolation  6 382.55       8.49   0.01 -184.17
    ## Treatments * Isolation 10 385.70      11.64   0.00 -179.62

No effects.

LRT for mixed models

``` r
mod_predators_SS4_treatments_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ treatments, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_predators_SS4_isolation_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ isolation, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_predators_SS4_isolation_treatments_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_predators_SS4_isolation_i_treatments_sum <- glmmTMB(sum_com_predators_ab_SS4 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS4_sum)

#Isolation
anova(mod_predators_SS4_treatments_sum, mod_predators_SS4_isolation_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_predators_SS4_treatments_sum: sum_com_predators_ab_SS4 ~ treatments, zi=~0, disp=~1
    ## mod_predators_SS4_isolation_treatments_sum: sum_com_predators_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_predators_SS4_treatments_sum            4 376.78 384.01 -184.39   368.78
    ## mod_predators_SS4_isolation_treatments_sum  6 380.34 391.18 -184.17   368.34
    ##                                             Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS4_treatments_sum                                   
    ## mod_predators_SS4_isolation_treatments_sum 0.4449      2     0.8006

``` r
#Treatments
anova(mod_predators_SS4_isolation_sum, mod_predators_SS4_isolation_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_predators_SS4_isolation_sum: sum_com_predators_ab_SS4 ~ isolation, zi=~0, disp=~1
    ## mod_predators_SS4_isolation_treatments_sum: sum_com_predators_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_predators_SS4_isolation_sum             4 377.35 384.58 -184.68   369.35
    ## mod_predators_SS4_isolation_treatments_sum  6 380.34 391.18 -184.17   368.34
    ##                                             Chisq Chi Df Pr(>Chisq)
    ## mod_predators_SS4_isolation_sum                                    
    ## mod_predators_SS4_isolation_treatments_sum 1.0128      2     0.6027

``` r
#Interaction
anova(mod_predators_SS4_isolation_treatments_sum, mod_predators_SS4_isolation_i_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_predators_SS4_isolation_treatments_sum: sum_com_predators_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_predators_SS4_isolation_i_treatments_sum: sum_com_predators_ab_SS4 ~ isolation * treatments, zi=~0, disp=~1
    ##                                              Df    AIC    BIC  logLik deviance
    ## mod_predators_SS4_isolation_treatments_sum    6 380.34 391.18 -184.17   368.34
    ## mod_predators_SS4_isolation_i_treatments_sum 10 379.23 397.30 -179.62   359.23
    ##                                               Chisq Chi Df Pr(>Chisq)  
    ## mod_predators_SS4_isolation_treatments_sum                             
    ## mod_predators_SS4_isolation_i_treatments_sum 9.1081      4    0.05845 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No effects. Actually, marginal effects of interaction.
