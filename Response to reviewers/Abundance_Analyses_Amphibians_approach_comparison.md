Abundance Analyses - Herbivore and Detritivores - Amphibians
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

# Amphibians

First, lets make vectors of the abundances (summing up the abundance of
all amphibian taxa in each sample):

``` r
sum_com_amphibian_ab_SS1 <- rowSums(sum_com_orig_SS1[,which(Trait_SS1_sum_orig$trait == "amphibian_consumer")])
sum_com_amphibian_ab_SS2 <- rowSums(sum_com_orig_SS2[,which(Trait_SS2_sum_orig$trait == "amphibian_consumer")])
sum_com_amphibian_ab_SS3 <- rowSums(sum_com_orig_SS3[,which(Trait_SS3_sum_orig$trait == "amphibian_consumer")])
sum_com_amphibian_ab_SS4 <- rowSums(sum_com_orig_SS4[,which(Trait_SS4_sum_orig$trait == "amphibian_consumer")])

com_amphibian_ab_SS1 <- rowSums(com_SS1_orig[,which(Trait_SS1_orig$trait == "amphibian_consumer")])
com_amphibian_ab_SS2 <- rowSums(com_SS2_orig[,which(Trait_SS2_orig$trait == "amphibian_consumer")])
com_amphibian_ab_SS3 <- rowSums(com_SS3_orig[,which(Trait_SS3_orig$trait == "amphibian_consumer")])
com_amphibian_ab_SS4 <- rowSums(com_SS4_orig[,which(Trait_SS4_orig$trait == "amphibian_consumer")])
```

# First Survey (20 day)

Model Selection for mixed models

``` r
data_SS1 <- data.frame(isolation = isolation_SS1, treatments = treatments_SS1)

mod_amph_SS1_no_effect <- glmmTMB(com_amphibian_ab_SS1 ~ 1, family = nbinom2(link = "log"), data = data_SS1)
mod_amph_SS1_treatments <- glmmTMB(com_amphibian_ab_SS1 ~ treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_amph_SS1_isolation <- glmmTMB(com_amphibian_ab_SS1 ~ isolation, family = nbinom2(link = "log"), data = data_SS1)
mod_amph_SS1_isolation_treatments <- glmmTMB(com_amphibian_ab_SS1 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_amph_SS1_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS1 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS1)

model_selection_amphibian_SS1 <- aictab(list(mod_amph_SS1_no_effect,
                                             mod_amph_SS1_treatments,
                                             mod_amph_SS1_isolation,
                                             mod_amph_SS1_isolation_treatments,
                                             mod_amph_SS1_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_amphibian_SS1
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 255.25       1.77   0.26 -125.48
    ## Treatments              4 258.43       4.94   0.05 -124.71
    ## Isolation               4 253.49       0.00   0.62 -122.24
    ## Treatments + Isolation  6 258.07       4.58   0.06 -121.93
    ## Treatments * Isolation 10 265.74      12.26   0.00 -119.64

Marginal effect of isolation

``` r
mod_amphibian_SS1_treatments <- glmmTMB(com_amphibian_ab_SS1 ~ treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_amphibian_SS1_isolation <- glmmTMB(com_amphibian_ab_SS1 ~ isolation, family = nbinom2(link = "log"), data = data_SS1)
mod_amphibian_SS1_isolation_treatments <- glmmTMB(com_amphibian_ab_SS1 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_amphibian_SS1_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS1 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS1)

#Isolation
anova(mod_amphibian_SS1_treatments, mod_amphibian_SS1_isolation_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_amphibian_SS1_treatments: com_amphibian_ab_SS1 ~ treatments, zi=~0, disp=~1
    ## mod_amphibian_SS1_isolation_treatments: com_amphibian_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS1_treatments            4 257.43 264.66 -124.72   249.43       
    ## mod_amphibian_SS1_isolation_treatments  6 255.86 266.70 -121.93   243.86 5.5741
    ##                                        Chi Df Pr(>Chisq)  
    ## mod_amphibian_SS1_treatments                              
    ## mod_amphibian_SS1_isolation_treatments      2     0.0616 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_amphibian_SS1_isolation, mod_amphibian_SS1_isolation_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_amphibian_SS1_isolation: com_amphibian_ab_SS1 ~ isolation, zi=~0, disp=~1
    ## mod_amphibian_SS1_isolation_treatments: com_amphibian_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS1_isolation             4 252.49 259.71 -122.24   244.49       
    ## mod_amphibian_SS1_isolation_treatments  6 255.86 266.70 -121.93   243.86 0.6298
    ##                                        Chi Df Pr(>Chisq)
    ## mod_amphibian_SS1_isolation                             
    ## mod_amphibian_SS1_isolation_treatments      2     0.7298

``` r
#Interaction
anova(mod_amphibian_SS1_isolation_treatments, mod_amphibian_SS1_isolation_i_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_amphibian_SS1_isolation_treatments: com_amphibian_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_amphibian_SS1_isolation_i_treatments: com_amphibian_ab_SS1 ~ isolation * treatments, zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS1_isolation_treatments    6 255.86 266.70 -121.93   243.86
    ## mod_amphibian_SS1_isolation_i_treatments 10 259.27 277.34 -119.64   239.27
    ##                                           Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS1_isolation_treatments                           
    ## mod_amphibian_SS1_isolation_i_treatments 4.5824      4     0.3329

Marginal effect of isolation

# Second Survey (40 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS2 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_amph_SS2_no_effect <- glmmTMB(com_amphibian_ab_SS2 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_amph_SS2_treatments <- glmmTMB(com_amphibian_ab_SS2 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_amph_SS2_isolation <- glmmTMB(com_amphibian_ab_SS2 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_amph_SS2_isolation_treatments <- glmmTMB(com_amphibian_ab_SS2 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_amph_SS2_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS2 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)

model_selection_amphibian_SS2 <- aictab(list(mod_amph_SS2_no_effect,
                                             mod_amph_SS2_treatments,
                                             mod_amph_SS2_isolation,
                                             mod_amph_SS2_isolation_treatments,
                                             mod_amph_SS2_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_amphibian_SS2
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               3 896.05       7.24   0.02 -444.95
    ## Treatments              5 900.04      11.23   0.00 -444.85
    ## Isolation               5 888.81       0.00   0.86 -439.23
    ## Treatments + Isolation  7 892.93       4.12   0.11 -439.14
    ## Treatments * Isolation 11 899.75      10.95   0.00 -438.09

Effect of isolation

LRT for mixed models

``` r
mod_amphibian_SS2_treatments <- glmmTMB(com_amphibian_ab_SS2 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_amphibian_SS2_isolation <- glmmTMB(com_amphibian_ab_SS2 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_amphibian_SS2_isolation_treatments <- glmmTMB(com_amphibian_ab_SS2 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_amphibian_SS2_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS2 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)

#Isolation
anova(mod_amphibian_SS2_treatments, mod_amphibian_SS2_isolation_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_amphibian_SS2_treatments: com_amphibian_ab_SS2 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS2_isolation_treatments: com_amphibian_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS2_treatments            5 899.69 915.66 -444.85   889.69       
    ## mod_amphibian_SS2_isolation_treatments  7 892.28 914.63 -439.14   878.28 11.415
    ##                                        Chi Df Pr(>Chisq)   
    ## mod_amphibian_SS2_treatments                               
    ## mod_amphibian_SS2_isolation_treatments      2   0.003322 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_amphibian_SS2_isolation, mod_amphibian_SS2_isolation_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_amphibian_SS2_isolation: com_amphibian_ab_SS2 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS2_isolation_treatments: com_amphibian_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS2_isolation             5 888.46 904.43 -439.23   878.46       
    ## mod_amphibian_SS2_isolation_treatments  7 892.28 914.63 -439.14   878.28 0.1827
    ##                                        Chi Df Pr(>Chisq)
    ## mod_amphibian_SS2_isolation                             
    ## mod_amphibian_SS2_isolation_treatments      2     0.9127

``` r
#Interaction
anova(mod_amphibian_SS2_isolation_treatments, mod_amphibian_SS2_isolation_i_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_amphibian_SS2_isolation_treatments: com_amphibian_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS2_isolation_i_treatments: com_amphibian_ab_SS2 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS2_isolation_treatments    7 892.28 914.63 -439.14   878.28
    ## mod_amphibian_SS2_isolation_i_treatments 11 898.18 933.30 -438.09   876.18
    ##                                          Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS2_isolation_treatments                          
    ## mod_amphibian_SS2_isolation_i_treatments 2.097      4     0.7179

Effect of isolation.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS2_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_amph_SS2_no_effect_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ 1 , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_amph_SS2_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ treatments , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_amph_SS2_isolation_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ isolation , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_amph_SS2_isolation_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_amph_SS2_isolation_i_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS2_sum)

model_selection_amphibian_SS2_sum <- aictab(list(mod_amph_SS2_no_effect_sum,
                                             mod_amph_SS2_treatments_sum,
                                             mod_amph_SS2_isolation_sum,
                                             mod_amph_SS2_isolation_treatments_sum,
                                             mod_amph_SS2_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_amphibian_SS2_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 370.52       3.15   0.15 -183.12
    ## Treatments              4 373.33       5.96   0.04 -182.17
    ## Isolation               4 367.38       0.00   0.72 -179.19
    ## Treatments + Isolation  6 371.48       4.11   0.09 -178.64
    ## Treatments * Isolation 10 380.30      12.92   0.00 -176.91

Effect of isolation

LRT for mixed models

``` r
mod_amphibian_SS2_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ treatments, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_amphibian_SS2_isolation_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ isolation, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_amphibian_SS2_isolation_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_amphibian_SS2_isolation_i_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS2 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS2_sum)

#Isolation
anova(mod_amphibian_SS2_treatments_sum, mod_amphibian_SS2_isolation_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_amphibian_SS2_treatments_sum: sum_com_amphibian_ab_SS2 ~ treatments, zi=~0, disp=~1
    ## mod_amphibian_SS2_isolation_treatments_sum: sum_com_amphibian_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS2_treatments_sum            4 372.33 379.56 -182.17   364.33
    ## mod_amphibian_SS2_isolation_treatments_sum  6 369.27 380.11 -178.63   357.27
    ##                                             Chisq Chi Df Pr(>Chisq)  
    ## mod_amphibian_SS2_treatments_sum                                     
    ## mod_amphibian_SS2_isolation_treatments_sum 7.0609      2    0.02929 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_amphibian_SS2_isolation_sum, mod_amphibian_SS2_isolation_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_amphibian_SS2_isolation_sum: sum_com_amphibian_ab_SS2 ~ isolation, zi=~0, disp=~1
    ## mod_amphibian_SS2_isolation_treatments_sum: sum_com_amphibian_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS2_isolation_sum             4 366.38 373.60 -179.19   358.38
    ## mod_amphibian_SS2_isolation_treatments_sum  6 369.27 380.11 -178.63   357.27
    ##                                            Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS2_isolation_sum                                   
    ## mod_amphibian_SS2_isolation_treatments_sum 1.105      2     0.5755

``` r
#Interaction
anova(mod_amphibian_SS2_isolation_treatments_sum, mod_amphibian_SS2_isolation_i_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_amphibian_SS2_isolation_treatments_sum: sum_com_amphibian_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_amphibian_SS2_isolation_i_treatments_sum: sum_com_amphibian_ab_SS2 ~ isolation * treatments, zi=~0, disp=~1
    ##                                              Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS2_isolation_treatments_sum    6 369.27 380.11 -178.63   357.27
    ## mod_amphibian_SS2_isolation_i_treatments_sum 10 373.83 391.89 -176.91   353.83
    ##                                               Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS2_isolation_treatments_sum                           
    ## mod_amphibian_SS2_isolation_i_treatments_sum 3.4439      4     0.4865

Effect of isolation

# Third Survey (80 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS3 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_amph_SS3_no_effect <- glmmTMB(com_amphibian_ab_SS3 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_amph_SS3_treatments <- glmmTMB(com_amphibian_ab_SS3 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_amph_SS3_isolation <- glmmTMB(com_amphibian_ab_SS3 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_amph_SS3_isolation_treatments <- glmmTMB(com_amphibian_ab_SS3 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_amph_SS3_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS3 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)

model_selection_amphibian_SS3 <- aictab(list(mod_amph_SS3_no_effect,
                                             mod_amph_SS3_treatments,
                                             mod_amph_SS3_isolation,
                                             mod_amph_SS3_isolation_treatments,
                                             mod_amph_SS3_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_amphibian_SS3
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K    AICc Delta_AICc AICcWt      LL
    ## No Effect               3 1202.62       4.39   0.08 -598.24
    ## Treatments              5 1198.24       0.00   0.71 -593.95
    ## Isolation               5 1205.50       7.26   0.02 -597.58
    ## Treatments + Isolation  7 1200.90       2.66   0.19 -593.12
    ## Treatments * Isolation 11 1207.96       9.72   0.01 -592.19

Effect of agrochemical treatments.

LRT for mixed models

``` r
mod_amphibian_SS3_treatments <- glmmTMB(com_amphibian_ab_SS3 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_amphibian_SS3_isolation <- glmmTMB(com_amphibian_ab_SS3 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_amphibian_SS3_isolation_treatments <- glmmTMB(com_amphibian_ab_SS3 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_amphibian_SS3_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS3 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)

#Isolation
anova(mod_amphibian_SS3_treatments, mod_amphibian_SS3_isolation_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_amphibian_SS3_treatments: com_amphibian_ab_SS3 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS3_isolation_treatments: com_amphibian_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS3_treatments            5 1197.9 1213.9 -593.95   1187.9       
    ## mod_amphibian_SS3_isolation_treatments  7 1200.2 1222.6 -593.12   1186.2 1.6479
    ##                                        Chi Df Pr(>Chisq)
    ## mod_amphibian_SS3_treatments                            
    ## mod_amphibian_SS3_isolation_treatments      2     0.4387

``` r
#Treatments
anova(mod_amphibian_SS3_isolation, mod_amphibian_SS3_isolation_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_amphibian_SS3_isolation: com_amphibian_ab_SS3 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS3_isolation_treatments: com_amphibian_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS3_isolation             5 1205.2 1221.1 -597.58   1195.2       
    ## mod_amphibian_SS3_isolation_treatments  7 1200.2 1222.6 -593.12   1186.2 8.9051
    ##                                        Chi Df Pr(>Chisq)  
    ## mod_amphibian_SS3_isolation                               
    ## mod_amphibian_SS3_isolation_treatments      2    0.01165 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_amphibian_SS3_isolation_treatments, mod_amphibian_SS3_isolation_i_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_amphibian_SS3_isolation_treatments: com_amphibian_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS3_isolation_i_treatments: com_amphibian_ab_SS3 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS3_isolation_treatments    7 1200.2 1222.6 -593.12   1186.2
    ## mod_amphibian_SS3_isolation_i_treatments 11 1206.4 1241.5 -592.19   1184.4
    ##                                           Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS3_isolation_treatments                           
    ## mod_amphibian_SS3_isolation_i_treatments 1.8608      4     0.7613

Effect of agrochemical treatments.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS3_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_amph_SS3_no_effect_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ 1 , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_amph_SS3_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ treatments , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_amph_SS3_isolation_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ isolation , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_amph_SS3_isolation_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_amph_SS3_isolation_i_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS3_sum)

model_selection_amphibian_SS3_sum <- aictab(list(mod_amph_SS3_no_effect_sum,
                                             mod_amph_SS3_treatments_sum,
                                             mod_amph_SS3_isolation_sum,
                                             mod_amph_SS3_isolation_treatments_sum,
                                             mod_amph_SS3_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_amphibian_SS3_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 448.69       6.26   0.04 -222.20
    ## Treatments              4 442.43       0.00   0.87 -216.72
    ## Isolation               4 453.38      10.95   0.00 -222.19
    ## Treatments + Isolation  6 447.12       4.69   0.08 -216.45
    ## Treatments * Isolation 10 457.78      15.35   0.00 -215.66

Effect of agrochemical treatments.

LRT for mixed models

``` r
mod_amphibian_SS3_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ treatments, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_amphibian_SS3_isolation_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ isolation, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_amphibian_SS3_isolation_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_amphibian_SS3_isolation_i_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS3 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS3_sum)

#Isolation
anova(mod_amphibian_SS3_treatments_sum, mod_amphibian_SS3_isolation_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_amphibian_SS3_treatments_sum: sum_com_amphibian_ab_SS3 ~ treatments, zi=~0, disp=~1
    ## mod_amphibian_SS3_isolation_treatments_sum: sum_com_amphibian_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS3_treatments_sum            4 441.43 448.66 -216.72   433.43
    ## mod_amphibian_SS3_isolation_treatments_sum  6 444.91 455.75 -216.45   432.91
    ##                                             Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS3_treatments_sum                                   
    ## mod_amphibian_SS3_isolation_treatments_sum 0.5244      2     0.7694

``` r
#Treatments
anova(mod_amphibian_SS3_isolation_sum, mod_amphibian_SS3_isolation_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_amphibian_SS3_isolation_sum: sum_com_amphibian_ab_SS3 ~ isolation, zi=~0, disp=~1
    ## mod_amphibian_SS3_isolation_treatments_sum: sum_com_amphibian_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS3_isolation_sum             4 452.38 459.61 -222.19   444.38
    ## mod_amphibian_SS3_isolation_treatments_sum  6 444.91 455.75 -216.45   432.91
    ##                                             Chisq Chi Df Pr(>Chisq)   
    ## mod_amphibian_SS3_isolation_sum                                       
    ## mod_amphibian_SS3_isolation_treatments_sum 11.473      2   0.003227 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_amphibian_SS3_isolation_treatments_sum, mod_amphibian_SS3_isolation_i_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_amphibian_SS3_isolation_treatments_sum: sum_com_amphibian_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_amphibian_SS3_isolation_i_treatments_sum: sum_com_amphibian_ab_SS3 ~ isolation * treatments, zi=~0, disp=~1
    ##                                              Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS3_isolation_treatments_sum    6 444.91 455.75 -216.45   432.91
    ## mod_amphibian_SS3_isolation_i_treatments_sum 10 451.31 469.38 -215.66   431.31
    ##                                               Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS3_isolation_treatments_sum                           
    ## mod_amphibian_SS3_isolation_i_treatments_sum 1.5947      4     0.8097

Effect of agrochemical treatments.

# Fourth Survey (160 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS4 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_amph_SS4_no_effect <- glmmTMB(com_amphibian_ab_SS4 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_amph_SS4_treatments <- glmmTMB(com_amphibian_ab_SS4 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_amph_SS4_isolation <- glmmTMB(com_amphibian_ab_SS4 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_amph_SS4_isolation_treatments <- glmmTMB(com_amphibian_ab_SS4 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_amph_SS4_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS4 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)

model_selection_amphibian_SS4 <- aictab(list(mod_amph_SS4_no_effect,
                                             mod_amph_SS4_treatments,
                                             mod_amph_SS4_isolation,
                                             mod_amph_SS4_isolation_treatments,
                                             mod_amph_SS4_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_amphibian_SS4
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               3 706.31       0.00   0.41 -350.09
    ## Treatments              5 708.61       2.30   0.13 -349.13
    ## Isolation               5 706.70       0.38   0.34 -348.18
    ## Treatments + Isolation  7 708.84       2.52   0.12 -347.09
    ## Treatments * Isolation 11 716.22       9.91   0.00 -346.32

No effects.

LRT for mixed models

``` r
mod_amphibian_SS4_treatments <- glmmTMB(com_amphibian_ab_SS4 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_amphibian_SS4_isolation <- glmmTMB(com_amphibian_ab_SS4 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_amphibian_SS4_isolation_treatments <- glmmTMB(com_amphibian_ab_SS4 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_amphibian_SS4_isolation_i_treatments <- glmmTMB(com_amphibian_ab_SS4 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)

#Isolation
anova(mod_amphibian_SS4_treatments, mod_amphibian_SS4_isolation_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_amphibian_SS4_treatments: com_amphibian_ab_SS4 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS4_isolation_treatments: com_amphibian_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS4_treatments            5 708.27 724.23 -349.13   698.27       
    ## mod_amphibian_SS4_isolation_treatments  7 708.18 730.53 -347.09   694.18 4.0832
    ##                                        Chi Df Pr(>Chisq)
    ## mod_amphibian_SS4_treatments                            
    ## mod_amphibian_SS4_isolation_treatments      2     0.1298

``` r
#Treatments
anova(mod_amphibian_SS4_isolation, mod_amphibian_SS4_isolation_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_amphibian_SS4_isolation: com_amphibian_ab_SS4 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS4_isolation_treatments: com_amphibian_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                        Df    AIC    BIC  logLik deviance  Chisq
    ## mod_amphibian_SS4_isolation             5 706.35 722.32 -348.18   696.35       
    ## mod_amphibian_SS4_isolation_treatments  7 708.18 730.53 -347.09   694.18 2.1702
    ##                                        Chi Df Pr(>Chisq)
    ## mod_amphibian_SS4_isolation                             
    ## mod_amphibian_SS4_isolation_treatments      2     0.3379

``` r
#Interaction
anova(mod_amphibian_SS4_isolation_treatments, mod_amphibian_SS4_isolation_i_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_amphibian_SS4_isolation_treatments: com_amphibian_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_amphibian_SS4_isolation_i_treatments: com_amphibian_ab_SS4 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                          Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS4_isolation_treatments    7 708.18 730.53 -347.09   694.18
    ## mod_amphibian_SS4_isolation_i_treatments 11 714.65 749.77 -346.32   692.65
    ##                                           Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS4_isolation_treatments                           
    ## mod_amphibian_SS4_isolation_i_treatments 1.5349      4     0.8205

No effects.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS4_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_amph_SS4_no_effect_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ 1 , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_amph_SS4_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ treatments , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_amph_SS4_isolation_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ isolation , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_amph_SS4_isolation_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_amph_SS4_isolation_i_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS4_sum)

model_selection_amphibian_SS4_sum <- aictab(list(mod_amph_SS4_no_effect_sum,
                                             mod_amph_SS4_treatments_sum,
                                             mod_amph_SS4_isolation_sum,
                                             mod_amph_SS4_isolation_treatments_sum,
                                             mod_amph_SS4_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_amphibian_SS4_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 305.37       0.00   0.55 -150.54
    ## Treatments              4 308.02       2.65   0.15 -149.51
    ## Isolation               4 306.93       1.56   0.25 -148.97
    ## Treatments + Isolation  6 310.43       5.07   0.04 -148.11
    ## Treatments * Isolation 10 319.09      13.72   0.00 -146.31

No effects.

LRT for mixed models

``` r
mod_amphibian_SS4_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ treatments, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_amphibian_SS4_isolation_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ isolation, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_amphibian_SS4_isolation_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_amphibian_SS4_isolation_i_treatments_sum <- glmmTMB(sum_com_amphibian_ab_SS4 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS4_sum)

#Isolation
anova(mod_amphibian_SS4_treatments_sum, mod_amphibian_SS4_isolation_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_amphibian_SS4_treatments_sum: sum_com_amphibian_ab_SS4 ~ treatments, zi=~0, disp=~1
    ## mod_amphibian_SS4_isolation_treatments_sum: sum_com_amphibian_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS4_treatments_sum            4 307.02 314.25 -149.51   299.02
    ## mod_amphibian_SS4_isolation_treatments_sum  6 308.22 319.06 -148.11   296.22
    ##                                             Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS4_treatments_sum                                   
    ## mod_amphibian_SS4_isolation_treatments_sum 2.7947      2     0.2473

``` r
#Treatments
anova(mod_amphibian_SS4_isolation_sum, mod_amphibian_SS4_isolation_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_amphibian_SS4_isolation_sum: sum_com_amphibian_ab_SS4 ~ isolation, zi=~0, disp=~1
    ## mod_amphibian_SS4_isolation_treatments_sum: sum_com_amphibian_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ##                                            Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS4_isolation_sum             4 305.93 313.16 -148.97   297.93
    ## mod_amphibian_SS4_isolation_treatments_sum  6 308.22 319.06 -148.11   296.22
    ##                                             Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS4_isolation_sum                                    
    ## mod_amphibian_SS4_isolation_treatments_sum 1.7081      2     0.4257

``` r
#Interaction
anova(mod_amphibian_SS4_isolation_treatments_sum, mod_amphibian_SS4_isolation_i_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_amphibian_SS4_isolation_treatments_sum: sum_com_amphibian_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_amphibian_SS4_isolation_i_treatments_sum: sum_com_amphibian_ab_SS4 ~ isolation * treatments, zi=~0, disp=~1
    ##                                              Df    AIC    BIC  logLik deviance
    ## mod_amphibian_SS4_isolation_treatments_sum    6 308.22 319.06 -148.11   296.22
    ## mod_amphibian_SS4_isolation_i_treatments_sum 10 312.62 330.69 -146.31   292.62
    ##                                               Chisq Chi Df Pr(>Chisq)
    ## mod_amphibian_SS4_isolation_treatments_sum                           
    ## mod_amphibian_SS4_isolation_i_treatments_sum 3.6038      4     0.4623

No effects.
