Abundance Analyses - Herbivore and Detritivores - Insects
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

# Herbivores and Detritivores

First, lets make vectors of the abundances (summing up the abundance of
all herb_det taxa in each sample):

``` r
sum_com_herb_det_ab_SS1 <- rowSums(sum_com_orig_SS1[,which(Trait_SS1_sum_orig$trait == "insect_consumer")])
sum_com_herb_det_ab_SS2 <- rowSums(sum_com_orig_SS2[,which(Trait_SS2_sum_orig$trait == "insect_consumer")])
sum_com_herb_det_ab_SS3 <- rowSums(sum_com_orig_SS3[,which(Trait_SS3_sum_orig$trait == "insect_consumer")])
sum_com_herb_det_ab_SS4 <- rowSums(sum_com_orig_SS4[,which(Trait_SS4_sum_orig$trait == "insect_consumer")])

com_herb_det_ab_SS1 <- rowSums(com_SS1_orig[,which(Trait_SS1_orig$trait == "insect_consumer")])
com_herb_det_ab_SS2 <- rowSums(com_SS2_orig[,which(Trait_SS2_orig$trait == "insect_consumer")])
com_herb_det_ab_SS3 <- rowSums(com_SS3_orig[,which(Trait_SS3_orig$trait == "insect_consumer")])
com_herb_det_ab_SS4 <- rowSums(com_SS4_orig[,which(Trait_SS4_orig$trait == "insect_consumer")])
```

# First Survey (20 day)

Model Selection for mixed models

``` r
data_SS1 <- data.frame(isolation = isolation_SS1, treatments = treatments_SS1)

mod_herb_det_SS1_no_effect <- glmmTMB(com_herb_det_ab_SS1 ~ 1, family = nbinom2(link = "log"), data = data_SS1)
mod_herb_det_SS1_treatments <- glmmTMB(com_herb_det_ab_SS1 ~ treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_herb_det_SS1_isolation <- glmmTMB(com_herb_det_ab_SS1 ~ isolation, family = nbinom2(link = "log"), data = data_SS1)
mod_herb_det_SS1_isolation_treatments <- glmmTMB(com_herb_det_ab_SS1 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_herb_det_SS1_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS1 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS1)

model_selection_herb_det_SS1 <- aictab(list(mod_herb_det_SS1_no_effect,
                                             mod_herb_det_SS1_treatments,
                                             mod_herb_det_SS1_isolation,
                                             mod_herb_det_SS1_isolation_treatments,
                                             mod_herb_det_SS1_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_herb_det_SS1
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 301.71       4.43   0.07 -148.71
    ## Treatments              4 304.28       7.00   0.02 -147.64
    ## Isolation               4 297.28       0.00   0.67 -144.14
    ## Treatments + Isolation  6 300.63       3.35   0.13 -143.21
    ## Treatments * Isolation 10 300.85       3.57   0.11 -137.19

Effect of isolation

``` r
mod_herb_det_SS1_treatments <- glmmTMB(com_herb_det_ab_SS1 ~ treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_herb_det_SS1_isolation <- glmmTMB(com_herb_det_ab_SS1 ~ isolation, family = nbinom2(link = "log"), data = data_SS1)
mod_herb_det_SS1_isolation_treatments <- glmmTMB(com_herb_det_ab_SS1 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS1)
mod_herb_det_SS1_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS1 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS1)

#Isolation
anova(mod_herb_det_SS1_treatments, mod_herb_det_SS1_isolation_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_herb_det_SS1_treatments: com_herb_det_ab_SS1 ~ treatments, zi=~0, disp=~1
    ## mod_herb_det_SS1_isolation_treatments: com_herb_det_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance  Chisq
    ## mod_herb_det_SS1_treatments            4 303.28 310.51 -147.64   295.28       
    ## mod_herb_det_SS1_isolation_treatments  6 298.42 309.26 -143.21   286.42 8.8652
    ##                                       Chi Df Pr(>Chisq)  
    ## mod_herb_det_SS1_treatments                              
    ## mod_herb_det_SS1_isolation_treatments      2    0.01188 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_herb_det_SS1_isolation, mod_herb_det_SS1_isolation_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_herb_det_SS1_isolation: com_herb_det_ab_SS1 ~ isolation, zi=~0, disp=~1
    ## mod_herb_det_SS1_isolation_treatments: com_herb_det_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance  Chisq
    ## mod_herb_det_SS1_isolation             4 296.28 303.50 -144.14   288.28       
    ## mod_herb_det_SS1_isolation_treatments  6 298.42 309.26 -143.21   286.42 1.8616
    ##                                       Chi Df Pr(>Chisq)
    ## mod_herb_det_SS1_isolation                             
    ## mod_herb_det_SS1_isolation_treatments      2     0.3942

``` r
#Interaction
anova(mod_herb_det_SS1_isolation_treatments, mod_herb_det_SS1_isolation_i_treatments)
```

    ## Data: data_SS1
    ## Models:
    ## mod_herb_det_SS1_isolation_treatments: com_herb_det_ab_SS1 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_herb_det_SS1_isolation_i_treatments: com_herb_det_ab_SS1 ~ isolation * treatments, zi=~0, disp=~1
    ##                                         Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS1_isolation_treatments    6 298.42 309.26 -143.21   286.42
    ## mod_herb_det_SS1_isolation_i_treatments 10 294.38 312.44 -137.19   274.38
    ##                                          Chisq Chi Df Pr(>Chisq)  
    ## mod_herb_det_SS1_isolation_treatments                             
    ## mod_herb_det_SS1_isolation_i_treatments 12.042      4    0.01704 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Effect of isolation, and a weird effect of interaction.

# Second Survey (40 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS2 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_herb_det_SS2_no_effect <- glmmTMB(com_herb_det_ab_SS2 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_herb_det_SS2_treatments <- glmmTMB(com_herb_det_ab_SS2 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_herb_det_SS2_isolation <- glmmTMB(com_herb_det_ab_SS2 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_herb_det_SS2_isolation_treatments <- glmmTMB(com_herb_det_ab_SS2 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_herb_det_SS2_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS2 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS2)

model_selection_herb_det_SS2 <- aictab(list(mod_herb_det_SS2_no_effect,
                                             mod_herb_det_SS2_treatments,
                                             mod_herb_det_SS2_isolation,
                                             mod_herb_det_SS2_isolation_treatments,
                                             mod_herb_det_SS2_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_herb_det_SS2
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K    AICc Delta_AICc AICcWt      LL
    ## No Effect               3 1072.62      67.53   0.00 -533.24
    ## Treatments              5 1012.24       7.15   0.02 -500.95
    ## Isolation               5 1074.77      69.68   0.00 -532.21
    ## Treatments + Isolation  7 1005.09       0.00   0.80 -495.22
    ## Treatments * Isolation 11 1008.07       2.98   0.18 -492.25

Effect of isolation and agrochemical treatments

LRT for mixed models

``` r
mod_herb_det_SS2_treatments <- glmmTMB(com_herb_det_ab_SS2 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_herb_det_SS2_isolation <- glmmTMB(com_herb_det_ab_SS2 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_herb_det_SS2_isolation_treatments <- glmmTMB(com_herb_det_ab_SS2 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)
mod_herb_det_SS2_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS2 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS2)

#Isolation
anova(mod_herb_det_SS2_treatments, mod_herb_det_SS2_isolation_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_herb_det_SS2_treatments: com_herb_det_ab_SS2 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS2_isolation_treatments: com_herb_det_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance Chisq
    ## mod_herb_det_SS2_treatments            5 1011.9 1027.9 -500.95  1001.90      
    ## mod_herb_det_SS2_isolation_treatments  7 1004.4 1026.8 -495.22   990.44 11.46
    ##                                       Chi Df Pr(>Chisq)   
    ## mod_herb_det_SS2_treatments                               
    ## mod_herb_det_SS2_isolation_treatments      2   0.003247 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_herb_det_SS2_isolation, mod_herb_det_SS2_isolation_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_herb_det_SS2_isolation: com_herb_det_ab_SS2 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS2_isolation_treatments: com_herb_det_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance  Chisq
    ## mod_herb_det_SS2_isolation             5 1074.4 1090.4 -532.21  1064.43       
    ## mod_herb_det_SS2_isolation_treatments  7 1004.4 1026.8 -495.22   990.44 73.986
    ##                                       Chi Df Pr(>Chisq)    
    ## mod_herb_det_SS2_isolation                                 
    ## mod_herb_det_SS2_isolation_treatments      2  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_herb_det_SS2_isolation_treatments, mod_herb_det_SS2_isolation_i_treatments)
```

    ## Data: data_SS2
    ## Models:
    ## mod_herb_det_SS2_isolation_treatments: com_herb_det_ab_SS2 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS2_isolation_i_treatments: com_herb_det_ab_SS2 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                         Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS2_isolation_treatments    7 1004.4 1026.8 -495.22   990.44
    ## mod_herb_det_SS2_isolation_i_treatments 11 1006.5 1041.6 -492.25   984.50
    ##                                          Chisq Chi Df Pr(>Chisq)
    ## mod_herb_det_SS2_isolation_treatments                           
    ## mod_herb_det_SS2_isolation_i_treatments 5.9396      4     0.2037

Effect of isolation and agrochemical treatments

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS2_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_herb_det_SS2_no_effect_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ 1 , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_herb_det_SS2_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ treatments , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_herb_det_SS2_isolation_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ isolation , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_herb_det_SS2_isolation_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS2_sum)
mod_herb_det_SS2_isolation_i_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS2_sum)

model_selection_herb_det_SS2_sum <- aictab(list(mod_herb_det_SS2_no_effect_sum,
                                             mod_herb_det_SS2_treatments_sum,
                                             mod_herb_det_SS2_isolation_sum,
                                             mod_herb_det_SS2_isolation_treatments_sum,
                                             mod_herb_det_SS2_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_herb_det_SS2_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 420.76      50.03   0.00 -208.24
    ## Treatments              4 389.54      18.81   0.00 -190.27
    ## Isolation               4 418.99      48.26   0.00 -204.99
    ## Treatments + Isolation  6 370.73       0.00   0.85 -178.26
    ## Treatments * Isolation 10 374.14       3.40   0.15 -173.83

Effect of isolation and agrochemical treatments

LRT for mixed models

``` r
mod_herb_det_SS2_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ treatments, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_herb_det_SS2_isolation_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ isolation, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_herb_det_SS2_isolation_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS2_sum)
mod_herb_det_SS2_isolation_i_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS2 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS2_sum)

#Isolation
anova(mod_herb_det_SS2_treatments_sum, mod_herb_det_SS2_isolation_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_herb_det_SS2_treatments_sum: sum_com_herb_det_ab_SS2 ~ treatments, zi=~0, disp=~1
    ## mod_herb_det_SS2_isolation_treatments_sum: sum_com_herb_det_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ##                                           Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS2_treatments_sum            4 388.54 395.77 -190.27   380.54
    ## mod_herb_det_SS2_isolation_treatments_sum  6 368.52 379.36 -178.26   356.52
    ##                                            Chisq Chi Df Pr(>Chisq)    
    ## mod_herb_det_SS2_treatments_sum                                       
    ## mod_herb_det_SS2_isolation_treatments_sum 24.018      2   6.09e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_herb_det_SS2_isolation_sum, mod_herb_det_SS2_isolation_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_herb_det_SS2_isolation_sum: sum_com_herb_det_ab_SS2 ~ isolation, zi=~0, disp=~1
    ## mod_herb_det_SS2_isolation_treatments_sum: sum_com_herb_det_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ##                                           Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS2_isolation_sum             4 417.99 425.21 -204.99   409.99
    ## mod_herb_det_SS2_isolation_treatments_sum  6 368.52 379.36 -178.26   356.52
    ##                                            Chisq Chi Df Pr(>Chisq)    
    ## mod_herb_det_SS2_isolation_sum                                        
    ## mod_herb_det_SS2_isolation_treatments_sum 53.466      2  2.454e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_herb_det_SS2_isolation_treatments_sum, mod_herb_det_SS2_isolation_i_treatments_sum)
```

    ## Data: data_SS2_sum
    ## Models:
    ## mod_herb_det_SS2_isolation_treatments_sum: sum_com_herb_det_ab_SS2 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_herb_det_SS2_isolation_i_treatments_sum: sum_com_herb_det_ab_SS2 ~ isolation * treatments, zi=~0, disp=~1
    ##                                             Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS2_isolation_treatments_sum    6 368.52 379.36 -178.26   356.52
    ## mod_herb_det_SS2_isolation_i_treatments_sum 10 367.67 385.73 -173.83   347.67
    ##                                              Chisq Chi Df Pr(>Chisq)  
    ## mod_herb_det_SS2_isolation_treatments_sum                             
    ## mod_herb_det_SS2_isolation_i_treatments_sum 8.8552      4    0.06482 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Effect of isolation and agrochemical treatments.

# Third Survey (80 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS3 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_herb_det_SS3_no_effect <- glmmTMB(com_herb_det_ab_SS3 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_herb_det_SS3_treatments <- glmmTMB(com_herb_det_ab_SS3 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_herb_det_SS3_isolation <- glmmTMB(com_herb_det_ab_SS3 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_herb_det_SS3_isolation_treatments <- glmmTMB(com_herb_det_ab_SS3 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_herb_det_SS3_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS3 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS3)

model_selection_herb_det_SS3 <- aictab(list(mod_herb_det_SS3_no_effect,
                                             mod_herb_det_SS3_treatments,
                                             mod_herb_det_SS3_isolation,
                                             mod_herb_det_SS3_isolation_treatments,
                                             mod_herb_det_SS3_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_herb_det_SS3
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K    AICc Delta_AICc AICcWt      LL
    ## No Effect               3 1746.04      13.02   0.00 -869.95
    ## Treatments              5 1745.35      12.32   0.00 -867.50
    ## Isolation               5 1735.42       2.40   0.21 -862.54
    ## Treatments + Isolation  7 1733.03       0.00   0.70 -859.19
    ## Treatments * Isolation 11 1737.13       4.11   0.09 -856.78

Effect of isolation and agrochemical treatments.

LRT for mixed models

``` r
mod_herb_det_SS3_treatments <- glmmTMB(com_herb_det_ab_SS3 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_herb_det_SS3_isolation <- glmmTMB(com_herb_det_ab_SS3 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_herb_det_SS3_isolation_treatments <- glmmTMB(com_herb_det_ab_SS3 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)
mod_herb_det_SS3_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS3 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS3)

#Isolation
anova(mod_herb_det_SS3_treatments, mod_herb_det_SS3_isolation_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_herb_det_SS3_treatments: com_herb_det_ab_SS3 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS3_isolation_treatments: com_herb_det_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance  Chisq
    ## mod_herb_det_SS3_treatments            5 1745.0 1761.0 -867.50   1735.0       
    ## mod_herb_det_SS3_isolation_treatments  7 1732.4 1754.7 -859.19   1718.4 16.626
    ##                                       Chi Df Pr(>Chisq)    
    ## mod_herb_det_SS3_treatments                                
    ## mod_herb_det_SS3_isolation_treatments      2  0.0002453 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_herb_det_SS3_isolation, mod_herb_det_SS3_isolation_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_herb_det_SS3_isolation: com_herb_det_ab_SS3 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS3_isolation_treatments: com_herb_det_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance  Chisq
    ## mod_herb_det_SS3_isolation             5 1735.1 1751.0 -862.54   1725.1       
    ## mod_herb_det_SS3_isolation_treatments  7 1732.4 1754.7 -859.19   1718.4 6.7015
    ##                                       Chi Df Pr(>Chisq)  
    ## mod_herb_det_SS3_isolation                               
    ## mod_herb_det_SS3_isolation_treatments      2    0.03506 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Interaction
anova(mod_herb_det_SS3_isolation_treatments, mod_herb_det_SS3_isolation_i_treatments)
```

    ## Data: data_SS3
    ## Models:
    ## mod_herb_det_SS3_isolation_treatments: com_herb_det_ab_SS3 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS3_isolation_i_treatments: com_herb_det_ab_SS3 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                         Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS3_isolation_treatments    7 1732.4 1754.7 -859.19   1718.4
    ## mod_herb_det_SS3_isolation_i_treatments 11 1735.6 1770.7 -856.78   1713.6
    ##                                          Chisq Chi Df Pr(>Chisq)
    ## mod_herb_det_SS3_isolation_treatments                           
    ## mod_herb_det_SS3_isolation_i_treatments 4.8134      4      0.307

Effect of isolation and agrochemical treatments.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS3_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_herb_det_SS3_no_effect_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ 1 , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_herb_det_SS3_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ treatments , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_herb_det_SS3_isolation_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ isolation , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_herb_det_SS3_isolation_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS3_sum)
mod_herb_det_SS3_isolation_i_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS3_sum)

model_selection_herb_det_SS3_sum <- aictab(list(mod_herb_det_SS3_no_effect_sum,
                                             mod_herb_det_SS3_treatments_sum,
                                             mod_herb_det_SS3_isolation_sum,
                                             mod_herb_det_SS3_isolation_treatments_sum,
                                             mod_herb_det_SS3_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_herb_det_SS3_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 577.49      12.85   0.00 -286.60
    ## Treatments              4 579.85      15.22   0.00 -285.43
    ## Isolation               4 564.64       0.00   0.80 -277.82
    ## Treatments + Isolation  6 567.64       3.01   0.18 -276.72
    ## Treatments * Isolation 10 572.56       7.92   0.02 -273.04

Effect of isolation only.

LRT for mixed models

``` r
mod_herb_det_SS3_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ treatments, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_herb_det_SS3_isolation_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ isolation, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_herb_det_SS3_isolation_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS3_sum)
mod_herb_det_SS3_isolation_i_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS3 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS3_sum)

#Isolation
anova(mod_herb_det_SS3_treatments_sum, mod_herb_det_SS3_isolation_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_herb_det_SS3_treatments_sum: sum_com_herb_det_ab_SS3 ~ treatments, zi=~0, disp=~1
    ## mod_herb_det_SS3_isolation_treatments_sum: sum_com_herb_det_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ##                                           Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS3_treatments_sum            4 578.85 586.08 -285.43   570.85
    ## mod_herb_det_SS3_isolation_treatments_sum  6 565.43 576.27 -276.72   553.43
    ##                                            Chisq Chi Df Pr(>Chisq)    
    ## mod_herb_det_SS3_treatments_sum                                       
    ## mod_herb_det_SS3_isolation_treatments_sum 17.419      2   0.000165 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Treatments
anova(mod_herb_det_SS3_isolation_sum, mod_herb_det_SS3_isolation_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_herb_det_SS3_isolation_sum: sum_com_herb_det_ab_SS3 ~ isolation, zi=~0, disp=~1
    ## mod_herb_det_SS3_isolation_treatments_sum: sum_com_herb_det_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ##                                           Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS3_isolation_sum             4 563.64 570.86 -277.82   555.64
    ## mod_herb_det_SS3_isolation_treatments_sum  6 565.43 576.27 -276.72   553.43
    ##                                            Chisq Chi Df Pr(>Chisq)
    ## mod_herb_det_SS3_isolation_sum                                    
    ## mod_herb_det_SS3_isolation_treatments_sum 2.2034      2     0.3323

``` r
#Interaction
anova(mod_herb_det_SS3_isolation_treatments_sum, mod_herb_det_SS3_isolation_i_treatments_sum)
```

    ## Data: data_SS3_sum
    ## Models:
    ## mod_herb_det_SS3_isolation_treatments_sum: sum_com_herb_det_ab_SS3 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_herb_det_SS3_isolation_i_treatments_sum: sum_com_herb_det_ab_SS3 ~ isolation * treatments, zi=~0, disp=~1
    ##                                             Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS3_isolation_treatments_sum    6 565.43 576.27 -276.72   553.43
    ## mod_herb_det_SS3_isolation_i_treatments_sum 10 566.09 584.15 -273.04   546.09
    ##                                              Chisq Chi Df Pr(>Chisq)
    ## mod_herb_det_SS3_isolation_treatments_sum                           
    ## mod_herb_det_SS3_isolation_i_treatments_sum 7.3478      4     0.1186

Effect of isolation only.

# Fourth Survey (160 day)

### Mixed models

Model Selection for mixed models

``` r
data_SS4 <- data.frame(ID = ID_SS2_3_4, isolation = isolation_SS2_3_4, treatments = treatments_SS2_3_4)


mod_herb_det_SS4_no_effect <- glmmTMB(com_herb_det_ab_SS4 ~ 1 + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_herb_det_SS4_treatments <- glmmTMB(com_herb_det_ab_SS4 ~ treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_herb_det_SS4_isolation <- glmmTMB(com_herb_det_ab_SS4 ~ isolation + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_herb_det_SS4_isolation_treatments <- glmmTMB(com_herb_det_ab_SS4 ~ isolation + treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_herb_det_SS4_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS4 ~ isolation * treatments + (1|ID), family = nbinom2(link = "log"), data = data_SS4)

model_selection_herb_det_SS4 <- aictab(list(mod_herb_det_SS4_no_effect,
                                             mod_herb_det_SS4_treatments,
                                             mod_herb_det_SS4_isolation,
                                             mod_herb_det_SS4_isolation_treatments,
                                             mod_herb_det_SS4_isolation_i_treatments),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_herb_det_SS4
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K    AICc Delta_AICc AICcWt      LL
    ## No Effect               3 1909.70       4.36   0.09 -951.78
    ## Treatments              5 1909.52       4.18   0.10 -949.59
    ## Isolation               5 1913.13       7.79   0.02 -951.39
    ## Treatments + Isolation  7 1913.02       7.68   0.02 -949.18
    ## Treatments * Isolation 11 1905.34       0.00   0.78 -940.88

Interactive effects.

LRT for mixed models

``` r
mod_herb_det_SS4_treatments <- glmmTMB(com_herb_det_ab_SS4 ~ treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_herb_det_SS4_isolation <- glmmTMB(com_herb_det_ab_SS4 ~ isolation+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_herb_det_SS4_isolation_treatments <- glmmTMB(com_herb_det_ab_SS4 ~ isolation + treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)
mod_herb_det_SS4_isolation_i_treatments <- glmmTMB(com_herb_det_ab_SS4 ~ isolation * treatments+ (1|ID), family = nbinom2(link = "log"), data = data_SS4)

#Isolation
anova(mod_herb_det_SS4_treatments, mod_herb_det_SS4_isolation_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_herb_det_SS4_treatments: com_herb_det_ab_SS4 ~ treatments + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS4_isolation_treatments: com_herb_det_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance  Chisq
    ## mod_herb_det_SS4_treatments            5 1909.2 1925.1 -949.59   1899.2       
    ## mod_herb_det_SS4_isolation_treatments  7 1912.4 1934.7 -949.18   1898.4 0.8066
    ##                                       Chi Df Pr(>Chisq)
    ## mod_herb_det_SS4_treatments                            
    ## mod_herb_det_SS4_isolation_treatments      2     0.6681

``` r
#Treatments
anova(mod_herb_det_SS4_isolation, mod_herb_det_SS4_isolation_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_herb_det_SS4_isolation: com_herb_det_ab_SS4 ~ isolation + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS4_isolation_treatments: com_herb_det_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ##                                       Df    AIC    BIC  logLik deviance  Chisq
    ## mod_herb_det_SS4_isolation             5 1912.8 1928.8 -951.39   1902.8       
    ## mod_herb_det_SS4_isolation_treatments  7 1912.4 1934.7 -949.18   1898.4 4.4132
    ##                                       Chi Df Pr(>Chisq)
    ## mod_herb_det_SS4_isolation                             
    ## mod_herb_det_SS4_isolation_treatments      2     0.1101

``` r
#Interaction
anova(mod_herb_det_SS4_isolation_treatments, mod_herb_det_SS4_isolation_i_treatments)
```

    ## Data: data_SS4
    ## Models:
    ## mod_herb_det_SS4_isolation_treatments: com_herb_det_ab_SS4 ~ isolation + treatments + (1 | ID), zi=~0, disp=~1
    ## mod_herb_det_SS4_isolation_i_treatments: com_herb_det_ab_SS4 ~ isolation * treatments + (1 | ID), zi=~0, disp=~1
    ##                                         Df    AIC    BIC  logLik deviance Chisq
    ## mod_herb_det_SS4_isolation_treatments    7 1912.4 1934.7 -949.18   1898.4      
    ## mod_herb_det_SS4_isolation_i_treatments 11 1903.8 1938.9 -940.88   1881.8  16.6
    ##                                         Chi Df Pr(>Chisq)   
    ## mod_herb_det_SS4_isolation_treatments                       
    ## mod_herb_det_SS4_isolation_i_treatments      4   0.002311 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Interactive effects.

### Models pooling samples together

Model Selection for mixed models

``` r
data_SS4_sum <- data.frame(ID = ID_SS1, isolation = isolation_SS1, treatments = treatments_SS1)

mod_herb_det_SS4_no_effect_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ 1 , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_herb_det_SS4_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ treatments , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_herb_det_SS4_isolation_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ isolation , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_herb_det_SS4_isolation_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ isolation + treatments , family = nbinom2(link = "log"), data = data_SS4_sum)
mod_herb_det_SS4_isolation_i_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ isolation * treatments , family = nbinom2(link = "log"), data = data_SS4_sum)

model_selection_herb_det_SS4_sum <- aictab(list(mod_herb_det_SS4_no_effect_sum,
                                             mod_herb_det_SS4_treatments_sum,
                                             mod_herb_det_SS4_isolation_sum,
                                             mod_herb_det_SS4_isolation_treatments_sum,
                                             mod_herb_det_SS4_isolation_i_treatments_sum),
                                        modnames = c("No Effect",
                                                     "Treatments",
                                                     "Isolation",
                                                     "Treatments + Isolation",
                                                     "Treatments * Isolation"), sort = FALSE)

model_selection_herb_det_SS4_sum
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##                         K   AICc Delta_AICc AICcWt      LL
    ## No Effect               2 604.07       1.40   0.24 -299.89
    ## Treatments              4 604.33       1.65   0.22 -297.66
    ## Isolation               4 608.29       5.62   0.03 -299.65
    ## Treatments + Isolation  6 609.26       6.58   0.02 -297.52
    ## Treatments * Isolation 10 602.68       0.00   0.49 -288.10

Interactive effects is the most plausible model, but the model with no
effects is also plausible.

LRT for mixed models.

``` r
mod_herb_det_SS4_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ treatments, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_herb_det_SS4_isolation_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ isolation, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_herb_det_SS4_isolation_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ isolation + treatments, family = nbinom2(link = "log"), data = data_SS4_sum)
mod_herb_det_SS4_isolation_i_treatments_sum <- glmmTMB(sum_com_herb_det_ab_SS4 ~ isolation * treatments, family = nbinom2(link = "log"), data = data_SS4_sum)

#Isolation
anova(mod_herb_det_SS4_treatments_sum, mod_herb_det_SS4_isolation_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_herb_det_SS4_treatments_sum: sum_com_herb_det_ab_SS4 ~ treatments, zi=~0, disp=~1
    ## mod_herb_det_SS4_isolation_treatments_sum: sum_com_herb_det_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ##                                           Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS4_treatments_sum            4 603.33 610.55 -297.66   595.33
    ## mod_herb_det_SS4_isolation_treatments_sum  6 607.05 617.89 -297.52   595.05
    ##                                            Chisq Chi Df Pr(>Chisq)
    ## mod_herb_det_SS4_treatments_sum                                   
    ## mod_herb_det_SS4_isolation_treatments_sum 0.2783      2     0.8701

``` r
#Treatments
anova(mod_herb_det_SS4_isolation_sum, mod_herb_det_SS4_isolation_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_herb_det_SS4_isolation_sum: sum_com_herb_det_ab_SS4 ~ isolation, zi=~0, disp=~1
    ## mod_herb_det_SS4_isolation_treatments_sum: sum_com_herb_det_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ##                                           Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS4_isolation_sum             4 607.29 614.52 -299.65   599.29
    ## mod_herb_det_SS4_isolation_treatments_sum  6 607.05 617.89 -297.52   595.05
    ##                                            Chisq Chi Df Pr(>Chisq)
    ## mod_herb_det_SS4_isolation_sum                                    
    ## mod_herb_det_SS4_isolation_treatments_sum 4.2459      2     0.1197

``` r
#Interaction
anova(mod_herb_det_SS4_isolation_treatments_sum, mod_herb_det_SS4_isolation_i_treatments_sum)
```

    ## Data: data_SS4_sum
    ## Models:
    ## mod_herb_det_SS4_isolation_treatments_sum: sum_com_herb_det_ab_SS4 ~ isolation + treatments, zi=~0, disp=~1
    ## mod_herb_det_SS4_isolation_i_treatments_sum: sum_com_herb_det_ab_SS4 ~ isolation * treatments, zi=~0, disp=~1
    ##                                             Df    AIC    BIC  logLik deviance
    ## mod_herb_det_SS4_isolation_treatments_sum    6 607.05 617.89 -297.52   595.05
    ## mod_herb_det_SS4_isolation_i_treatments_sum 10 596.20 614.27 -288.10   576.20
    ##                                              Chisq Chi Df Pr(>Chisq)    
    ## mod_herb_det_SS4_isolation_treatments_sum                               
    ## mod_herb_det_SS4_isolation_i_treatments_sum 18.844      4  0.0008433 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Interactive effects.
