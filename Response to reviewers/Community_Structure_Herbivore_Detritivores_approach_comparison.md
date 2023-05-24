Community Structure - Herbivore and Detritivores
================
Rodolfo Pelinson
27/12/2022

``` r
library(vegan)
library(gllvm)
```

Before anything, here are three functions that will automate most of the
model selection procedures that we will be doing. The script would be
gigantic without them.

``` r
setwd("C:/Users/rodol/OneDrive/Trabalho/Papers/Analysis/IsolationAgrochemicals")
source("Auxiliary scripts/gllvm_AICc_tab.R")
source("Auxiliary scripts/run_multiple_lv.R")
source("Auxiliary scripts/run_multiple_gllvm.R")
source("Auxiliary scripts/get_scaled_lvs.R")
```

# Herbivores and Detritivores

Lets make tables with only herbivore and insectivores. Note that here we
excludes species that occur in less than 4 samples in each survey
because are less informative to community patterns and complicate
parameter estimation.

``` r
com_herb_det_SS1 <- com_SS1[,Trait_SS1$trophic == "consumer"]
com_herb_det_SS2 <- com_SS2[,Trait_SS2$trophic == "consumer"]
com_herb_det_SS3 <- com_SS3[,Trait_SS3$trophic == "consumer"]
com_herb_det_SS4 <- com_SS4[,Trait_SS4$trophic == "consumer"]

sum_com_herb_det_SS2 <- sum_com_SS2[,Trait_SS2_sum$trophic == "consumer"]
sum_com_herb_det_SS3 <- sum_com_SS3[,Trait_SS3_sum$trophic == "consumer"]
sum_com_herb_det_SS4 <- sum_com_SS4[,Trait_SS4_sum$trophic == "consumer"]
```

## First Survey (20 day)

We always used a negative binomial distribution to model species
abundances.

``` r
SS1_predictors <- data.frame(treatments = treatments_SS1,
                             isolation = isolation_SS1)
```

Number of latent variables:

``` r
n_latent_tab_SS1 <- run_multiple_lv(formula = ~ treatments * isolation,
                                    num.lv = c(0,1,2,3),
                                    y = com_herb_det_SS1, X = SS1_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10, starting.val = 0)

n_latent_tab_SS1$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 759.1294    0.00000 40  180
    ## 2     1 772.8654   13.73596 44  180
    ## 3     2 783.7139   24.58444 47  180
    ## 4     3 791.2244   32.09493 49  180

It looks that zero latent variables is the way to go. This will repeat
itself for the next cases.

Model selection of effects:

``` r
model_herb_det_selection_SS1 <- run_multiple_gllvm(formulas = list(~ treatments,
                                                                    ~ isolation,
                                                                    ~ treatments + isolation,
                                                                    ~ treatments * isolation),
                                                    names = c("Treatments",
                                                              "Isolation",
                                                              "Isolation + Treatments",
                                                              "Isolation * Treatments"),
                                                    no_effect = TRUE,
                                                    num.lv = 0,
                                                    X = SS1_predictors,
                                                    y = com_herb_det_SS1,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10, starting.val = 0)

model_herb_det_selection_SS1$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 741.9469   17.07903  8  180
    ## 2             Treatments 753.2960   28.42813 16  180
    ## 3              Isolation 724.8679    0.00000 16  180
    ## 4 Isolation + Treatments 738.1281   13.26021 24  180
    ## 5 Isolation * Treatments 759.1294   34.26158 40  180

It seems that isolation is important. For model selection with mixed
models.

From a frequentist perspective, doing LRT.

``` r
anova(model_herb_det_selection_SS1$models$Treatments, model_herb_det_selection_SS1$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1      164  0.00000       0            
    ## 2      156 35.57243       8 2.10241e-05

``` r
anova(model_herb_det_selection_SS1$models$Isolation, model_herb_det_selection_SS1$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      164 0.000000       0         
    ## 2      156 7.144299       8 0.521147

``` r
anova(model_herb_det_selection_SS1$models$`Isolation + Treatments`, model_herb_det_selection_SS1$models$`Isolation * Treatments`)
```

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff   P.value
    ## 1      156  0.00000       0          
    ## 2      140 26.85382      16 0.0431368

Here the interaction seems to be marginally important. But normally we
should not even consider this possibility as agrochemical treatments
didn’t even started yet.

### I would say that LRT and AICc model selection agree on this one.

## Second Survey (40 day)

``` r
SS2_predictors <- data.frame(ID = ID_SS2_3_4,
                             treatments = treatments_SS2_3_4,
                             isolation = isolation_SS2_3_4)
```

Number of latent variables:

``` r
n_latent_tab_SS2 <- run_multiple_lv(num.lv = c(0,1,2,3),
                                    formula = ~ treatments * isolation,
                                    y = com_herb_det_SS2, X = SS2_predictors,
                                    row.eff = ~ (1|ID),
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


n_latent_tab_SS2$AICc_tab
```

    ##   model     AICc delta_AICc  df nobs
    ## 1     0 3459.311    0.00000  91 1620
    ## 2     1 3479.650   20.33895 100 1620
    ## 3     2 3497.933   38.62252 108 1620
    ## 4     3 3514.091   54.78012 115 1620

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
model_herb_det_selection_SS2 <- run_multiple_gllvm(formulas = list(~ treatments,
                                                                    ~ isolation,
                                                                    ~ treatments + isolation,
                                                                    ~ treatments * isolation),
                                                    names = c("Treatments",
                                                              "Isolation",
                                                              "Isolation + Treatments",
                                                              "Isolation * Treatments"),
                                                    no_effect = TRUE,
                                                    row.eff = ~ (1|ID),
                                                    X_row = data.frame(ID = ID_SS2_3_4),
                                                    num.lv = 0,
                                                    X = SS2_predictors,
                                                    y = com_herb_det_SS2,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 11:20)

model_herb_det_selection_SS2$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 3590.279  150.16775 19 1620
    ## 2             Treatments 3477.169   37.05774 37 1620
    ## 3              Isolation 3542.616  102.50456 37 1620
    ## 4 Isolation + Treatments 3440.112    0.00000 55 1620
    ## 5 Isolation * Treatments 3459.311   19.19902 91 1620

From a frequentist perspective:

``` r
#Isolation
anova(model_herb_det_selection_SS2$models$Treatments, model_herb_det_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1583  0.00000       0            
    ## 2     1565 75.21886      18 5.80101e-09

``` r
#Treatments
anova(model_herb_det_selection_SS2$models$Isolation, model_herb_det_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff P.value
    ## 1     1583   0.0000       0        
    ## 2     1565 140.6657      18       0

``` r
#Interaction
anova(model_herb_det_selection_SS2$models$`Isolation + Treatments`, model_herb_det_selection_SS2$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(model_herb_det_selection_SS2$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff    P.value
    ## 1     1565  0.00000       0           
    ## 2     1529 59.82047      36 0.00758119

Here the interaction also seems to be important. But pay attention to
the warning message.

#### Now, lets repeat this analysis pooling all samples from the same pond together.

``` r
sum_SS2_predictors <- data.frame(ID = ID_SS1,
                             treatments = treatments_SS1,
                             isolation = isolation_SS1)
```

Number of latent variables:

``` r
sum_n_latent_tab_SS2 <- run_multiple_lv(num.lv = c(0,1,2,3),
                                    formula = ~ treatments * isolation,
                                    y = sum_com_herb_det_SS2, X = sum_SS2_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


sum_n_latent_tab_SS2$AICc_tab
```

    ##   model     AICc delta_AICc  df nobs
    ## 1     0 1577.580    0.00000  90  405
    ## 2     1 1608.332   30.75216  99  405
    ## 3     2 1637.232   59.65231 107  405
    ## 4     3 1663.828   86.24792 114  405

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
sum_model_herb_det_selection_SS2 <- run_multiple_gllvm(formulas = list(~ treatments,
                                                                    ~ isolation,
                                                                    ~ treatments + isolation,
                                                                    ~ treatments * isolation),
                                                    names = c("Treatments",
                                                              "Isolation",
                                                              "Isolation + Treatments",
                                                              "Isolation * Treatments"),
                                                    no_effect = TRUE,
                                                    num.lv = 0,
                                                    X = sum_SS2_predictors,
                                                    y = sum_com_herb_det_SS2,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 11:20)

sum_model_herb_det_selection_SS2$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 1592.845   77.78224 18  405
    ## 2             Treatments 1526.222   11.15931 36  405
    ## 3              Isolation 1591.497   76.43472 36  405
    ## 4 Isolation + Treatments 1515.063    0.00000 54  405
    ## 5 Isolation * Treatments 1577.580   62.51734 90  405

Again, additive effects of isolation and treatment.

From a frequentist perspective:

``` r
#Isolation
anova(sum_model_herb_det_selection_SS2$models$Treatments, sum_model_herb_det_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff    P.value
    ## 1      369  0.00000       0           
    ## 2      351 56.89161      18 6.4396e-06

``` r
#Treatments
anova(sum_model_herb_det_selection_SS2$models$Isolation, sum_model_herb_det_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df       D Df.diff P.value
    ## 1      369   0.000       0        
    ## 2      351 122.167      18       0

``` r
#Interaction
anova(sum_model_herb_det_selection_SS2$models$`Isolation + Treatments`, sum_model_herb_det_selection_SS2$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(sum_model_herb_det_selection_SS2$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      351  0.00000       0         
    ## 2      315 44.67684      36 0.152083

Now the interaction seems to not be important according to LRT.

#### Conclusion

Model selection by AICc agreed in only additive effects of isolation and
agrochemical treatments both using mixed models and pooling all samples
from the same ponds together. However, LRT for mixed models pointed out
towards effects of interaction between treatments and isolation. Because
AICc tends to penalize more heavily models with more parameters than LRT
do, I believe we only have enough evidence for additive effects of
treatments in this case.

## Third Survey (80 day)

``` r
SS3_predictors <- data.frame(ID = ID_SS2_3_4,
                             treatments = treatments_SS2_3_4,
                             isolation = isolation_SS2_3_4)
```

Number of latent variables:

``` r
n_latent_tab_SS3 <- run_multiple_lv(num.lv = c(0,1,2,3),
                                    formula = ~ treatments * isolation,
                                    y = com_herb_det_SS3, X = SS3_predictors,
                                    row.eff = ~ (1|ID),
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


n_latent_tab_SS3$AICc_tab
```

    ##   model     AICc delta_AICc  df nobs
    ## 1     0 5990.175    0.00000 101 1800
    ## 2     1 6012.770   22.59567 111 1800
    ## 3     2 6032.330   42.15493 120 1800
    ## 4     3 6050.812   60.63763 128 1800

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
model_herb_det_selection_SS3 <- run_multiple_gllvm(formulas = list(~ treatments,
                                                                    ~ isolation,
                                                                    ~ treatments + isolation,
                                                                    ~ treatments * isolation),
                                                    names = c("Treatments",
                                                              "Isolation",
                                                              "Isolation + Treatments",
                                                              "Isolation * Treatments"),
                                                    no_effect = TRUE,
                                                    row.eff = ~ (1|ID),
                                                    X_row = data.frame(ID = ID_SS2_3_4),
                                                    num.lv = 0,
                                                    X = SS3_predictors,
                                                    y = com_herb_det_SS3,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 11:20)

model_herb_det_selection_SS3$AICc_tab
```

    ##                    model     AICc delta_AICc  df nobs
    ## 1              No Effect 6147.149 156.973830  21 1800
    ## 2             Treatments 6021.149  30.974248  41 1800
    ## 3              Isolation 6112.839 122.664558  41 1800
    ## 4 Isolation + Treatments 5992.850   2.675272  61 1800
    ## 5 Isolation * Treatments 5990.175   0.000000 101 1800

It seems like there is evidence for interactive effects of agrochemicals
and isolation here

From a frequentist perspective:

``` r
#Isolation
anova(model_herb_det_selection_SS3$models$Treatments, model_herb_det_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1759  0.00000       0            
    ## 2     1739 70.69206      20 1.40356e-07

``` r
#Treatments
anova(model_herb_det_selection_SS3$models$Isolation, model_herb_det_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff P.value
    ## 1     1759   0.0000       0        
    ## 2     1739 162.3824      20       0

``` r
#Interaction
anova(model_herb_det_selection_SS3$models$`Isolation + Treatments`, model_herb_det_selection_SS3$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(model_herb_det_selection_SS3$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1739  0.00000       0            
    ## 2     1699 90.45742      40 8.93819e-06

P-values for LRT go to the same direction

#### Now, lets repeat this analysis pooling all samples from the same pond together.

``` r
sum_SS3_predictors <- data.frame(ID = ID_SS1,
                             treatments = treatments_SS1,
                             isolation = isolation_SS1)
```

Number of latent variables:

``` r
sum_n_latent_tab_SS3 <- run_multiple_lv(num.lv = c(0,1,2,3),
                                    formula = ~ treatments * isolation,
                                    y = sum_com_herb_det_SS3, X = sum_SS3_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


sum_n_latent_tab_SS3$AICc_tab
```

    ##   model     AICc delta_AICc  df nobs
    ## 1     0 2502.644    0.00000 100  450
    ## 2     1 2524.019   21.37511 110  450
    ## 3     2 2556.529   53.88516 119  450
    ## 4     3 2586.953   84.30865 127  450

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
sum_model_herb_det_selection_SS3 <- run_multiple_gllvm(formulas = list(~ treatments,
                                                                    ~ isolation,
                                                                    ~ treatments + isolation,
                                                                    ~ treatments * isolation),
                                                    names = c("Treatments",
                                                              "Isolation",
                                                              "Isolation + Treatments",
                                                              "Isolation * Treatments"),
                                                    no_effect = TRUE,
                                                    num.lv = 0,
                                                    X = sum_SS3_predictors,
                                                    y = sum_com_herb_det_SS3,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10)

sum_model_herb_det_selection_SS3$AICc_tab
```

    ##                    model     AICc delta_AICc  df nobs
    ## 1              No Effect 2491.435  55.085925  20  450
    ## 2             Treatments 2439.629   3.280483  40  450
    ## 3              Isolation 2482.105  45.756184  40  450
    ## 4 Isolation + Treatments 2436.349   0.000000  60  450
    ## 5 Isolation * Treatments 2502.644  66.295689 100  450

Now it seems we only have evidence for additive effects.

From a frequentist perspective:

``` r
#Isolation
anova(sum_model_herb_det_selection_SS3$models$Treatments, sum_model_herb_det_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df       D Df.diff     P.value
    ## 1      410  0.0000       0            
    ## 2      390 54.0784      20 5.63194e-05

``` r
#Treatments
anova(sum_model_herb_det_selection_SS3$models$Isolation, sum_model_herb_det_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df       D Df.diff     P.value
    ## 1      410  0.0000       0            
    ## 2      390 96.5541      20 5.18396e-12

``` r
#Interaction
anova(sum_model_herb_det_selection_SS3$models$`Isolation + Treatments`, sum_model_herb_det_selection_SS3$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(sum_model_herb_det_selection_SS3$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff   P.value
    ## 1      390  0.00000       0          
    ## 2      350 52.76649      40 0.0851066

LRT in this case also points towards only additive effects. But the
effect of the interaction could be considered marginally significant.

#### Conclusion

Here mixed models points towards interactive effects of isolation and
treatments. However, pooling samples together only indicate additive
effects. We chose here to with the effects of mixed models as the data
set used in those models may carry more statistical information than
when we pool all samples from the same pond together. Bur we generally
agree that the evidence for interactive effects here is not as strong as
we will se in the next case.

## Fourth Survey (160 day)

``` r
SS4_predictors <- data.frame(ID = ID_SS2_3_4,
                             treatments = treatments_SS2_3_4,
                             isolation = isolation_SS2_3_4)
```

Number of latent variables:

``` r
n_latent_tab_SS4 <- run_multiple_lv(num.lv = c(0,1,2,3),
                                    formula = ~ treatments * isolation,
                                    y = com_herb_det_SS4, X = SS4_predictors,
                                    row.eff = ~ (1|ID),
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


n_latent_tab_SS4$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 4583.193    0.00000 71 1260
    ## 2     1 4599.022   15.82921 78 1260
    ## 3     2 4612.740   29.54720 84 1260
    ## 4     3 4624.279   41.08631 89 1260

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
model_herb_det_selection_SS4 <- run_multiple_gllvm(formulas = list(~ treatments,
                                                                    ~ isolation,
                                                                    ~ treatments + isolation,
                                                                    ~ treatments * isolation),
                                                    names = c("Treatments",
                                                              "Isolation",
                                                              "Isolation + Treatments",
                                                              "Isolation * Treatments"),
                                                    no_effect = TRUE,
                                                    row.eff = ~ (1|ID),
                                                    X_row = data.frame(ID = ID_SS2_3_4),
                                                    num.lv = 0,
                                                    X = SS4_predictors,
                                                    y = com_herb_det_SS4,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10)

model_herb_det_selection_SS4$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 4631.047  47.854593 15 1260
    ## 2             Treatments 4608.898  25.704879 29 1260
    ## 3              Isolation 4599.931  16.738328 29 1260
    ## 4 Isolation + Treatments 4585.070   1.877098 43 1260
    ## 5 Isolation * Treatments 4583.193   0.000000 71 1260

It seems like there is weak evidence for interactive effects of
agrochemicals and isolation here.

From a frequentist perspective:

``` r
#Isolation
anova(model_herb_det_selection_SS4$models$Treatments, model_herb_det_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1231  0.00000       0            
    ## 2     1217 53.52499      14 1.54841e-06

``` r
#Treatments
anova(model_herb_det_selection_SS4$models$Isolation, model_herb_det_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1231  0.00000       0            
    ## 2     1217 44.55844      14 4.80997e-05

``` r
#Interaction
anova(model_herb_det_selection_SS4$models$`Isolation + Treatments`, model_herb_det_selection_SS4$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(model_herb_det_selection_SS4$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1217  0.00000       0            
    ## 2     1189 63.37132      28 0.000148481

P-values for LRT points towards interactive effects, again.

#### Now, lets repeat this analysis pooling all samples from the same pond together.

``` r
sum_SS4_predictors <- data.frame(ID = ID_SS1,
                             treatments = treatments_SS1,
                             isolation = isolation_SS1)
```

Number of latent variables:

``` r
sum_n_latent_tab_SS4 <- run_multiple_lv(num.lv = c(0,1,2,3),
                                    formula = ~ treatments * isolation,
                                    y = sum_com_herb_det_SS4, X = sum_SS4_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


sum_n_latent_tab_SS4$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 1804.377    0.00000 70  315
    ## 2     1 1828.323   23.94587 77  315
    ## 3     2 1850.003   45.62597 83  315
    ## 4     3 1868.949   64.57206 88  315

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
sum_model_herb_det_selection_SS4 <- run_multiple_gllvm(formulas = list(~ treatments,
                                                                    ~ isolation,
                                                                    ~ treatments + isolation,
                                                                    ~ treatments * isolation),
                                                    names = c("Treatments",
                                                              "Isolation",
                                                              "Isolation + Treatments",
                                                              "Isolation * Treatments"),
                                                    no_effect = TRUE,
                                                    num.lv = 0,
                                                    row.eff = FALSE,
                                                    X = sum_SS4_predictors,
                                                    y = sum_com_herb_det_SS4,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10)

sum_model_herb_det_selection_SS4$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 1783.451  6.4863018 14  315
    ## 2             Treatments 1777.148  0.1829907 28  315
    ## 3              Isolation 1778.979  2.0145949 28  315
    ## 4 Isolation + Treatments 1776.965  0.0000000 42  315
    ## 5 Isolation * Treatments 1804.377 27.4121352 70  315

Model selection here points towards only effects of agrochemical
treatments.

From a frequentist perspective:

``` r
#Isolation
anova(sum_model_herb_det_selection_SS4$models$Treatments, sum_model_herb_det_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff    P.value
    ## 1      287  0.00000       0           
    ## 2      273 35.78408      14 0.00112398

``` r
#Treatments
anova(sum_model_herb_det_selection_SS4$models$Isolation, sum_model_herb_det_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff    P.value
    ## 1      287  0.00000       0           
    ## 2      273 37.61568      14 0.00059483

``` r
#Interaction
anova(sum_model_herb_det_selection_SS4$models$`Isolation + Treatments`, sum_model_herb_det_selection_SS4$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(sum_model_herb_det_selection_SS4$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff   P.value
    ## 1      273  0.00000       0          
    ## 2      245 56.04616      28 0.0012698

LRT in this case also points towards interactive effects.

#### Conclusion

Here the only case where interactive effects seems to not be important
happens when we pool all samples from the same pond together and do
model selection through AICc, which penalizes complex models more
strongly.

## Gerneral conclusion

Generally, model selection with samples pooled together seems to be to
conservative, whereas LRTs seems to be the oposit. Both LRT and model
selection through AICc have pros and cons.

- Likelihood Ratio Tests generically don’t penalize for the number of
  parameters in the model and are likely to consider that more
  parameters provide a better fit to the models. Also, they are not
  indicated for gllvms because of the large difference in degrees of
  freedom.

- Model Selection through AICc have a problem of not accurately estimate
  degrees of freedom in mixed models and heavely penalize models with
  too many parameters. However this seems to not be a problem with the
  mixed model approach.

Finally, given these results and the ones for predators, we chose to
interpret the results for model selection the mixed models.
