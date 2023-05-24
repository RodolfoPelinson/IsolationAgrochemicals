Community Structure - Predators
================
Rodolfo Pelinson
27/12/2022

``` r
library(vegan); packageDescription("vegan", fields = "Version")
```

    ## [1] "2.6-2"

``` r
library(gllvm); packageDescription("gllvm", fields = "Version")
```

    ## [1] "1.3.1"

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
com_pred_SS1 <- com_SS1[,Trait_SS1$trophic == "predator"]
com_pred_SS2 <- com_SS2[,Trait_SS2$trophic == "predator"]
com_pred_SS3 <- com_SS3[,Trait_SS3$trophic == "predator"]
com_pred_SS4 <- com_SS4[,Trait_SS4$trophic == "predator"]

sum_com_pred_SS2 <- sum_com_SS2[,Trait_SS2_sum$trophic == "predator"]
sum_com_pred_SS3 <- sum_com_SS3[,Trait_SS3_sum$trophic == "predator"]
sum_com_pred_SS4 <- sum_com_SS4[,Trait_SS4_sum$trophic == "predator"]
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
                                    y = com_pred_SS1, X = SS1_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10, starting.val = 0)

n_latent_tab_SS1$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 498.8127    0.00000 40  180
    ## 2     1 512.5489   13.73620 44  180
    ## 3     2 523.3974   24.58466 47  180
    ## 4     3 530.9079   32.09519 49  180

It looks that zero latent variables is the way to go. This will repeat
itself for the next cases.

Model selection of effects:

``` r
model_pred_selection_SS1 <- run_multiple_gllvm(formulas = list(~ treatments,
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
                                                    y = com_pred_SS1,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10, starting.val = 0)

model_pred_selection_SS1$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 442.4666   0.000000  8  180
    ## 2             Treatments 450.0960   7.629371 16  180
    ## 3              Isolation 447.9068   5.440168 16  180
    ## 4 Isolation + Treatments 456.8692  14.402625 24  180
    ## 5 Isolation * Treatments 498.8127  56.346117 40  180

No effects according to model selection

From a frequentist perspective, doing LRT.

``` r
anova(model_pred_selection_SS1$models$Treatments, model_pred_selection_SS1$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff   P.value
    ## 1      164  0.00000       0          
    ## 2      156 13.63126       8 0.0918974

``` r
anova(model_pred_selection_SS1$models$Isolation, model_pred_selection_SS1$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      164  0.00000       0         
    ## 2      156 11.44206       8 0.177887

``` r
anova(model_pred_selection_SS1$models$`Isolation + Treatments`, model_pred_selection_SS1$models$`Isolation * Treatments`)
```

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      156 0.000000       0         
    ## 2      140 5.911695      16 0.989022

No effects according to LRT as well.

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
                                    y = com_pred_SS2, X = SS2_predictors,
                                    row.eff = ~ (1|ID),
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10, control = list(optim.method = "L-BFGS-B"))


n_latent_tab_SS2$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 1762.985    0.00000 71 1260
    ## 2     1 1778.801   15.81640 78 1260
    ## 3     2 1792.519   29.53458 84 1260
    ## 4     3 1804.058   41.07355 89 1260

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
model_pred_selection_SS2 <- run_multiple_gllvm(formulas = list(~ treatments,
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
                                                    y = com_pred_SS2,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10, control = list(optim.method = "L-BFGS-B"))

model_pred_selection_SS2$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 1761.332  41.483816 15 1260
    ## 2             Treatments 1719.848   0.000000 29 1260
    ## 3              Isolation 1771.902  52.053702 29 1260
    ## 4 Isolation + Treatments 1729.207   9.358933 43 1260
    ## 5 Isolation * Treatments 1762.985  43.136511 71 1260

Only effects of treatments.

From a frequentist perspective:

``` r
#Isolation
anova(model_pred_selection_SS2$models$Treatments, model_pred_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1     1231  0.00000       0         
    ## 2     1217 20.33827      14 0.119832

``` r
#Treatments
anova(model_pred_selection_SS2$models$Isolation, model_pred_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1231  0.00000       0            
    ## 2     1217 72.39198      14 7.09537e-10

``` r
#Interaction
anova(model_pred_selection_SS2$models$`Isolation + Treatments`, model_pred_selection_SS2$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(model_pred_selection_SS2$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1     1217  0.00000       0         
    ## 2     1189 27.71664      28 0.479537

Also only effects of treatments.

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
                                    y = sum_com_pred_SS2, X = sum_SS2_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


sum_n_latent_tab_SS2$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 919.7595    0.00000 70  315
    ## 2     1 943.7053   23.94585 77  315
    ## 3     2 965.3854   45.62594 83  315
    ## 4     3 984.3315   64.57204 88  315

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
sum_model_pred_selection_SS2 <- run_multiple_gllvm(formulas = list(~ treatments,
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
                                                    y = sum_com_pred_SS2,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 11:20)

sum_model_pred_selection_SS2$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 937.0239   88.59773 14  315
    ## 2             Treatments 848.4262    0.00000 28  315
    ## 3              Isolation 955.6005  107.17431 28  315
    ## 4 Isolation + Treatments 866.6334   18.20717 42  315
    ## 5 Isolation * Treatments 919.7595   71.33327 70  315

Only effects of treatments.

From a frequentist perspective:

``` r
#Isolation
anova(sum_model_pred_selection_SS2$models$Treatments, sum_model_pred_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      287  0.00000       0         
    ## 2      273 17.39392      14 0.235793

``` r
#Treatments
anova(sum_model_pred_selection_SS2$models$Isolation, sum_model_pred_selection_SS2$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff P.value
    ## 1      287   0.0000       0        
    ## 2      273 124.5682      14       0

``` r
#Interaction
anova(sum_model_pred_selection_SS2$models$`Isolation + Treatments`, sum_model_pred_selection_SS2$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(sum_model_pred_selection_SS2$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      273  0.00000       0         
    ## 2      245 30.33219      28 0.347517

Also only effects of treatments

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
                                    y = com_pred_SS3, X = SS3_predictors,
                                    row.eff = ~ (1|ID),
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


n_latent_tab_SS3$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 3317.771    0.00000 71 1260
    ## 2     1 3332.137   14.36644 78 1260
    ## 3     2 3345.855   28.08441 84 1260
    ## 4     3 3357.395   39.62353 89 1260

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
model_pred_selection_SS3 <- run_multiple_gllvm(formulas = list(~ treatments,
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
                                                    y = com_pred_SS3,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10)

model_pred_selection_SS3$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 3419.553  130.47324 15 1260
    ## 2             Treatments 3313.546   24.46642 29 1260
    ## 3              Isolation 3397.888  108.80859 29 1260
    ## 4 Isolation + Treatments 3289.079    0.00000 43 1260
    ## 5 Isolation * Treatments 3317.771   28.69171 71 1260

It seems like there is evidence for additive effects of agrochemicals
and isolation here

From a frequentist perspective:

``` r
#Isolation
anova(model_pred_selection_SS3$models$Treatments, model_pred_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1231  0.00000       0            
    ## 2     1217 54.16363      14 1.20445e-06

``` r
#Treatments
anova(model_pred_selection_SS3$models$Isolation, model_pred_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff P.value
    ## 1     1231   0.0000       0        
    ## 2     1217 138.5058      14       0

``` r
#Interaction
anova(model_pred_selection_SS3$models$`Isolation + Treatments`, model_pred_selection_SS3$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(model_pred_selection_SS3$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1     1217  0.00000       0         
    ## 2     1189 32.80251      28 0.243093

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
                                    y = sum_com_pred_SS3, X = sum_SS3_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


sum_n_latent_tab_SS3$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 1440.160    0.00000 70  315
    ## 2     1 1473.934   33.77482 77  315
    ## 3     2 1495.615   55.45494 83  315
    ## 4     3 1514.561   74.40105 88  315

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
sum_model_pred_selection_SS3 <- run_multiple_gllvm(formulas = list(~ treatments,
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
                                                    y = sum_com_pred_SS3,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10)

sum_model_pred_selection_SS3$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 1460.540  63.336955 14  315
    ## 2             Treatments 1397.203   0.000000 28  315
    ## 3              Isolation 1468.474  71.270970 28  315
    ## 4 Isolation + Treatments 1401.327   4.124174 42  315
    ## 5 Isolation * Treatments 1440.160  42.956920 70  315

Now we only have evidence for the effect of agrochemical treatments.

From a frequentist perspective:

``` r
#Isolation
anova(sum_model_pred_selection_SS3$models$Treatments, sum_model_pred_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff   P.value
    ## 1      287  0.00000       0          
    ## 2      273 31.47692      14 0.0047509

``` r
#Treatments
anova(sum_model_pred_selection_SS3$models$Isolation, sum_model_pred_selection_SS3$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1      287   0.0000       0            
    ## 2      273 102.7479      14 1.44329e-15

``` r
#Interaction
anova(sum_model_pred_selection_SS3$models$`Isolation + Treatments`, sum_model_pred_selection_SS3$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(sum_model_pred_selection_SS3$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      273  0.00000       0         
    ## 2      245 44.62555      28 0.024067

LRT, however says we aditive effects AND interactive effects.

#### Conclusion

For mixed models both LRT and model selection agree that we have
additive effects of isolation and agrochemicals. However, when we pool
samples from the same pond together only agrochemicals are important for
model selection and interactive effects are important doing LRT. The
mixed model approach seems to be more reliable.

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
                                    y = com_pred_SS4, X = SS4_predictors,
                                    row.eff = ~ (1|ID),
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


n_latent_tab_SS4$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 3115.011    0.00000 61 1080
    ## 2     1 3128.585   13.57370 67 1080
    ## 3     2 3140.019   25.00867 72 1080
    ## 4     3 3149.250   34.23874 76 1080

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
model_pred_selection_SS4 <- run_multiple_gllvm(formulas = list(~ treatments,
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
                                                    y = com_pred_SS4,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10)

model_pred_selection_SS4$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 3124.320  28.126356 13 1080
    ## 2             Treatments 3119.504  23.310786 25 1080
    ## 3              Isolation 3103.177   6.983416 25 1080
    ## 4 Isolation + Treatments 3096.194   0.000000 37 1080
    ## 5 Isolation * Treatments 3115.011  18.817292 61 1080

Evidence for additive effects of agrochemicals and isolation.

From a frequentist perspective:

``` r
#Isolation
anova(model_pred_selection_SS4$models$Treatments, model_pred_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff     P.value
    ## 1     1055  0.00000       0            
    ## 2     1043 48.77605      12 2.28895e-06

``` r
#Treatments
anova(model_pred_selection_SS4$models$Isolation, model_pred_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff    P.value
    ## 1     1055  0.00000       0           
    ## 2     1043 32.44868      12 0.00117939

``` r
#Interaction
anova(model_pred_selection_SS4$models$`Isolation + Treatments`, model_pred_selection_SS4$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(model_pred_selection_SS4$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff   P.value
    ## 1     1043  0.00000       0          
    ## 2     1019 33.91431      24 0.0862036

Evidence for additive effects as well.

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
                                    y = sum_com_pred_SS4, X = sum_SS4_predictors,
                                    family = "negative.binomial",
                                    method = "VA",
                                    n.init = 10, seed = 1:10)


sum_n_latent_tab_SS4$AICc_tab
```

    ##   model     AICc delta_AICc df nobs
    ## 1     0 1433.050    0.00000 60  270
    ## 2     1 1452.261   19.21087 66  270
    ## 3     2 1470.331   37.28073 71  270
    ## 4     3 1485.457   52.40725 75  270

It looks that zero latent variables is the way to go.

Model selection of effects:

``` r
sum_model_pred_selection_SS4 <- run_multiple_gllvm(formulas = list(~ treatments,
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
                                                    y = sum_com_pred_SS4,
                                                    family = "negative.binomial",
                                                    method = "VA",
                                                    n.init = 10, seed = 1:10)

sum_model_pred_selection_SS4$AICc_tab
```

    ##                    model     AICc delta_AICc df nobs
    ## 1              No Effect 1383.716  0.0000000 12  270
    ## 2             Treatments 1391.799  8.0833725 24  270
    ## 3              Isolation 1384.595  0.8791621 24  270
    ## 4 Isolation + Treatments 1392.869  9.1528734 36  270
    ## 5 Isolation * Treatments 1433.050 49.3339596 60  270

No effects here.

From a frequentist perspective:

``` r
#Isolation
anova(sum_model_pred_selection_SS4$models$Treatments, sum_model_pred_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ treatments 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff    P.value
    ## 1      246  0.00000       0           
    ## 2      234 29.46602      12 0.00335796

``` r
#Treatments
anova(sum_model_pred_selection_SS4$models$Isolation, sum_model_pred_selection_SS4$models$`Isolation + Treatments`)
```

    ## Model  1 :  ~ isolation 
    ## Model  2 :  ~ treatments + isolation

    ##   Resid.Df        D Df.diff   P.value
    ## 1      246  0.00000       0          
    ## 2      234 22.26181      12 0.0346882

``` r
#Interaction
anova(sum_model_pred_selection_SS4$models$`Isolation + Treatments`, sum_model_pred_selection_SS4$models$`Isolation * Treatments`)
```

    ## Warning in anova.gllvm(sum_model_pred_selection_SS4$models$`Isolation + Treatments`, : This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.

    ## Model  1 :  ~ treatments + isolation 
    ## Model  2 :  ~ treatments * isolation

    ##   Resid.Df        D Df.diff  P.value
    ## 1      234  0.00000       0         
    ## 2      210 31.40936      24 0.142327

LRT, however, seems to show evidence for additive effects.

#### Conclusion

Here only model selection for the non-mixed model seems say that there
is no effects of treatments. Again, the model selection with the mixed
model approach seems to be more reliable.

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

Finally, given these results and the ones for herbivores and
detritivores, we chose to interpret the results for model selection the
mixed models.
