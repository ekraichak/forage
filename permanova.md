PERMANOVA on understory plant data
================
Ekaphan Kraichak
2022-11-01

## load the libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6      ✔ purrr   0.3.4 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.0 
    ## ✔ readr   2.1.1      ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-7

## import data

``` r
env <- read.csv("Environmental_matrix.csv")
comm <- read.csv("UnderGrowth_matrix.csv")
```

## clean data

combine env and comm to select only the plot that experience fire to
make sure that they are in the same order. Filter only the plot with
fire

``` r
comm.env <- env %>% 
  select(Plot, Fire, Fire.appear, forest.type) %>%
  mutate(forest.type = recode(forest.type, "MDF,bamboo" = "MDFB", "OpenArea" = "OA")) %>% 
  inner_join(comm, by = c("Plot", "Fire")) %>% 
  filter(Fire.appear == "Fire") %>% 
  filter(!(Plot %in% c(123, 124, 115))) ## these plots are OA and way outliers. Plot 115 is also acting weird

comm.env.E0 <- comm.env[,1:4] ## environmental matrix
comm.env.C0 <- comm.env[,5:484] 
```

after all the filtering, some species disappear from the matrix. We will
select only species that occur in more than 20 plots (remove rare
species)

``` r
main.species <- names(which(colSums(comm.env.C0 > 0) > 20))


comm.env.C1 <- comm.env.C0[,main.species]
## find not empty row
not.empty <- which(rowSums(comm.env.C1) > 0)

comm.env.C <- wisconsin(comm.env.C1[not.empty, main.species])
## community matrix with Wisconsin double standardization
comm.env.E <- comm.env.E0[not.empty, ]
```

## run PERMANOVA using adonis

``` r
perm <- how(nperm = 1999)
setBlocks(perm) <- with(comm.env.E, forest.type)
adonis2(comm.env.C ~ Fire*forest.type, data = comm.env.E, permutations = perm)
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  with(comm.env.E, forest.type) 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = comm.env.C ~ Fire * forest.type, data = comm.env.E, permutations = perm)
    ##                   Df SumOfSqs      R2      F Pr(>F)    
    ## Fire               1    1.320 0.01597 3.3549 0.0005 ***
    ## forest.type        2    6.458 0.07813 8.2094 0.0005 ***
    ## Fire:forest.type   2    0.929 0.01124 1.1810 0.1565    
    ## Residual         188   73.951 0.89466                  
    ## Total            193   82.658 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## visualize the data via MDS

``` r
mod <- metaMDS(comm.env.C, trymax = 100)
```

    ## Run 0 stress 0.286066 
    ## Run 1 stress 0.2909642 
    ## Run 2 stress 0.289493 
    ## Run 3 stress 0.2900714 
    ## Run 4 stress 0.2914391 
    ## Run 5 stress 0.2916312 
    ## Run 6 stress 0.2917299 
    ## Run 7 stress 0.2892094 
    ## Run 8 stress 0.2926075 
    ## Run 9 stress 0.2866541 
    ## Run 10 stress 0.2911262 
    ## Run 11 stress 0.2969354 
    ## Run 12 stress 0.2955816 
    ## Run 13 stress 0.2869332 
    ## Run 14 stress 0.2916523 
    ## Run 15 stress 0.292242 
    ## Run 16 stress 0.2858833 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03719113  max resid 0.1617921 
    ## Run 17 stress 0.2929706 
    ## Run 18 stress 0.2973698 
    ## Run 19 stress 0.2962667 
    ## Run 20 stress 0.2904767 
    ## Run 21 stress 0.2941225 
    ## Run 22 stress 0.2891671 
    ## Run 23 stress 0.2915529 
    ## Run 24 stress 0.2890668 
    ## Run 25 stress 0.2931322 
    ## Run 26 stress 0.2885492 
    ## Run 27 stress 0.2891056 
    ## Run 28 stress 0.2914826 
    ## Run 29 stress 0.2881071 
    ## Run 30 stress 0.2896346 
    ## Run 31 stress 0.290271 
    ## Run 32 stress 0.2866413 
    ## Run 33 stress 0.2918016 
    ## Run 34 stress 0.2937366 
    ## Run 35 stress 0.2921608 
    ## Run 36 stress 0.2901316 
    ## Run 37 stress 0.28974 
    ## Run 38 stress 0.2887213 
    ## Run 39 stress 0.2894043 
    ## Run 40 stress 0.2938599 
    ## Run 41 stress 0.2922708 
    ## Run 42 stress 0.292961 
    ## Run 43 stress 0.293631 
    ## Run 44 stress 0.2867398 
    ## Run 45 stress 0.2890861 
    ## Run 46 stress 0.2887262 
    ## Run 47 stress 0.2936358 
    ## Run 48 stress 0.2892041 
    ## Run 49 stress 0.2911707 
    ## Run 50 stress 0.2911797 
    ## Run 51 stress 0.2869517 
    ## Run 52 stress 0.2884241 
    ## Run 53 stress 0.2929101 
    ## Run 54 stress 0.2880671 
    ## Run 55 stress 0.2870691 
    ## Run 56 stress 0.2911053 
    ## Run 57 stress 0.2879051 
    ## Run 58 stress 0.2882536 
    ## Run 59 stress 0.2945162 
    ## Run 60 stress 0.2889885 
    ## Run 61 stress 0.2851685 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03222848  max resid 0.1614575 
    ## Run 62 stress 0.2884186 
    ## Run 63 stress 0.2917296 
    ## Run 64 stress 0.2881475 
    ## Run 65 stress 0.2881977 
    ## Run 66 stress 0.2929859 
    ## Run 67 stress 0.2887668 
    ## Run 68 stress 0.2910785 
    ## Run 69 stress 0.2904185 
    ## Run 70 stress 0.289787 
    ## Run 71 stress 0.2896176 
    ## Run 72 stress 0.2910058 
    ## Run 73 stress 0.2937205 
    ## Run 74 stress 0.294921 
    ## Run 75 stress 0.2878659 
    ## Run 76 stress 0.2910399 
    ## Run 77 stress 0.2895541 
    ## Run 78 stress 0.2865824 
    ## Run 79 stress 0.2902583 
    ## Run 80 stress 0.2928861 
    ## Run 81 stress 0.2927935 
    ## Run 82 stress 0.2881601 
    ## Run 83 stress 0.2918459 
    ## Run 84 stress 0.2907056 
    ## Run 85 stress 0.2884631 
    ## Run 86 stress 0.2931895 
    ## Run 87 stress 0.2928485 
    ## Run 88 stress 0.2950689 
    ## Run 89 stress 0.2888873 
    ## Run 90 stress 0.295729 
    ## Run 91 stress 0.2888139 
    ## Run 92 stress 0.2889446 
    ## Run 93 stress 0.2946312 
    ## Run 94 stress 0.292748 
    ## Run 95 stress 0.2915824 
    ## Run 96 stress 0.2887879 
    ## Run 97 stress 0.2950357 
    ## Run 98 stress 0.2869257 
    ## Run 99 stress 0.2926437 
    ## Run 100 stress 0.2901003 
    ## *** No convergence -- monoMDS stopping criteria:
    ##      9: no. of iterations >= maxit
    ##     91: stress ratio > sratmax

``` r
mod_scores <- scores(mod)

## without hull
mod_score_env <- cbind(comm.env.E, mod_scores) %>% 
  mutate(forest = paste(forest.type, Fire, sep = " ")) %>% 
  mutate(forest = factor(forest, levels = c("DDF before", "DDF after", "MDF before", "MDF after", "MDFB before", "MDFB after"))) %>% 
  mutate(Fire = factor(Fire, levels = c("before","after")))

## with hull
mod_score_hull <- mod_score_env %>% 
  group_by(Fire, forest.type) %>% 
  slice(chull(NMDS1, NMDS2)) ## create convex hull by groups

p1 <- ggplot() +
  geom_point(aes(x = NMDS1, y = NMDS2, col = forest), data = mod_score_env) +
  geom_polygon(aes(x = NMDS1, y = NMDS2, fill = forest, color = forest), data = mod_score_hull, alpha = 0.3) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired")

p1
```

![](permanova_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## facet by forest types and fire

``` r
p1 +
  facet_grid(~ Fire)
```

![](permanova_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Just the effect of fire

``` r
adonis2(comm.env.C ~ Fire, data = comm.env.E, permutations = perm)
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  with(comm.env.E, forest.type) 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = comm.env.C ~ Fire, data = comm.env.E, permutations = perm)
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## Fire       1    1.320 0.01597 3.1151  5e-04 ***
    ## Residual 192   81.339 0.98403                  
    ## Total    193   82.658 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod_score_hull2 <- mod_score_env %>% 
  group_by(Fire) %>% 
  slice(chull(NMDS1, NMDS2))

ggplot() +
  geom_point(aes(x = NMDS1, y = NMDS2, col = Fire), data = mod_score_env) +
  geom_polygon(aes(x = NMDS1, y = NMDS2, fill = Fire, color = Fire), data = mod_score_hull2, alpha = 0.3)
```

![](permanova_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## PERMANOVA by forest types

``` r
comm.env2 <- cbind(comm.env.E, comm.env.C)

ado_type <- function(TYPE){
  temp.df <- filter(comm.env2, forest.type == TYPE)
  adonis2(temp.df[,5:67] ~ Fire, data = temp.df)
}
```

DDF

``` r
ado_type("DDF")
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = temp.df[, 5:67] ~ Fire, data = temp.df)
    ##          Df SumOfSqs      R2      F Pr(>F)  
    ## Fire      1    0.622 0.01944 1.5861  0.036 *
    ## Residual 80   31.373 0.98056                
    ## Total    81   31.995 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

DDF

``` r
ado_type("MDF")
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = temp.df[, 5:67] ~ Fire, data = temp.df)
    ##          Df SumOfSqs      R2      F Pr(>F)   
    ## Fire      1   0.8316 0.04497 2.0717  0.002 **
    ## Residual 44  17.6616 0.95503                 
    ## Total    45  18.4932 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

DDF

``` r
ado_type("MDFB")
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = temp.df[, 5:67] ~ Fire, data = temp.df)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Fire      1   0.7952 0.03093 2.0425  0.001 ***
    ## Residual 64  24.9164 0.96907                  
    ## Total    65  25.7116 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p1 +
  facet_grid(~ forest.type)
```

![](permanova_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## simper analysis

``` r
ss <- simper(comm.env.C, comm.env.E$Fire, permutations = 99)
summary(ss)
```

    ## 
    ## Contrast: before_after 
    ## 
    ##          average      sd  ratio       ava      avb  cumsum    p   
    ## COMMPA  0.047641 0.07860 0.6061 0.0449602 0.070261 0.05146 0.10 . 
    ## STREJU  0.037959 0.05887 0.6448 0.0435979 0.046953 0.09247 1.00   
    ## CROTHU  0.028539 0.06006 0.4752 0.0157847 0.048675 0.12330 0.03 * 
    ## CYRTAC  0.027901 0.06104 0.4571 0.0550998 0.001505 0.15344 0.01 **
    ## FABACE6 0.024760 0.05644 0.4387 0.0200489 0.033930 0.18019 0.23   
    ## CURCSP1 0.022779 0.04881 0.4667 0.0406495 0.008951 0.20479 0.01 **
    ## ACACCO  0.022560 0.06591 0.3423 0.0227342 0.025163 0.22916 0.69   
    ## OTTONO  0.021659 0.04148 0.5222 0.0291480 0.020719 0.25256 0.69   
    ## PHYLIN  0.021333 0.05081 0.4199 0.0232932 0.022567 0.27561 0.64   
    ## HELIEL  0.021176 0.05392 0.3927 0.0231855 0.022449 0.29848 0.88   
    ## CYPETR  0.020842 0.04117 0.5062 0.0221080 0.024710 0.32100 0.78   
    ## THUNFR  0.020785 0.04712 0.4411 0.0249292 0.020261 0.34345 0.79   
    ## GLOBSP1 0.020214 0.04522 0.4471 0.0348035 0.008766 0.36529 0.02 * 
    ## DIOSSP1 0.018804 0.06232 0.3017 0.0126635 0.026421 0.38560 0.45   
    ## ALBIOD  0.018157 0.05280 0.3438 0.0062094 0.031915 0.40521 0.03 * 
    ## UVARDU  0.017679 0.04115 0.4296 0.0227453 0.015395 0.42431 0.34   
    ## CISSPA  0.017653 0.03733 0.4728 0.0274382 0.011860 0.44338 0.06 . 
    ## COMBPU  0.017135 0.03794 0.4516 0.0161033 0.021327 0.46189 0.67   
    ## CAYRTR  0.017092 0.03793 0.4506 0.0119004 0.025697 0.48035 0.01 **
    ## POLYDE  0.016941 0.05846 0.2898 0.0190434 0.017168 0.49865 0.84   
    ## DALBOL  0.016812 0.03689 0.4557 0.0228681 0.013613 0.51681 0.31   
    ## OCHNIN  0.016241 0.04214 0.3854 0.0214650 0.013188 0.53436 0.30   
    ## OPLICO  0.016070 0.04001 0.4017 0.0281364 0.005552 0.55172 0.01 **
    ## BREYAS  0.015431 0.04521 0.3413 0.0143109 0.018277 0.56839 0.97   
    ## DENDUM  0.014790 0.04631 0.3194 0.0162176 0.014924 0.58436 0.99   
    ## ACRARA  0.014518 0.04475 0.3244 0.0000000 0.029035 0.60005 0.01 **
    ## CROTMA  0.014410 0.03913 0.3682 0.0155852 0.015090 0.61561 0.82   
    ## STERPE  0.014078 0.04910 0.2867 0.0067313 0.022644 0.63082 0.17   
    ## CHROOD  0.014029 0.03466 0.4048 0.0200316 0.010205 0.64598 0.05 * 
    ## MURDME  0.014003 0.04095 0.3419 0.0207835 0.008490 0.66110 0.04 * 
    ## STERGU  0.013065 0.04090 0.3194 0.0087097 0.018786 0.67522 0.42   
    ## MURDLO  0.013027 0.03706 0.3515 0.0089120 0.018654 0.68929 0.16   
    ## CISSRE  0.012984 0.03752 0.3461 0.0118130 0.016357 0.70331 0.95   
    ## KAEMSP1 0.012480 0.02973 0.4197 0.0124685 0.014852 0.71680 0.29   
    ## SHORSI  0.012396 0.03983 0.3112 0.0170469 0.009314 0.73019 0.64   
    ## PTERMA  0.011903 0.03756 0.3169 0.0087019 0.016173 0.74304 0.54   
    ## HELIIS  0.011699 0.03485 0.3357 0.0128309 0.011771 0.75568 0.90   
    ## AMPHMA  0.011550 0.04304 0.2684 0.0151541 0.008547 0.76816 0.01 **
    ## BAUHSP1 0.011490 0.02965 0.3876 0.0121320 0.012919 0.78057 0.84   
    ## UVARCH  0.011277 0.02881 0.3915 0.0122691 0.012072 0.79275 0.53   
    ## HELIAN  0.010843 0.03506 0.3092 0.0190689 0.003391 0.80447 0.02 * 
    ## HEDYOV  0.010385 0.03125 0.3323 0.0156042 0.006317 0.81569 0.26   
    ## PAEDPI  0.010283 0.03532 0.2911 0.0074460 0.014494 0.82679 0.57   
    ## MILLBR  0.010154 0.02684 0.3784 0.0049832 0.016420 0.83776 0.06 . 
    ## TEPHVE  0.010092 0.04680 0.2157 0.0094723 0.011102 0.84866 0.36   
    ## PSEULA  0.009979 0.03138 0.3180 0.0164523 0.004625 0.85944 0.07 . 
    ## TERMNI  0.009882 0.03238 0.3052 0.0063805 0.014226 0.87012 0.34   
    ## RUELPR  0.009387 0.02943 0.3189 0.0074838 0.012095 0.88026 0.42   
    ## RUELRE  0.009286 0.02611 0.3557 0.0139843 0.005441 0.89029 0.04 * 
    ## SCHLOL  0.009163 0.02622 0.3495 0.0079295 0.011662 0.90019 0.53   
    ## SETAPA  0.009002 0.02126 0.4234 0.0076173 0.012147 0.90991 0.07 . 
    ## DIOSDE  0.008725 0.01828 0.4774 0.0087933 0.010538 0.91934 0.87   
    ## XYLIXY  0.008372 0.03738 0.2240 0.0078207 0.009377 0.92838 0.83   
    ## GMELSP1 0.008129 0.03000 0.2710 0.0085788 0.008298 0.93716 0.97   
    ## AMORSP2 0.008049 0.02370 0.3396 0.0080903 0.008852 0.94586 0.86   
    ## CYANCR  0.007792 0.02969 0.2624 0.0000000 0.015585 0.95428 0.01 **
    ## AGERCO  0.007693 0.02458 0.3130 0.0066565 0.009456 0.96259 0.72   
    ## DIOSAL  0.007193 0.04630 0.1553 0.0127806 0.001808 0.97036 0.35   
    ## TRIGRE  0.006657 0.02175 0.3060 0.0063033 0.007592 0.97755 0.50   
    ## CYANBU  0.006448 0.02173 0.2967 0.0078918 0.005600 0.98451 0.43   
    ## CHUKTA  0.006419 0.01699 0.3777 0.0049507 0.008696 0.99145 0.10 . 
    ## ZINGZE  0.005504 0.02357 0.2335 0.0083171 0.003019 0.99739 0.26   
    ## HARRPE  0.002413 0.01764 0.1368 0.0007808 0.004104 1.00000 0.52   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## Permutation: free
    ## Number of permutations: 99

make new data.frame from the results

``` r
ss_df <- data.frame(species = ss$before_after$species,
                    contribution = round(ss$before_after$average,3),
                    abun_before = round(ss$before_after$ava,3), 
                    abun_after = round(ss$before_after$avb,3),
                    p_value = ss$before_after$p)

ss_df %>% 
  filter(p_value <= 0.05)
```

    ##         species contribution abun_before abun_after p_value
    ## AMPHMA   AMPHMA        0.012       0.015      0.009    0.01
    ## CROTHU   CROTHU        0.029       0.016      0.049    0.03
    ## CURCSP1 CURCSP1        0.023       0.041      0.009    0.01
    ## CAYRTR   CAYRTR        0.017       0.012      0.026    0.01
    ## OPLICO   OPLICO        0.016       0.028      0.006    0.01
    ## HELIAN   HELIAN        0.011       0.019      0.003    0.02
    ## GLOBSP1 GLOBSP1        0.020       0.035      0.009    0.02
    ## CHROOD   CHROOD        0.014       0.020      0.010    0.05
    ## ACRARA   ACRARA        0.015       0.000      0.029    0.01
    ## ALBIOD   ALBIOD        0.018       0.006      0.032    0.03
    ## CYANCR   CYANCR        0.008       0.000      0.016    0.01
    ## RUELRE   RUELRE        0.009       0.014      0.005    0.04
    ## MURDME   MURDME        0.014       0.021      0.008    0.04
    ## CYRTAC   CYRTAC        0.028       0.055      0.002    0.01
