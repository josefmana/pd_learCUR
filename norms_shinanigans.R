
rm( list = ls() )

library("here")
library("foreign")
library("tidyverse")
library("openxlsx")
library("bayestestR")
#library("olsrr")
#library("car")

if ( !dir.exists("tabs") ) dir.create("tabs")

norm.dat <- read.spss( here("_data","all_20-04-21.sav"), to.data.frame = T)
norm.dat <- norm.dat[which(norm.dat$without_3SD == "included"), ]

## variables getready
measure <- c("SS_T1_3SD", "SS_T2_3SD", "SS_T3_3SD", "SS_Total_3SD", "SS_DR_3SD")

form <- c(" ~ age",
          " ~ age + age_squared",
          " ~ age + education",
          " ~ age + gender",
          " ~ age + gender + education",
          " ~ age + age_squared + education",
          " ~ age + age_squared + gender",
          " ~ age + age_squared + gender + education")

## fit regressions
mod <- lapply(1:length(measure), function(i)
  lapply(1:length(form), function(j)
    lm(as.formula(paste0(measure[i], form[j])), data = norm.dat)))

names(mod) <- measure

## residual plots
#dia <- c("residuals-v-fitted", "qq", "scale-location", "cook", "residuals-v-leverage")
#for (i in 1:length(mod)) {
#  for (j in 1:length(form)) {
#    for (k in 1:length(dia)) {
#      jpeg(filename = paste0("/res/lm_ss_dia/", dia[k], "/",
#                             measure[i], form[j], ".jpg"),
#           quality = 75)
#      plot(mod[[i]][[j]], which = k)
#      dev.off()
#    }
#  }
#}

## extract bic, aic and r2
comp <- c("aic", "bic", "r2") 

aic <- lapply(1:length(measure), function(i)
  lapply(1:length(form), function(j)
    AIC(mod[[i]][[j]])))

bic <- lapply(1:length(measure), function(i)
  lapply(1:length(form), function(j)
    BIC(mod[[i]][[j]])))

r2 <- lapply(1:length(measure), function(i)
  lapply(1:length(form), function(j)
    summary(mod[[i]][[j]])$r.squared))

comp.res <- list(aic, bic, r2)

for (i in 1:length(comp.res)) {
  for (j in 1:length(comp.res[[1]])) {
    names(comp.res) <- comp
    names(comp.res[[i]]) <- measure
    names(comp.res[[i]][[j]]) <- form
  }
}

res <- lapply( setNames( names(comp.res), names(comp.res) ),
               function(i)
                 as.data.frame(do.call(cbind, comp.res[[i]])) %>%
                 rownames_to_column("model")
               )

## calculate bayes factor approximations
bf <- lapply(1:length(mod), function(i)
  bayesfactor_models(mod[[i]][[2]], mod[[i]][[3]], mod[[i]][[4]],
                     mod[[i]][[5]], mod[[i]][[6]], mod[[i]][[7]],
                     mod[[i]][[8]], denominator = mod[[i]][[1]]))

## write bayes factor approximations to an excel list
bf.tab <- bf[[1]][,2]
for (i in 2:length(bf)) {
  bf.tab <- cbind(bf.tab, bf[[i]][,2])
}

bf.tab <-
  bf.tab %>%
  as.data.frame() %>%
  `rownames<-`( bf[[1]][,1] ) %>%
  `colnames<-`(measure) %>%
  rownames_to_column("model")

res$bf <- bf.tab

## write aic, bic and r2 to excel lists
write.xlsx( res, here("tabs","Tab-comp.xlsx") )

## breusch pagan homoscedasticity tests
#for (i in 1:length(mod)) {
#  sink(paste0("res/breusch-pagan-", measure[i], ".txt"))
#  for (j in 1:length(form)) {
#    print(ols_test_breusch_pagan(mod[[i]][[j]], rhs = T,
#                                 multiple = T, p.adj = "bonferroni"))
#  }
#  sink()
#}

## qqp
#for (i in 1:length(mod)) {
#  for (j in 1:length(form)) {
#    jpeg(filename = paste0(dir, "/res/lm_ss_dia/qqp/", measure[i], form[j], ".jpg"),
#         quality = 200)
#    qqp(mod[[i]][[j]]$residuals)
#    dev.off()
#  }
#}

## get coefficients
coef <- lapply(1:length(mod), function(i)
  mod[[i]][[8]]$coefficients)
names(coef) <- measure

## get r2
rtwo <- lapply(1:length(mod), function(i)
  summary(mod[[i]][[8]])$r.squared)
names(rtwo) <- measure

## get residuals sd
sig <- lapply(1:length(mod), function(i)
  sigma(mod[[i]][[8]]))
names(sig) <- measure

## make tab
tab4 <-
  as.data.frame(do.call(rbind, coef)) %>%
  mutate( R2 = unlist(rtwo) ) %>%
  mutate( sigma = unlist(sig) ) %>%
  rownames_to_column("Index")

## save tab
write.table( tab4, here("tabs","Tab-3.csv"), sep = ",", row.names = F, quote = F )
