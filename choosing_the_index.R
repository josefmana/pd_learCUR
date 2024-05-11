
rm( list = ls() )

library("here")
library("tidyverse")
library("dplyr")
library("brms")
library("bayestestR")
library("cowplot")
library("ggpubr")

if( !dir.exists("figs") ) dir.create("figs")

## reading data ----
dat <- read.table( here("_data","Match_16-01-20.txt"), header = T)

## creating indexes ----
dat$BVMT_L1 <- pmax(dat$BVMT_T2, dat$BVMT_T3) - dat$BVMT_T1
dat$BVMT_L2 <- dat$BVMT_T3/dat$BVMT_T1
dat$BVMT_L3 <- (dat$BVMT_T1 + dat$BVMT_T2 + dat$BVMT_T3) - 3*dat$BVMT_T1
dat$BVMT_L4 <- dat$BVMT_Total * pmax(dat$BVMT_T3-dat$BVMT_T1, dat$BVMT_T2-dat$BVMT_T1)
x <- (1^2+2^2+3^2)-(6^2/3)
dat$BVMT_L5 <- ((dat$BVMT_T1*1 + dat$BVMT_T2*2 + dat$BVMT_T3*3)-((6*dat$BVMT_Total)/3))/x
dat$BVMT_L6 <- dat$BVMT_Total

## recoding
dat$MCI_Control <- recode_factor(dat$MCI_Control, `0` = "Con", `1` = "PD-NC", `2` = "PD-MCI")
dat$pohlavi <- recode(dat$pohlm0z1, "0" = "1", "1" = "0")

## creating sss ----
dum <- seq(0, 36, 1)
ss_1 <- c(2, 2, 4, 5, 6, 8, 9, 9, 10, 12, 13, 15, 17)
ss_2 <- c(2, 2, 2, 2, 4, 5, 6, 7, 8, 9, 10, 12, 14)
ss_3 <- c(2, 2, 2, 2, 2, 4, 5, 6, 6, 7, 9, 10, 13)
ss_tot <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 9,
            9, 9, 10, 10, 11, 12, 12, 13, 14, 15, 17)
ss_dr <- c(2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13)
for (i in 1:nrow(dat)) {
    dat$ss_1[i] <- ss_1[which(dum == dat$BVMT_T1[i])]
    dat$ss_2[i] <- ss_2[which(dum == dat$BVMT_T2[i])]
    dat$ss_3[i] <- ss_3[which(dum == dat$BVMT_T3[i])]
    dat$ss_tot[i] <- ss_tot[which(dum == dat$BVMT_Total[i])]
    dat$ss_dr[i] <- ss_dr[which(dum == dat$BVMT_DR[i])]
}

## calculating z
tab4 <- read.csv( here("tabs","Tab-3.csv") ) %>% column_to_rownames("Index")
dat$pohlavi <- as.numeric(as.character(dat$pohlavi))
dat$BVMT_z_t1 <- (dat$ss_1-(tab4[1,1]+tab4[1,2]*dat$vek_r+tab4[1,3]*dat$vek_r_naDruhou+
                              tab4[1,4]*dat$pohlavi+tab4[1,5]*dat$vzd_r))/tab4[1,7]
dat$BVMT_z_t2 <- (dat$ss_2-(tab4[2,1]+tab4[2,2]*dat$vek_r+tab4[2,3]*dat$vek_r_naDruhou+
                           tab4[2,4]*dat$pohlavi+tab4[2,5]*dat$vzd_r))/tab4[2,7]
dat$BVMT_z_t3 <- (dat$ss_3-(tab4[3,1]+tab4[3,2]*dat$vek_r+tab4[3,3]*dat$vek_r_naDruhou+
                           tab4[3,4]*dat$pohlavi+tab4[3,5]*dat$vzd_r))/tab4[3,7]
dat$BVMT_z_tot <- (dat$ss_tot-(tab4[4,1]+tab4[4,2]*dat$vek_r+
                                 tab4[4,3]*dat$vek_r_naDruhou+tab4[4,4]*dat$pohlavi+
                                 tab4[4,5]*dat$vzd_r))/tab4[4,7]
dat$BVMT_z_dr <- (dat$ss_dr-(tab4[5,1]+tab4[5,2]*dat$vek_r+tab4[5,3]*dat$vek_r_naDruhou+
                             tab4[5,4]*dat$pohlavi+tab4[5,5]*dat$vzd_r))/tab4[5,7]


## multionomial models ----
options(mc.cores = parallel::detectCores()) ## all cores
n.sample = 5000
n.warmup = 5000
n.iter = n.warmup + n.sample
ind <- c("L1", "L3", "L4", "L5", "L6", "z_t1", "z_t2", "z_t3", "z_tot", "z_dr")
mod <- lapply(1:length(ind), function(i)
  brm(as.formula(paste0("MCI_Control ~ BVMT_", ind[i])), 
      family = categorical(refcat = "PD-NC"),
      warmup = n.warmup, iter = n.iter,
      data = dat,
      control = list(adapt_delta = .99),
      seed = 87542))
names(mod) <- ind

## save priors
#for (i in 1:length(mod)) {
#  fname <- paste0(date, "-priors-", ind[i], ".txt")
#  sink(paste(dir, "priors", fname, sep = "/"))
#  print(prior_summary(mod[[i]]))
#  sink()
#}

## add waic and loo information criteria to each model
#for (i in 1:length(ind)) {
#  mod[[i]] <- add_criterion(mod[[i]], criterion = c("loo", "waic"))
#}
## compare models according to loo (same results as waic)
## ordered from the 'best fit' model
## as a heuristic, if elpd_diff > 2*se_diff, the difference between fits is 'significant'
#sink(paste0(dir, "/loo/", date,"-loo.txt"))
#print(loo(mod$L1, mod$L3, mod$L4, mod$L5, mod$L6,
#          mod$z_t1, mod$z_t2, mod$z_t3, mod$z_tot, mod$z_dr))
#sink()
## get looic estimates and se for table
looest <- loo(mod$L1, mod$L3, mod$L4, mod$L5, mod$L6,
              mod$z_t1, mod$z_t2, mod$z_t3, mod$z_tot, mod$z_dr)
tab.loo <- list()
for (i in 1:length(looest[[1]])) {
  tab.loo[[i]] <- looest[[1]][[i]]$estimates[3,]
}
tab.loo <- as.data.frame(do.call(rbind, tab.loo))

## get odds ratios for 'slope' parameters
oddrats <- lapply(1:length(mod), function(i)
  exp(describe_posterior(mod[[i]])[3:4, c("Median", "CI_low", "CI_high")]))
tab.odds <- as.data.frame(do.call(rbind, oddrats))
tab.odds <- cbind(tab.odds[-seq(2, nrow(tab.odds), 2), ],
                  tab.odds[seq(2, nrow(tab.odds), 2), ])

## get full tab5
tab4 <- cbind(tab.odds, tab.loo)
colnames(tab4) <- c("OR.Md.PDNC.vs.Con", "HDI_low", "HDI_high",
                    "OR.Md.PDNC.vs.PDMCI", "HDI_low", "HDI_high",
                    "LOOIC", "LOOCI.SE")
tab4$Index <- ind
tab4 <- tab4[ , c( ncol(tab4), 2:(ncol(tab4)-1) ) ]

write.table( tab4, here("tabs","Tab-4.csv"), sep = ",", quote = F, row.names = F )


## fig. 1, multinomial models ----

xnames <- c(expression("Z"["trial1"]), expression("Z"["trial2"]),
            expression("Z"["trial3"]), expression("Z"["total"]), expression("Z"["DR"]),
            "L1", "L3", "L4", "L5", "L6")

margeffs <- lapply(1:length(mod), function(i)
  plot(marginal_effects(mod[[i]], categorical = T)))

fig.1 <- list()
for (i in 1:length(margeffs)) {
  fig.1[[i]] <- margeffs[[i]][[1]] +
    scale_x_continuous(name = xnames[[i]]) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2)) +
    scale_color_discrete(name = "Group", labels = c("HC", "PD-NC", "PD-MCI")) +
    scale_fill_discrete(name = "Group", labels = c("HC", "PD-NC", "PD-MCI"))
}

# extract the legend
leg <- get_legend(fig.1[[1]])
leg <- as_ggplot(leg)

# arrange (3x3 plot)
ggarrange(fig.1[[1]],
          fig.1[[2]] + theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.ticks.y = element_blank()),
          fig.1[[3]]+theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank()),
          
          fig.1[[4]],
          fig.1[[5]] + theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.ticks.y = element_blank()),
          fig.1[[6]]+theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank()),
          
          fig.1[[7]],
          fig.1[[8]] + theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.ticks.y = element_blank()),
          fig.1[[9]]+theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank()),
          ncol = 3, nrow = 3, common.legend = T, legend = "bottom",
          widths = c(1.2, 1, 1, 1.2, 1, 1, 1.2, 1, 1))

# print (use ca a4 page size)
ggsave(filename = "fig.1.jpg",
       height = 11, width = 8,
       path = here("figs"),
       device = "jpeg")
