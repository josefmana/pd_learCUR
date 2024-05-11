
rm( list = ls() )

library("here")
library("tidyverse")
library("dplyr")
library("brms")
library("tidybayes")
library("bayestestR")
library("emmeans")
library("ggpubr")

n.sample = 5000
n.warmup = 5000
n.iter = n.sample + n.warmup
family <- "skew_normal"
link <- "identity"

dat <- read.table( here("_data","Learning-Curve-Data.txt"), header = T) ## load data
dat <- dat[,c(1:which(colnames(dat) == "updrs_III_on"))] 

## check variable types
str(dat)  

## recode
dat$id1 <- as.factor(dat$id1)
dat$pohlavi <- as.factor(dat$pohlavi)
dat$MCI_Control <- recode_factor(dat$MCI_Control, `0` = "Con", `1` = "PD-NC",
                                 `2` = "PD-MCI")
dat$lat <- recode_factor(dat$lat, `0` = "R", `1` = "L", `2` = "O")

## set-up variables
data <- list()
data[["RAVLT"]] <- dat[which(dat$Test == "RAVLT"),]
data[["BVMT"]] <- dat[which(dat$Test == "BVMT"),]
outcome <- c("ravlt", "bvmt")

## outcome distribution
fam <- brmsfamily(family, link = link)

## technical thingies
options(mc.cores = parallel::detectCores())
emm_options(ref_grid = list(level = .89))
options(contrasts = c('contr.bayes', 'contr.poly'))


## models ----
mod <- lapply( setNames( names(data), names(data) ) , function(i)
  brm(Score ~ log(Trial) * MCI_Control + (1|id1),
      iter = n.iter, warmup = n.warmup,
      family = fam,
      data = data[[i]],
      seed = 87542,
      control = list(adapt_delta = .99)))

## save pp-checks
#for (i in 1:length(mod)) {
#  pp_check(mod[[i]], nsamples = 100) +
#    xlab(outcome[i]) +
#    theme(legend.position = 'none')
#  ggsave(filename = paste0(date, "-pp_check-", outcome[i], ".jpg"),
#         path = here("res"),
#         device = "jpeg")
#}

## contrasts ----
slopes <- lapply( setNames( names(mod), names(mod) ), function(i) emtrends(mod[[i]], ~ MCI_Control | Trial, var = "log(Trial)") )
two.way <- lapply( setNames( names(slopes), names(slopes) ), function(i) contrast(slopes[[i]], interaction = "pairwise") )
means <- lapply( setNames( names(mod), names(mod) ), function(i) emmeans(mod[[i]], ~ MCI_Control) )
main <- lapply( setNames( names(means), names(means) ), function(i) pairs(means[[i]]) )
simple <- lapply(setNames( names(mod), names(mod) ), function(i) ref_grid(mod[[i]], at = list(Trial = 1, cov.reduce = MCI_Control ~ Trial)) )
simple <- lapply(setNames( names(simple), names(simple) ), function(i) pairs(simple[[i]]) )

c_all <- list()
for (i in names(mod) ) {
  c_all[[i]] <- list( two.way = two.way[[i]], main = main[[i]], simple = simple[[i]], means = means[[i]], slopes = slopes[[i]] )
}


## results ----

rope <- list( RAVLT = c(-0.223, 0.223), BVMT = c(-0.210, 0.210) )

res <-
  
  lapply(
    
    setNames( names(mod), names(mod) ),
    function(i)
      
      lapply(
        setNames( names(c_all[[i]]), names(c_all[[i]]) ),
        function(j)
          describe_posterior(c_all[[i]][[j]],
                             centrality = c("median", "mean"), dispersion = T,
                             ci = .89, ci_method = "hdi",
                             test = c("p_direction", "rope"),
                             rope_range = rope[[i]],
                             rope_ci = 1
                             ) %>%
          
          mutate( outcome = i, effect = j, .before = 1)
  
    ) %>%
      
      do.call( rbind.data.frame, . )
) %>%
  
  do.call( rbind.data.frame, . )

# save it
write.table( res, here("tabs","Tab-5.csv"), sep = ",", quote = F, row.names = F )


## fig. 2, learning curves ----

outcome.cap <- c("RAVLT", "BVMT-R")

fig.2 <- list()

for (i in 1:length(mod)) {
  mod_grid <- ref_grid(mod[[i]],
                       at = list(Trial = seq(1, max(data[[i]]$Trial), 1),
                                 cov.reduce = MCI_Control ~ Trial))
  mod_emm <- emmeans(mod_grid, ~ Trial * MCI_Control)
  fig.2[[i]] <- as.data.frame(mod_emm) %>% # from here on = plot making
    ggplot(aes(x = Trial, y = emmean, colour = MCI_Control)) +
    geom_point(aes(y = emmean), # points
               position = position_dodge(width = 0.2),
               alpha = 0.4) +
    geom_pointrange(aes(ymax = upper.HPD, ymin = lower.HPD), # whiskers/CIs
                    position = position_dodge(width = 0.2)) +
    geom_line(aes(y = emmean), position = position_dodge(width = 0.2)) +
    scale_color_discrete(name = "Group",
                         labels = as_labeller(c("Con" = "HC",
                                                "PD-NC" = "PD-NC",
                                                "PD-MCI" = "PD-MCI"))) +
    scale_y_continuous(name = paste0(outcome.cap[[i]], " score"),
                       limits = c(0, max(data[[i]]$Score)),
                       breaks = seq(0, max(data[[i]]$Score), 3)) +
    scale_x_continuous(name = "Trial", breaks = seq(1, max(data[[i]]$Trial), 1))
}

# arrange into grid and print
ggarrange(fig.2[[2]], fig.2[[1]], common.legend = T, legend = "right")
ggsave(plot = last_plot(),
       filename = "fig.2.jpg",
       path = here("figs"),
       device = "jpeg")
