library(rstan)
library(dplyr)
library(here)
here::i_am("stanfit/stanfit_default_all_wide.R")

source(here::here("script", "script_diversity.R"))

fitdata <- lapply(split(data_processed_all, data_processed_all[,c("replication", "treatment", "condition")]), function(x) {
  if (nrow(x) > 0) {
    
    n_lm_0 <- filter(x, day==0, species=="lm")$area
    n_lp_0 <- filter(x, day==0, species=="lp")$area
    n_sp_0 <- filter(x, day==0, species=="sp")$area
    
    n_lm_3 <- filter(x, day==3, species=="lm")$area
    n_lp_3 <- filter(x, day==3, species=="lp")$area
    n_sp_3 <- filter(x, day==3, species=="sp")$area
    
    n_lm_5 <- filter(x, day==5, species=="lm")$area
    n_lp_5 <- filter(x, day==5, species=="lp")$area
    n_sp_5 <- filter(x, day==5, species=="sp")$area
    
    r_3 <- log(c(n_lm_3/n_lm_0, n_lp_3/n_lp_0, n_sp_3/n_sp_0))/3
    r_5 <- log(c(n_lm_5/n_lm_3, n_lp_5/n_lp_3, n_sp_5/n_sp_3))/2
    
    data.frame(
      r=c(r_3, r_5),
      species=c(1, 2, 3, 1, 2, 3),
      time=c(1, 1, 1, 2, 2, 2),
      n1=c(n_lm_0, n_lm_0, n_lm_0, n_lm_3, n_lm_3, n_lm_3),
      n2=c(n_lp_0, n_lp_0, n_lp_0, n_lp_3, n_lp_3, n_lp_3),
      n3=c(n_sp_0, n_sp_0, n_sp_0, n_sp_3, n_sp_3, n_sp_3),
      replication=x$replication[1],
      treatment=x$treatment[1],
      condition=x$condition[1]
    )  
  }
}) %>%
  bind_rows %>%
  mutate(
    group=as.numeric(factor(interaction(treatment, condition)))
  ) %>%
  filter(is.finite(r)) 

model <- stan_model(here::here("stanmodel", "model_wide.stan"))

standata <- list(
  N=nrow(fitdata),
  logr=fitdata$r,
  species=fitdata$species,
  n1=fitdata$n1,
  n2=fitdata$n2,
  n3=fitdata$n3,
  Ngroup=length(unique(fitdata$group)),
  group=fitdata$group
)

stanfit_default_all_wide <- sampling(model,
                                 data = standata,
                                 seed=102,
                                 chain=4,
                                 cores=4,
                                 iter=2000)

save(
  "stanfit_default_all_wide",
  file = here::here("stanfit", "stanfit_default_all_wide.rda")
)

ss <- rstan::summary(stanfit_default_all_wide)

min(ss$summary[,"n_eff"]) # 848
max(ss$summary[,"Rhat"]) # 1.0006
