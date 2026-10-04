library(rstan)
library(dplyr)
library(here)
here::i_am("stanfit/stanfit_default_all_wide.R")

source(here::here("script", "script_diversity.R"))

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
