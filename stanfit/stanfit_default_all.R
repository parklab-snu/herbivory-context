library(rstan)
library(dplyr)
source("../script/script_diversity.R")

fitdata <- lapply(split(data_processed_all, data_processed_all[,c("replication", "treatment", "condition")]), function(x) {
  
  r <- log(tail(arrange(x, species, day)$area, -1)/head(arrange(x, species, day)$area, -1))/diff(arrange(x, species, day)$day)
  
  data.frame(
    r=r[c(1, 2, 4, 5, 7, 8)],
    species=c(1, 1, 2, 2, 3, 3),
    time=c(1, 2, 1, 2, 1, 2),
    n1=rep(filter(x, species=="lm")$area[1:2], 3),
    n2=rep(filter(x, species=="lp")$area[1:2], 3),
    n3=rep(filter(x, species=="sp")$area[1:2], 3),
    replication=x$replication[1],
    treatment=x$treatment[1],
    condition=x$condition[1]
  )  
}) %>%
  bind_rows %>%
  mutate(
    group=as.numeric(factor(interaction(treatment, condition)))
  ) %>%
  filter(is.finite(r)) 

model <- stan_model("../stanmodel/model_default.stan")

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

stanfit_default_all <- sampling(model,
                                 data = standata,
                                 seed=102,
                                 chain=4,
                                 cores=4,
                                 iter=2000)

save("stanfit_default_all", file="stanfit_default_all.rda")

ss <- rstan::summary(stanfit_default_all)

min(ss$summary[,"n_eff"])
max(ss$summary[,"Rhat"])
