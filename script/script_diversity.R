library(dplyr)
library(readxl)
library(here)
here::i_am("script/script_diversity.R")

data_processed_singular <- read.csv(
  here::here(
    "data_processed",
    "data_processed_singular_area.csv"
  )
)

data_processed_comb <- read.csv(
  here::here(
    "data_processed",
    "data_processed_comb_area.csv"
  )
)

data_processed_arrival <- read.csv(
  here::here(
    "data_processed",
    "data_processed_arrival_area.csv"
  )
)

nonzero_min <- min(data_processed_singular$area[data_processed_singular$area!=0])
# min(data_processed_comb$area[data_processed_comb$area!=0])
# min(data_processed_arrival$area[data_processed_arrival$area!=0])

data_diversity_singular <- data_processed_singular %>%
  group_by(day, condition, treatment, replication) %>%
  filter(sum(area) != 0) %>%
  mutate(
    area=pmax(nonzero_min, area),
    prop=area/sum(area)
  ) %>%
  summarize(
    diversity=-sum(prop[prop!=0]*log(prop[prop!=0]))
  )

data_diversity_comb <- data_processed_comb %>%
  group_by(day, condition, treatment, replication) %>%
  filter(sum(area) != 0) %>%
  mutate(
    area=pmax(nonzero_min, area),
    prop=area/sum(area)
  ) %>%
  summarize(
    diversity=-sum(prop[prop!=0]*log(prop[prop!=0]))
  )

data_diversity_arrival <- data_processed_arrival  %>%
  group_by(day, condition, treatment, replication) %>%
  filter(sum(area) != 0) %>%
  mutate(
    area=pmax(nonzero_min, area),
    prop=area/sum(area)
  ) %>%
  summarize(
    diversity=-sum(prop[prop!=0]*log(prop[prop!=0]))
  )

data_processed_all <- bind_rows(
  data_processed_singular %>% mutate(replication=replication+6),
  data_processed_comb
) %>%
  mutate(
    area=pmax(nonzero_min, area)
  )

data_diversity_all <- bind_rows(
  data_diversity_singular %>% mutate(replication=replication+6),
  data_diversity_comb
)

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
