library(tidyr)
library(dplyr)
library(readxl)

data_raw <- read_xlsx("../data_raw/Choi et al_duckweed and herbivore feeding guild_raw_data.xlsx", sheet=9)

data_comb <- bind_rows(
  data.frame(
    treatment="hyper",
    replication=c(data_raw[,2]),
    sp=c(data_raw[,3]),
    lp=c(data_raw[,4]),
    lm=c(data_raw[,5]),
    type="nest"
  ),
  data.frame(
    treatment="hyper",
    replication=c(data_raw[,8]),
    sp=c(data_raw[,9]),
    lp=c(data_raw[,10]),
    lm=c(data_raw[,11]),
    type="feeding"
  ),
  data.frame(
    treatment="meso",
    replication=c(data_raw[,15]),
    sp=c(data_raw[,16]),
    lp=c(data_raw[,17]),
    lm=c(data_raw[,18]),
    type="nest"
  ),
  data.frame(
    treatment="meso",
    replication=c(data_raw[,21]),
    sp=c(data_raw[,22]),
    lp=c(data_raw[,23]),
    lm=c(data_raw[,24]),
    type="feeding"
  )
)

data_processed_eta_pref <- data_comb %>%
  gather(key, value, -treatment, -type,
         -replication) %>%
  mutate(
    value=!is.na(value)
  ) %>%
  group_by(key, treatment, type) %>%
  summarize(
    preference=sum(value)
  )

write.csv(data_processed_eta_pref,
          file="data_processed_eta_pref.csv")
