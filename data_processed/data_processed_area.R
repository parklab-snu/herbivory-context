library(dplyr)
library(readxl)
library(here)
here::i_am("data_processed/data_processed_area.R")

raw_file <- here::here(
  "data_raw",
  "Choi et al_duckweed and herbivore feeding guild_raw_data.xlsx"
)

data_raw  <- readxl::read_xlsx(raw_file, sheet = 1)
data_raw2 <- readxl::read_xlsx(raw_file, sheet = 2)
data_raw3 <- readxl::read_xlsx(raw_file, sheet = 3)
data_raw4 <- readxl::read_xlsx(raw_file, sheet = 4)
data_raw5 <- readxl::read_xlsx(raw_file, sheet = 5)

data_singular_area <- data.frame(
  day=rep(c(0, 3, 5), each=nrow(data_raw)),
  condition=c(gsub("_.*", "", unname(unlist(data_raw[,1]))), gsub("_.*", "", unname(unlist(data_raw[,7]))), 
              gsub("_.*", "", unname(unlist(data_raw[,13])))),
  treatment=c(gsub(".*_", "", unname(unlist(data_raw[,1]))), gsub(".*_", "", unname(unlist(data_raw[,7]))), 
              gsub(".*_", "", unname(unlist(data_raw[,13])))),
  replication=c(unname(unlist(data_raw[,2])), unname(unlist(data_raw[,8])), unname(unlist(data_raw[,14]))),
  species=c(unname(unlist(data_raw[,3])), unname(unlist(data_raw[,9])), unname(unlist(data_raw[,15]))),
  area=c(unname(unlist(data_raw[,4])), unname(unlist(data_raw[,10])), unname(unlist(data_raw[,16])))
)

data_comb_area_hyper <- data.frame(
  day=rep(c(0, 3, 5), each=nrow(data_raw2)),
  condition="hyper",
  treatment=c(unname(unlist(data_raw2[,1])), unname(unlist(data_raw2[,7])), unname(unlist(data_raw2[,13]))),
  replication=c(unname(unlist(data_raw2[,2])), unname(unlist(data_raw2[,8])), unname(unlist(data_raw2[,14]))),
  species=c(unname(unlist(data_raw2[,3])), unname(unlist(data_raw2[,9])), unname(unlist(data_raw2[,15]))),
  area=c(unname(unlist(data_raw2[,4])), unname(unlist(data_raw2[,10])), unname(unlist(data_raw2[,16])))
)

data_comb_area_meso <- data.frame(
  day=rep(c(0, 3, 5), each=nrow(data_raw3)),
  condition="meso",
  treatment=c(unname(unlist(data_raw3[,1])), unname(unlist(data_raw3[,7])), unname(unlist(data_raw3[,13]))),
  replication=c(unname(unlist(data_raw3[,2])), unname(unlist(data_raw3[,8])), unname(unlist(data_raw3[,14]))),
  species=c(unname(unlist(data_raw3[,3])), unname(unlist(data_raw3[,9])), unname(unlist(data_raw3[,15]))),
  area=c(unname(unlist(data_raw3[,4])), unname(unlist(data_raw3[,10])), unname(unlist(data_raw3[,16])))
)

data_comb_area <- bind_rows(
  data_comb_area_hyper,
  data_comb_area_meso
)

data_arrival_area_hyper <- data.frame(
  day=rep(c(0, 3, 5), each=nrow(data_raw4)),
  condition="hyper",
  treatment=c(unname(unlist(data_raw4[,1])), unname(unlist(data_raw4[,7])), unname(unlist(data_raw4[,13]))),
  replication=c(unname(unlist(data_raw4[,2])), unname(unlist(data_raw4[,8])), unname(unlist(data_raw4[,14]))),
  species=c(unname(unlist(data_raw4[,3])), unname(unlist(data_raw4[,9])), unname(unlist(data_raw4[,15]))),
  area=c(unname(unlist(data_raw4[,4])), unname(unlist(data_raw4[,10])), unname(unlist(data_raw4[,16])))
)

data_arrival_area_meso <- data.frame(
  day=rep(c(0, 3, 5), each=nrow(data_raw5)),
  condition="meso",
  treatment=c(unname(unlist(data_raw5[,1])), unname(unlist(data_raw5[,7])), unname(unlist(data_raw5[,13]))),
  replication=c(unname(unlist(data_raw5[,2])), unname(unlist(data_raw5[,8])), unname(unlist(data_raw5[,14]))),
  species=c(unname(unlist(data_raw5[,3])), unname(unlist(data_raw5[,9])), unname(unlist(data_raw5[,15]))),
  area=c(unname(unlist(data_raw5[,4])), unname(unlist(data_raw5[,10])), unname(unlist(data_raw5[,16])))
)

data_arrival_area <- bind_rows(
  data_arrival_area_hyper,
  data_arrival_area_meso
)

write.csv(
  data_singular_area,
  file = here::here(
    "data_processed",
    "data_processed_singular_area.csv"
  )
)

write.csv(
  data_comb_area,
  file = here::here(
    "data_processed",
    "data_processed_comb_area.csv"
  )
)

write.csv(
  data_arrival_area,
  file = here::here(
    "data_processed",
    "data_processed_arrival_area.csv"
  )
)
