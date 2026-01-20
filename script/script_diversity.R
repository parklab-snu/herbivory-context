library(dplyr)
library(readxl)

data_processed_singular <- read.csv("../data_processed/data_processed_singular_area.csv")
data_processed_comb <- read.csv("../data_processed/data_processed_comb_area.csv")
data_processed_arrival <- read.csv("../data_processed/data_processed_arrival_area.csv")

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
