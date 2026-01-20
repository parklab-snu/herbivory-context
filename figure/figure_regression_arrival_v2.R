library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(lme4)
library(ggtext)
source("../script/script_diversity.R")

data_diversity_all_reg <- bind_rows(
  data_diversity_all %>%
    filter(day!=0, !is.na(diversity)) %>%
    mutate(
      ETA=grepl("ETA", treatment),
      TLE=grepl("TLE", treatment),
      RNE=grepl("RNE", treatment),
      arrival=FALSE
    ),
    data_diversity_arrival %>%
      mutate(
        treatment=factor(treatment,
                         levels=c("no herbivore-ETA", 
                                  "RNE-ETA", "TLE-ETA",
                                  "TLExRNE-ETA")),
        replication=replication+12
      ) %>%
      mutate(
        ETA=FALSE,
        TLE=grepl("TLE", treatment),
        RNE=grepl("RNE", treatment),
        arrival=ifelse(day==5, TRUE, FALSE)
      )
  )

lfit_hyp <- lm(diversity~as.factor(day)+(ETA+TLE+RNE)^2+arrival, data=filter(data_diversity_all_reg, condition=="hyper"))
lfit_mes <- lm(diversity~as.factor(day)+(ETA+TLE+RNE)^2+arrival, data=filter(data_diversity_all_reg, condition=="meso"))

data_reg_coef <- bind_rows(
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_hyp))[-(1:3),1],
    lwr=confint(lfit_hyp)[-(1:3),1],
    upr=confint(lfit_hyp)[-(1:3),2],
    condition="eutrophic"
  ),
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_mes))[-(1:3),1],
    lwr=confint(lfit_mes)[-(1:3),1],
    upr=confint(lfit_mes)[-(1:3),2],
    condition="mesotrophic"
  )
) %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "arrival"))
  )

g4 <- ggplot(data_reg_coef) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(aes(est, treatment, col=condition, shape=condition), stat="identity", size=2,
           position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin=lwr, xmax=upr, y=treatment, group=condition, col=condition),
                 height=0.2,
                 position=position_dodge(width=0.5)) +
  scale_x_continuous("Effect on Shannon diversity") +
  scale_y_discrete("Treatment",limits=rev) +
  scale_fill_manual(values=c("purple", "orange")) +
  scale_color_manual(values=c("purple", "orange")) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.border = element_rect(linewidth=1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown()
  )

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
  filter(is.finite(r)) %>%
  mutate(
    ETA=grepl("ETA", treatment),
    TLE=grepl("TLE", treatment),
    RNE=grepl("RNE", treatment),
    arrival=FALSE
  )

fitdata_arrival <- lapply(split(data_processed_arrival, data_processed_arrival[,c("replication", "treatment", "condition")]), function(x) {
  
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
    condition=x$condition[1],
    arrival=c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
  )  
}) %>%
  bind_rows %>%
  mutate(
    group=as.numeric(factor(interaction(treatment, condition)))
  ) %>%
  filter(is.finite(r)) %>%
  mutate(
    ETA=FALSE,
    TLE=grepl("TLE", treatment),
    RNE=grepl("RNE", treatment)
  )

fitdata_comb <- bind_rows(
  fitdata,
  fitdata_arrival %>% mutate(replication=replication+12)
)

lfit_hyp_r_lm <- lm(r~as.factor(time)+(ETA+TLE+RNE)^2+arrival, data=filter(fitdata_comb, condition=="hyper", species==1))
lfit_hyp_r_lp <- lm(r~as.factor(time)+(ETA+TLE+RNE)^2+arrival, data=filter(fitdata_comb, condition=="hyper", species==2))
lfit_hyp_r_sp <- lm(r~as.factor(time)+(ETA+TLE+RNE)^2+arrival, data=filter(fitdata_comb, condition=="hyper", species==3))

lfit_mes_r_lm <- lm(r~as.factor(time)+(ETA+TLE+RNE)^2+arrival, data=filter(fitdata_comb, condition=="meso", species==1))
lfit_mes_r_lp <- lm(r~as.factor(time)+(ETA+TLE+RNE)^2+arrival, data=filter(fitdata_comb, condition=="meso", species==2))
lfit_mes_r_sp <- lm(r~as.factor(time)+(ETA+TLE+RNE)^2+arrival, data=filter(fitdata_comb, condition=="meso", species==3))

data_reg_coef_r <- bind_rows(
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_hyp_r_lm))[-(1:2),1],
    lwr=confint(lfit_hyp_r_lm)[-(1:2),1],
    upr=confint(lfit_hyp_r_lm)[-(1:2),2],
    species="lm",
    condition="eutrophic"
  ),
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_hyp_r_lp))[-(1:2),1],
    lwr=confint(lfit_hyp_r_lp)[-(1:2),1],
    upr=confint(lfit_hyp_r_lp)[-(1:2),2],
    species="lp",
    condition="eutrophic"
  ),
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_hyp_r_sp))[-(1:2),1],
    lwr=confint(lfit_hyp_r_sp)[-(1:2),1],
    upr=confint(lfit_hyp_r_sp)[-(1:2),2],
    species="sp",
    condition="eutrophic"
  ),
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_mes_r_lm))[-(1:2),1],
    lwr=confint(lfit_mes_r_lm)[-(1:2),1],
    upr=confint(lfit_mes_r_lm)[-(1:2),2],
    species="lm",
    condition="mesotrophic"
  ),
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_mes_r_lp))[-(1:2),1],
    lwr=confint(lfit_mes_r_lp)[-(1:2),1],
    upr=confint(lfit_mes_r_lp)[-(1:2),2],
    species="lp",
    condition="mesotrophic"
  ),
  data.frame(
    treatment=c("ETA", "TLE", "RNE", "arrival", "ETAxTLE", "ETAxRNE", "TLExRNE"),
    est=coef(summary(lfit_mes_r_sp))[-(1:2),1],
    lwr=confint(lfit_mes_r_sp)[-(1:2),1],
    upr=confint(lfit_mes_r_sp)[-(1:2),2],
    species="sp",
    condition="mesotrophic"
  )
) %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "arrival"))
  )

g1 <- ggplot(data_reg_coef_r %>% filter(species=="lm")) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(aes(est, treatment, col=condition, shape=condition), stat="identity", size=2,
           position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin=lwr, xmax=upr, y=treatment, group=condition, col=condition),
                 height=0.2,
                 position=position_dodge(width=0.5)) +
  scale_x_continuous("Growth rate effect (1/day)", limits=c(-1.3, 1.3)) +
  scale_y_discrete("Treatment",limits=rev) +
  scale_fill_manual(values=c("purple", "orange")) +
  scale_color_manual(values=c("purple", "orange")) +
  ggtitle(expression(italic(L.~minor))) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.border = element_rect(linewidth=1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown()
  )

g2 <- ggplot(data_reg_coef_r %>% filter(species=="lp")) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(aes(est, treatment, col=condition, shape=condition), stat="identity", size=2,
             position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin=lwr, xmax=upr, y=treatment, group=condition, col=condition),
                 height=0.2,
                 position=position_dodge(width=0.5)) +
  scale_x_continuous("Growth rate effect (1/day)", limits=c(-1.3, 1.3)) +
  scale_y_discrete("Treatment",limits=rev) +
  scale_fill_manual(values=c("purple", "orange")) +
  scale_color_manual(values=c("purple", "orange")) +
  ggtitle(expression(italic(L.~perpusilla))) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.border = element_rect(linewidth=1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown()
  )

g3 <- ggplot(data_reg_coef_r %>% filter(species=="sp")) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(aes(est, treatment, col=condition, shape=condition), stat="identity", size=2,
             position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin=lwr, xmax=upr, y=treatment, group=condition, col=condition),
                 height=0.2,
                 position=position_dodge(width=0.5)) +
  scale_x_continuous("Growth rate effect (1/day)", limits=c(-1.3, 1.3)) +
  scale_y_discrete("Treatment",limits=rev) +
  scale_fill_manual(values=c("purple", "orange")) +
  scale_color_manual(values=c("purple", "orange")) +
  ggtitle(expression(italic(S.~polyrhiza))) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.border = element_rect(linewidth=1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown()
  )

gcomb <- ggarrange(g1, g2, g3, g4, nrow=1, labels=c("A", "B", "C", "D"))

ggsave("figure_regression_arrival_v2.pdf", gcomb, width=12, height=4)
ggsave("figure_regression_arrival_v2.png", gcomb, width=12, height=4)

data_reg_coef_r_eta <- data_reg_coef_r %>%
  filter(treatment=="ETA")

data_processed_eta_pref <- read.csv("../data_processed/data_processed_eta_pref.csv") %>%
  rename(
    species=key,
    condition=treatment
  ) %>%
  mutate(
    condition=factor(condition, levels=c("meso", "hyper"),
                     labels=c("mesotrophic", "eutrophic"))
  )

data_reg_coef_r_eta_pref <- merge(data_reg_coef_r_eta,
      data_processed_eta_pref) %>%
  mutate(
    species=factor(species,
                   levels=c("lm", "lp", "sp"),
                   labels=c("L. minor", "L. perpusilla",
                            "S. polyrhiza"))
  )

g5 <- ggplot(data_reg_coef_r_eta_pref %>% filter(type=="nest", condition=="eutrophic")) +
  geom_bar(aes(species, preference), stat="identity", color="black", fill=NA) +
  geom_point(aes(species, -est*30), col="red") +
  geom_errorbar(aes(species, ymin=-lwr*30, ymax=-upr*30), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Nesting preference", expand=c(0, 0), limits=c(0, 34),
                     sec.axis = sec_axis(transform=~-./30,
                                         "Mean effect on growth rate")) +
  ggtitle("Nesting preference, high-nutrient") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    axis.line.y.right = element_line(color="red"),
    axis.ticks.y.right = element_line(color="red"),
    axis.text.y.right = element_text(color="red"),
    axis.title.y.right = element_text(color="red"),
    axis.text.x = element_text(face = "italic")
  )

g6 <- ggplot(data_reg_coef_r_eta_pref %>% filter(type=="nest", condition=="mesotrophic")) +
  geom_bar(aes(species, preference), stat="identity", color="black", fill=NA) +
  geom_point(aes(species, -est*10), col="red") +
  geom_errorbar(aes(species, ymin=-lwr*10, ymax=-upr*10), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Nesting preference", expand=c(0, 0), limits=c(0, 34),
                     sec.axis = sec_axis(transform=~-./10,
                                         "Mean effect on growth rate")) +
  ggtitle("Nesting preference, low-nutrient") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    axis.line.y.right = element_line(color="red"),
    axis.ticks.y.right = element_line(color="red"),
    axis.text.y.right = element_text(color="red"),
    axis.title.y.right = element_text(color="red"),
    axis.text.x = element_text(face = "italic")
  )

g7 <- ggplot(data_reg_coef_r_eta_pref %>% filter(type=="feeding", condition=="eutrophic")) +
  geom_bar(aes(species, preference), stat="identity", color="black", fill=NA) +
  geom_point(aes(species, -est*30), col="red") +
  geom_errorbar(aes(species, ymin=-lwr*30, ymax=-upr*30), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Feeding preference", expand=c(0, 0), limits=c(0, 34),
                     sec.axis = sec_axis(transform=~-./30,
                                         "Mean effect on growth rate")) +
  ggtitle("Feeding preference, high-nutrient") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    axis.line.y.right = element_line(color="red"),
    axis.ticks.y.right = element_line(color="red"),
    axis.text.y.right = element_text(color="red"),
    axis.title.y.right = element_text(color="red"),
    axis.text.x = element_text(face = "italic")
  )

g8 <- ggplot(data_reg_coef_r_eta_pref %>% filter(type=="feeding", condition=="mesotrophic")) +
  geom_bar(aes(species, preference), stat="identity", color="black", fill=NA) +
  geom_point(aes(species, -est*10), col="red") +
  geom_errorbar(aes(species, ymin=-lwr*10, ymax=-upr*10), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Feeding preference", expand=c(0, 0), limits=c(0, 34),
                     sec.axis = sec_axis(transform=~-./10,
                                         "Mean effect on growth rate")) +
  ggtitle("Feeding preference, low-nutrient") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    axis.line.y.right = element_line(color="red"),
    axis.ticks.y.right = element_line(color="red"),
    axis.text.y.right = element_text(color="red"),
    axis.title.y.right = element_text(color="red"),
    axis.text.x = element_text(face = "italic")
  )

gpref <- ggarrange(g5, g6, g7, g8, nrow=2,
          labels=c("A", "B", "C", "D"))

ggsave("figure_regression_arrival_pref.pdf", gpref, width=6*1.1, height=4*1.1)
ggsave("figure_regression_arrival_pref.png", gpref, width=6*1.1, height=4*1.1)
