library(rstan)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(here)
here::i_am("figure_stanfit/figure_stanfit_default_all_fit.R")

source(here::here("script", "script_diversity.R"))

load(here::here("stanfit", "stanfit_default_all.rda"))

fitdata_group <- fitdata %>%
  group_by(group) %>%
  filter(1:n()==1)

ss <- summary(stanfit_default_all)

plotdata <- data.frame(
  logr=fitdata$r,
  species=fitdata$species,
  est=log(ss$summary[grepl("expected", rownames(ss$summary)),6]),
  lwr=log(ss$summary[grepl("expected", rownames(ss$summary)),4]),
  upr=log(ss$summary[grepl("expected", rownames(ss$summary)),8]),
  treatment=fitdata$treatment,
  condition=fitdata$condition
) %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore", 
                              "ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "ETAxTLExRNE")),
    species=factor(species,
                   levels=1:3,
                   labels=c("italic(L.~minor)", "italic(L.~perpusilla)", "italic(S.~polyrhiza)")),
    condition=factor(condition, 
                     levels=c("hyper", "meso"),
                     labels=c("high-nutrient", "low-nutrient"))
  )

g1 <- ggplot(plotdata) +
  geom_abline(intercept=0, slope=1, lty=2) +
  geom_point(aes(logr, est, col=species, shape=species)) +
  geom_errorbar(aes(ymin=lwr, ymax=upr, x=logr, col=species), alpha=0.5, width=0) +
  scale_x_continuous("Observed growth rates (1/day)") +
  scale_y_continuous("Predicted growth rates (1/day)") +
  scale_color_viridis_d(end=0.8,labels=scales::parse_format()) +
  facet_grid(condition~treatment) +
  scale_shape_discrete(labels=scales::parse_format()) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(
  here::here("figure_stanfit", "figure_stanfit_default_all_fit.pdf"),
  g1,
  width=8,
  height=4
)
ggsave(
  here::here("figure_stanfit", "figure_stanfit_default_all_fit.png"),
  g1,
  width=8,
  height=4
)
