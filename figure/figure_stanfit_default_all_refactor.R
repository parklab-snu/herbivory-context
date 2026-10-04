library(rstan)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(here)
here::i_am("figure/figure_stanfit_default_all_refactor.R")

source(here::here("script", "script_diversity.R"))

load(here::here("stanfit", "stanfit_default_all_refactor.rda"))

fitdata_group <- fitdata %>%
  group_by(group) %>%
  filter(1:n()==1)

ss <- summary(stanfit_default_all_refactor)

plotdata <- data.frame(
  species=rep(c("italic(L.~minor)", "italic(L.~perpusilla)", "italic(S.~polyrhiza)"), each=16*3),
  xvar=rep(rep(c("italic(L.~minor)", "italic(L.~perpusilla)", "italic(S.~polyrhiza)"), each=16), 3),
  treatment=rep(fitdata_group$treatment, 9),
  condition=rep(fitdata_group$condition, 9),
  est=ss$summary[grepl("alpha", rownames(ss$summary)),6],
  lwr=ss$summary[grepl("alpha", rownames(ss$summary)),4],
  upr=ss$summary[grepl("alpha", rownames(ss$summary)),8]
) %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore", 
                              "ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "ETAxTLExRNE")),
    text=paste0(xvar,"~-~", species)
  )

plotdata_r <- data.frame(
  species=rep(c("italic(L.~minor)", "italic(L.~perpusilla)", "italic(S.~polyrhiza)"), each=16),
  treatment=rep(fitdata_group$treatment, 3),
  condition=rep(fitdata_group$condition, 3),
  est=ss$summary[grepl("lambda", rownames(ss$summary)),6],
  lwr=ss$summary[grepl("lambda", rownames(ss$summary)),4],
  upr=ss$summary[grepl("lambda", rownames(ss$summary)),8]
)%>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore", 
                              "ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "ETAxTLExRNE"))
  )

g2 <- ggplot(plotdata) +
  geom_point(aes(treatment, est, col=condition), position = position_dodge(width=0.3)) +
  geom_errorbar(aes(treatment, ymin=lwr, ymax=upr, col=condition), width=0.2, position = position_dodge(width=0.3)) +
  geom_text(data=plotdata %>% group_by(text) %>% filter(1:n()==1), aes(0.5, 4.2, label=text), hjust=0, family="Times",
            size=3, parse=TRUE) +
  scale_x_discrete("Treatment") +
  scale_y_continuous("Compeitition coefficient", limits=c(0, NA)) +
  facet_grid(species~xvar) +
  scale_color_manual(values=c("purple", "orange")) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_blank()
  )

g1 <- ggplot(plotdata_r) +
  geom_point(aes(treatment, est, col=condition), position = position_dodge(width=0.3)) +
  geom_errorbar(aes(treatment, ymin=lwr, ymax=upr, col=condition), width=0.2, position = position_dodge(width=0.3)) +
  geom_text(data=plotdata %>% group_by(species) %>% filter(1:n()==1), aes(0.5, 18, label=species), hjust=0, family="Times",
            size=3, parse=TRUE) +
  scale_x_discrete("Treatment") +
  scale_y_continuous("Intrinsic growth rate (1/day)") +
  facet_grid(species~.) +
  scale_color_manual(values=c("purple", "orange")) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

gcomb <- ggarrange(g1, g2, nrow=1, widths=c(1, 3),
                   labels=c("A", "B"))

ggsave(here::here("figure", "figure_stanfit_default_all_refactor.pdf"), gcomb, width=8, height=6)
