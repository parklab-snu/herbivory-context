library(rstan)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(ggdensity)
source("../script/script_diversity.R")

load("../stanfit/stanfit_default_all.rda")

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

fitdata_group <- fitdata %>%
  group_by(group) %>%
  filter(1:n()==1)

ss <- summary(stanfit_default_all)
ee <- rstan::extract(stanfit_default_all)

dim(ee$ndiff)

reslist_i <- vector('list', 4000)
for (i in 1:4000) {
  
  reslist_j <- vector('list', 3)
  for (j in 1:3) {
    
    reslist_k <- vector('list', 16)
    for (k in 1:16) {
      print(paste0(i, ", ", j, ", ", k))
      rho <- ee$rho[i,j,k]
      fdiff <- ee$fdiff[i,j,k]
      
      fdiff <- pmax(fdiff, 1/fdiff)
      
      if (rho < fdiff & fdiff < 1/rho) {
        outcome <- "coexistence"
      } else if (rho > fdiff & fdiff > 1/rho) {
        outcome <- "priority effect"
      } else {
        outcome <- "competitive exclusion"
      }
      
      reslist_k[[k]] <- data.frame(
        sample=i,
        pair=j,
        group=k,
        outcome=outcome,
        rho=rho,
        fdiff=fdiff
      )
      
    }
    
    reslist_j[[j]] <- reslist_k %>%
      bind_rows
  }
  reslist_i[[i]] <- reslist_j %>%
    bind_rows
}

reslist_all <- reslist_i %>%
  bind_rows

summlist_all <- reslist_all %>%
  group_by(pair, group, outcome) %>%
  summarize(
    tot=n()
  ) %>%
  merge(select(fitdata_group, treatment, condition, group)) %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore", 
                              "ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "ETAxTLExRNE")),
    pair=factor(pair,
                labels=c("italic(L.~minor)~vs~italic(L.~perpusilla)", 
                         "italic(L.~minor)~vs~italic(S.~polyrhiza)", 
                         "italic(L.~perpusilla)~vs~italic(S.~polyrhiza)")),
    outcome=factor(outcome,
                   levels=rev(c("coexistence", "priority effect",
                            "competitive exclusion"))),
    condition=factor(condition,
                     levels=c("hyper", "meso"),
                     labels=c("eutrophic", "mesotrophic"))
  )

g1 <- ggplot(summlist_all) +
  geom_bar(aes(treatment, tot/4000, fill=outcome), stat="identity") +
  scale_y_continuous("Probability of competition outcome", expand=c(0, 0)) +
  scale_fill_viridis_d(end=0.8) +
  facet_grid(condition~pair, labeller = label_parsed) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave("figure_stanfit_default_niche.pdf", g1, width=6, height=6)
ggsave("figure_stanfit_default_niche.png", g1, width=6, height=6)

reslist_all2 <- reslist_all %>%
  merge(select(fitdata_group, treatment, condition, group)) %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore", 
                              "ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "ETAxTLExRNE")),
    pair=factor(pair,
                labels=c("italic(L.~minor)~vs~italic(L.~perpusilla)", 
                         "italic(L.~minor)~vs~italic(S.~polyrhiza)", 
                         "italic(L.~perpusilla)~vs~italic(S.~polyrhiza)")),
    outcome=factor(outcome,
                   levels=rev(c("coexistence", "priority effect",
                                "competitive exclusion"))),
    condition=factor(condition,
                     levels=c("hyper", "meso"),
                     labels=c("high-nutrient", "low-nutrient"))
  )

plotdata1 <- data.frame(
  x=seq(0, 9, length.out=21),
  ymax=seq(0, 9, length.out=21),
  ymin=0
)

plotdata2 <- data.frame(
  x=-seq(0, 9, length.out=21),
  ymax=seq(0, 9, length.out=21),
  ymin=0
)

g2 <- ggplot(reslist_all2 %>% filter(pair=="italic(L.~minor)~vs~italic(L.~perpusilla)")) +
  geom_ribbon(data=plotdata1, aes(x, ymin=ymin, ymax=ymax), fill="gray70", alpha=0.5) +
  geom_ribbon(data=plotdata2, aes(x, ymin=ymin, ymax=ymax), fill="gray90", alpha=0.5) +
  geom_hdr(aes(log(1/rho), log(fdiff), fill = after_stat(probs)), alpha=1)  +
  annotate("text", x=0, y=7, label="Competitive\nexclusion", size=2) +
  annotate("text", x=6.5, y=1, label="Coexistence", size=2) +
  annotate("text", x=-6.5, y=1, label="Priority\neffect", size=2) +
  geom_abline(intercept=0, slope=1) +
  geom_abline(intercept=0, slope=-1) +
  scale_x_continuous("Niche difference") +
  scale_fill_viridis_d("Probability") +
  scale_y_continuous("Fitness difference", expand=c(0, 0), limits=c(0, NA)) +
  coord_cartesian(xlim=c(-8, 8)) +
  scale_color_viridis_d(end=0.8) +
  facet_grid(condition~treatment) +
  ggtitle("L. perpusilla vs S. polyrhiza") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "right"
  )

g3 <- ggplot(reslist_all2 %>% filter(pair=="italic(L.~minor)~vs~italic(S.~polyrhiza)")) +
  geom_ribbon(data=plotdata1, aes(x, ymin=ymin, ymax=ymax), fill="gray70", alpha=0.5) +
  geom_ribbon(data=plotdata2, aes(x, ymin=ymin, ymax=ymax), fill="gray90", alpha=0.5) +
  geom_hdr(aes(log(1/rho), log(fdiff), fill = after_stat(probs)), alpha=1)  +
  annotate("text", x=0, y=7, label="Competitive\nexclusion", size=2) +
  annotate("text", x=6.5, y=1, label="Coexistence", size=2) +
  annotate("text", x=-6.5, y=1, label="Priority\neffect", size=2) +
  geom_abline(intercept=0, slope=1) +
  geom_abline(intercept=0, slope=-1) +
  scale_x_continuous("Niche difference") +
  scale_fill_viridis_d("Probability") +
  scale_y_continuous("Fitness difference", expand=c(0, 0), limits=c(0, NA)) +
  coord_cartesian(xlim=c(-8, 8)) +
  scale_color_viridis_d(end=0.8) +
  facet_grid(condition~treatment) +
  ggtitle("L. perpusilla vs S. polyrhiza") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "right"
  )

g4 <- ggplot(reslist_all2 %>% filter(pair=="italic(L.~perpusilla)~vs~italic(S.~polyrhiza)")) +
  geom_ribbon(data=plotdata1, aes(x, ymin=ymin, ymax=ymax), fill="gray70", alpha=0.5) +
  geom_ribbon(data=plotdata2, aes(x, ymin=ymin, ymax=ymax), fill="gray90", alpha=0.5) +
  geom_hdr(aes(log(1/rho), log(fdiff), fill = after_stat(probs)), alpha=1)  +
  annotate("text", x=0, y=7, label="Competitive\nexclusion", size=2) +
  annotate("text", x=6.5, y=1, label="Coexistence", size=2) +
  annotate("text", x=-6.5, y=1, label="Priority\neffect", size=2) +
  geom_abline(intercept=0, slope=1) +
  geom_abline(intercept=0, slope=-1) +
  scale_x_continuous("Niche difference") +
  scale_fill_viridis_d("Probability") +
  scale_y_continuous("Fitness difference", expand=c(0, 0), limits=c(0, NA)) +
  coord_cartesian(xlim=c(-8, 8)) +
  scale_color_viridis_d(end=0.8) +
  facet_grid(condition~treatment) +
  ggtitle("L. perpusilla vs S. polyrhiza") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "right"
  )

gcomb <- ggarrange(g2, g3, g4, nrow=3, draw=FALSE)

ggsave("figure_stanfit_default_niche2.pdf", gcomb, width=18, height=12)
ggsave("figure_stanfit_default_niche2.png", gcomb, width=18, height=12)

g5_v1 <- ggplot(reslist_all2 %>% filter(pair=="italic(L.~perpusilla)~vs~italic(S.~polyrhiza)") %>%
         filter(treatment %in% c("no herbivore", "ETA", "ETAxTLExRNE"))) +
  geom_ribbon(data=plotdata1, aes(x, ymin=ymin, ymax=ymax), fill="gray70", alpha=0.5) +
  geom_ribbon(data=plotdata2, aes(x, ymin=ymin, ymax=ymax), fill="gray90", alpha=0.5) +
  geom_point(aes(log(1/rho), log(fdiff), col=outcome), size=0.5) +
  annotate("text", x=0, y=7, label="Competitive\nexclusion", size=2) +
  annotate("text", x=5, y=1, label="Coexistence", size=2) +
  annotate("text", x=-6, y=1, label="Priority\neffect", size=2) +
  geom_abline(intercept=0, slope=1) +
  geom_abline(intercept=0, slope=-1) +
  scale_x_continuous("Niche difference") +
  scale_y_continuous("Fitness difference", expand=c(0, 0), limits=c(0, NA)) +
  coord_cartesian(xlim=c(-8, 8)) +
  scale_color_viridis_d(end=0.8) +
  facet_grid(condition~treatment) +
  ggtitle("L. perpusilla vs S. polyrhiza") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

ggsave("figure_stanfit_default_niche_subset_v1.pdf", g5_v1, width=8, height=4)

reslist_all3 <- reslist_all2 %>% filter(pair=="italic(L.~perpusilla)~vs~italic(S.~polyrhiza)") %>%
  filter(treatment %in% c("no herbivore", "ETA", "ETAxTLExRNE"))

g5_v2 <- ggplot(reslist_all2 %>% filter(pair=="italic(L.~perpusilla)~vs~italic(S.~polyrhiza)") %>%
         filter(treatment %in% c("no herbivore", "ETA", "ETAxTLExRNE"))) +
  geom_ribbon(data=plotdata1, aes(x, ymin=ymin, ymax=ymax), fill="gray70", alpha=0.5) +
  geom_ribbon(data=plotdata2, aes(x, ymin=ymin, ymax=ymax), fill="gray90", alpha=0.5) +
  geom_hdr(aes(log(1/rho), log(fdiff), fill = after_stat(probs)), alpha=1) +
  annotate("text", x=0, y=7, label="Competitive\nexclusion", size=2) +
  annotate("text", x=6.5, y=1, label="Coexistence", size=2) +
  annotate("text", x=-6.5, y=1, label="Priority\neffect", size=2) +
  geom_abline(intercept=0, slope=1) +
  geom_abline(intercept=0, slope=-1) +
  scale_x_continuous("Niche difference") +
  scale_fill_viridis_d("Probability") +
  scale_y_continuous("Fitness difference", expand=c(0, 0), limits=c(0, NA)) +
  coord_cartesian(xlim=c(-8, 8)) +
  scale_color_viridis_d(end=0.8) +
  facet_grid(condition~treatment) +
  ggtitle("L. perpusilla vs S. polyrhiza") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "right"
  )

ggsave("figure_stanfit_default_niche_subset_v2.pdf", g5_v2, width=8, height=4)
