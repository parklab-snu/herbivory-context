library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(ggtext)
library(egg)
library(here)
here::i_am("figure/figure_data.R")

source(
  here::here(
    "script",
    "script_diversity.R"
  )
)

data_processed_all2 <- data_processed_all %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore", 
                              "ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "ETAxTLExRNE")),
    species=factor(species,
                   levels=c("lm", "lp", "sp"),
                   labels=c("italic(L.~minor)", "italic(L.~perpusilla)",
                            "italic(S.~polyrhiza)"))
  )

g1 <- ggplot(data_processed_all2) +
  geom_path(aes(factor(day), area, group=interaction(replication, condition),
                col=condition, lty=condition)) +
  geom_point(aes(factor(day), area, group=interaction(replication, condition),
                 col=condition, shape=condition)) +
  scale_x_discrete("Day") +
  scale_y_continuous("Area") +
  scale_color_manual(values=c("purple", "orange")) +
  facet_grid(species~treatment, 
             labeller = labeller(
               species = label_parsed,
               treatment = label_value
  )) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    strip.background = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_markdown()
  )

data_processed_arrival2 <- data_processed_arrival %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore-ETA", 
                              "RNE-ETA", "TLE-ETA",
                              "TLExRNE-ETA")),
    species=factor(species,
                   levels=c("lm", "lp", "sp"),
                   labels=c("italic(L.~minor)", "italic(L.~perpusilla)",
                            "italic(S.~polyrhiza)"))
  )

g2 <- ggplot(data_processed_arrival2) +
  geom_path(aes(factor(day), area, group=interaction(replication, condition),
                col=condition, lty=condition)) +
  geom_point(aes(factor(day), area, group=interaction(replication, condition),
                 col=condition, shape=condition)) +
  scale_x_discrete("Day") +
  scale_y_continuous("Area") +
  scale_color_manual(values=c("purple", "orange")) +
  facet_grid(species~treatment, 
             labeller = labeller(
               species = label_parsed,
               treatment = label_value
             )) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    strip.background = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_markdown()
  )

data_diversity_all2 <- data_diversity_all %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore", 
                              "ETA", "TLE", "RNE", 
                              "ETAxRNE", "ETAxTLE", "TLExRNE",
                              "ETAxTLExRNE"))
  )

g3 <- ggplot(data_diversity_all2) +
  geom_path(aes(factor(day), diversity, group=interaction(replication, condition),
                col=condition, lty=condition)) +
  geom_point(aes(factor(day), diversity, group=interaction(replication, condition),
                 col=condition, shape=condition))  +
  scale_x_discrete("Day") +
  scale_y_continuous("Shannon diversity") +
  scale_color_manual(values=c("purple", "orange")) +
  facet_grid(.~treatment) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    strip.background = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    strip.text.x = element_markdown()
  )

data_diversity_arrival2 <- data_diversity_arrival %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("no herbivore-ETA", 
                              "RNE-ETA", "TLE-ETA",
                              "TLExRNE-ETA"))
  )

g4 <- ggplot(data_diversity_arrival2) +
  geom_path(aes(factor(day), diversity, group=interaction(replication, condition),
                col=condition, lty=condition)) +
  geom_point(aes(factor(day), diversity, group=interaction(replication, condition),
                 col=condition, shape=condition))  +
  scale_x_discrete("Day") +
  scale_y_continuous("Shannon diversity") +
  scale_color_manual(values=c("purple", "orange")) +
  facet_grid(.~treatment) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=0.7),
    strip.background = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    strip.text.x = element_markdown()
  )

gcomb <- ggarrange(g1, g2, g3, g4, nrow=2, heights=c(3, 1),
          widths=c(8, 4),
          labels=c("A", "B", "C", "D"))

ggsave(
  here::here(
    "figure",
    "figure_data.pdf"
  ),
  gcomb,
  width=15,
  height=6
)

ggsave(
  here::here(
    "figure",
    "figure_data.png"
  ),
  gcomb,
  width=15,
  height=6
)
