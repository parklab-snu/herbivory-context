library(rstan)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(here)
here::i_am("figure_stanfit/figure_stanfit_default_all_post.R")

source(here::here("script", "script_diversity.R"))

load(here::here("stanfit", "stanfit_default_all_refactor.rda"))

fitdata_group <- fitdata %>%
  group_by(group) %>%
  filter(1:n()==1) %>%
  mutate(
    treatment=gsub(" ", "~", treatment)
  )

posterior_array <- as.array(stanfit_default_all_refactor)
species_names <- c("italic(L.~minor)", "italic(L.~perpusilla)", "italic(S.~polyrhiza)")

treatment_levels <- c("no~herbivore", 
                      "ETA", "TLE", "RNE", 
                      "ETAxRNE", "ETAxTLE", "TLExRNE",
                      "ETAxTLExRNE")

posterior_chain <- data.frame(
  draw=rep(1:dim(posterior_array)[1],
           times=dim(posterior_array)[2]*dim(posterior_array)[3]),
  chain=rep(rep(1:dim(posterior_array)[2],
                each=dim(posterior_array)[1]),
            times=dim(posterior_array)[3]),
  parameter=rep(dimnames(posterior_array)$parameters,
                each=dim(posterior_array)[1]*dim(posterior_array)[2]),
  value=as.vector(posterior_array)
)

posterior_lambda <- posterior_chain %>%
  filter(grepl("^lambda\\[", parameter)) %>%
  mutate(
    index=gsub(".*\\[", "", parameter),
    index=gsub("\\]", "", index),
    species_index=as.numeric(gsub(",.*", "", index)),
    group=as.numeric(gsub(".*,", "", index))
  ) %>%
  mutate(
    species=species_names[species_index],
    treatment=fitdata_group$treatment[group],
    condition=fitdata_group$condition[group],
    treatment=factor(treatment, levels=treatment_levels)
  )

g_lambda_hyper <- ggplot(posterior_lambda %>% filter(condition=="hyper")) +
  geom_density(aes(value, col=as.factor(chain)), linewidth=0.4) +
  facet_grid(species~treatment, labeller = label_parsed, scales="free_y") +
  scale_x_continuous("Intrinsic growth rate (1/day)") +
  scale_y_continuous("Posterior density") +
  scale_color_discrete("Posterior chains") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_blank(),
    legend.position = "top"
  )

g_lambda_meso <- ggplot(posterior_lambda %>% filter(condition=="meso")) +
  geom_density(aes(value, col=as.factor(chain)), linewidth=0.4) +
  facet_grid(species~treatment, labeller = label_parsed, scales="free_y") +
  scale_x_continuous("Intrinsic growth rate (1/day)") +
  scale_y_continuous("Posterior density") +
  scale_color_discrete("Posterior chains") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_blank(),
    legend.position = "top"
  )

posterior_alpha <- posterior_chain %>%
  filter(grepl("^alpha\\[", parameter)) %>%
  mutate(
    index=gsub(".*\\[", "", parameter),
    index=gsub("\\]", "", index),
    species_index=as.numeric(gsub(",.*", "", index)),
    xvar_index=gsub("^[^,]*,", "", index),
    xvar_index=as.numeric(gsub(",.*", "", xvar_index)),
    group=as.numeric(gsub(".*,", "", index)),
    species=species_names[species_index],
    xvar=species_names[xvar_index],
    treatment=fitdata_group$treatment[group],
    condition=fitdata_group$condition[group],
    treatment=factor(treatment, levels=treatment_levels),
    text=paste0("atop(", xvar, ",", species, ")")
  ) 

g_alpha_hyper <- ggplot(posterior_alpha %>% filter(condition=="hyper")) +
  geom_density(aes(value, col=as.factor(chain)), linewidth=0.4) +
  facet_grid(text~treatment, labeller = label_parsed, scales="free_y") +
  scale_x_continuous("Competition coefficient") +
  scale_y_continuous("Posterior density") +
  scale_color_discrete("Posterior chains") +
  coord_cartesian(xlim=c(0, 3)) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_blank(),
    legend.position = "top"
  )

g_alpha_meso <- ggplot(posterior_alpha %>% filter(condition=="meso")) +
  geom_density(aes(value, col=as.factor(chain)), linewidth=0.4) +
  facet_grid(text~treatment, labeller = label_parsed, scales="free_y") +
  scale_x_continuous("Competition coefficient") +
  scale_y_continuous("Posterior density") +
  scale_color_discrete("Posterior chains") +
  coord_cartesian(xlim=c(0, 3)) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_blank(),
    legend.position = "top"
  )

ggsave(
  here::here("figure_stanfit", "figure_lambda_hyper.png"),
  g_lambda_hyper,
  width=8,
  height=6
)
ggsave(
  here::here("figure_stanfit", "figure_lambda_meso.png"),
  g_lambda_meso,
  width=8,
  height=6
)
ggsave(
  here::here("figure_stanfit", "figure_alpha_hyper.png"),
  g_alpha_hyper,
  width=8,
  height=8
)
ggsave(
  here::here("figure_stanfit", "figure_alpha_meso.png"),
  g_alpha_meso,
  width=8,
  height=6
)
