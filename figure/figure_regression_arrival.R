library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(lme4)
library(ggtext)
library(flextable)
library(here)
here::i_am("figure/figure_regression_main.R")

source(
  here::here(
    "script",
    "script_diversity.R"
  )
)

data_diversity_all_reg <- bind_rows(
  data_diversity_all %>%
    filter(day==5) %>%
    mutate(
      ETA=grepl("ETA", treatment),
      TLE=grepl("TLE", treatment),
      RNE=grepl("RNE", treatment),
      arrival=FALSE
    )
) %>%
  mutate(
    nutrient = factor(
      condition,
      levels = c("meso", "hyper"),
      labels = c("low-nutrient", "high-nutrient")
    )
  )

lfit_all <- lm(diversity~nutrient*(ETA + TLE + RNE)^2,
               data=data_diversity_all_reg)

coef_order <- c("ETATRUE", "TLETRUE", "RNETRUE", 
                "ETATRUE:TLETRUE", "ETATRUE:RNETRUE", "TLETRUE:RNETRUE", 
                "nutrienthigh-nutrient", "nutrienthigh-nutrient:ETATRUE", 
                "nutrienthigh-nutrient:TLETRUE", "nutrienthigh-nutrient:RNETRUE", 
                "nutrienthigh-nutrient:ETATRUE:TLETRUE", "nutrienthigh-nutrient:ETATRUE:RNETRUE", 
                "nutrienthigh-nutrient:TLETRUE:RNETRUE")

data_reg_coef <- bind_rows(
  data.frame(
    treatment=c("ETA", "TLE", "RNE", 
                "ETA:TLE", "ETA:RNE","TLE:RNE",
                "high-nutrient", "high-nutrient:ETA",
                "high-nutrient:TLE", "high-nutrient:RNE",
                "high-nutrient:ETA:TLE", "high-nutrient:ETA:RNE", 
                "high-nutrient:TLE:RNE"),
    est=coef(lfit_all)[coef_order],
    lwr=confint(lfit_all)[coef_order,1],
    upr=confint(lfit_all)[coef_order,2]
  )
) %>%
  mutate(
    treatment=factor(treatment,
                     levels=c("ETA", "TLE", "RNE", 
                              "ETA:TLE", "ETA:RNE","TLE:RNE",
                              "high-nutrient", "high-nutrient:ETA",
                              "high-nutrient:TLE", "high-nutrient:RNE",
                              "high-nutrient:ETA:TLE", "high-nutrient:ETA:RNE", 
                              "high-nutrient:TLE:RNE"))
  )

g4 <- ggplot(data_reg_coef) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(aes(est, treatment), size=2) +
  geom_errorbarh(aes(xmin=lwr, xmax=upr, y=treatment),
                 height=0.2) +
  scale_x_continuous("Effect on Shannon diversity") +
  scale_y_discrete("Treatment",limits=rev) +
  theme(
    panel.border = element_rect(linewidth=1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

fitdata_comb <- data_processed_all %>%
  group_by(condition, treatment, replication, species) %>%
  summarise(
    area0 = area[day == 0],
    area5 = area[day == 5]
  ) %>%
  ungroup %>%
  mutate(
    r = (log(area5) - log(area0)) / 5,
    ETA = grepl("ETA", treatment),
    TLE = grepl("TLE", treatment),
    RNE = grepl("RNE", treatment),
    nutrient = factor(
      condition,
      levels = c("meso", "hyper"),
      labels = c("low-nutrient", "high-nutrient")
    )
  )

growth_formula <- r ~ nutrient*(ETA + TLE + RNE)^2

lfit_all_r_lm <- lm(
  growth_formula,
  data = filter(fitdata_comb, species == "lm")
)

lfit_all_r_sp <- lm(
  growth_formula,
  data = filter(fitdata_comb, species == "sp")
)

lfit_all_r_lp <- lm(
  growth_formula,
  data = filter(fitdata_comb, species == "lp")
)

growth_labels <- levels(data_reg_coef$treatment)

extract_growth <- function(model, species_code) {
  ci <- confint(model)
  ss <- summary(model)
  
  data.frame(
    treatment = factor(growth_labels, levels = growth_labels),
    est = unname(coef(model)[coef_order]),
    lwr = unname(ci[coef_order, 1]),
    upr = unname(ci[coef_order, 2]),
    p=unname(ss$coefficients[coef_order,4]),
    species = species_code
  )
}

data_reg_coef_r <- bind_rows(
  extract_growth(lfit_all_r_lm, "lm"),
  extract_growth(lfit_all_r_lp, "lp"),
  extract_growth(lfit_all_r_sp, "sp")
)

data_reg_coef_r_pretty <- data_reg_coef_r %>%
  mutate(
    est=signif(est,2),
    lwr=signif(lwr,2),
    upr=signif(upr,2),
    ci=paste0(lwr, " - ", upr),
    p=signif(p,2)
  ) %>%
  select(species, treatment, est, ci, p) %>%
  mutate(
    p=ifelse(p < 1e-3, "p < 0.001", p),
    p=ifelse(p < 1e-2, "p < 0.01", p),
    p=ifelse(p < 0.05, "p < 0.05", p)
  )

flextable(data_reg_coef_r_pretty)

g1 <- ggplot(filter(data_reg_coef_r, species == "lm")) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(est, treatment), size = 2) +
  geom_errorbarh(
    aes(xmin = lwr, xmax = upr, y = treatment),
    height = 0.2
  ) +
  scale_x_continuous(
    name = expression("Effect on growth rate" ~ (day^{-1})),
    limits=c(-1.25, 1), expand=c(0, 0)
  ) +
  scale_y_discrete("Treatment", limits = rev) +
  ggtitle(expression(italic(L.~minor))) +
  theme(
    panel.border = element_rect(linewidth = 1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  )

g2 <- ggplot(filter(data_reg_coef_r, species == "lp")) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(est, treatment), size = 2) +
  geom_errorbarh(
    aes(xmin = lwr, xmax = upr, y = treatment),
    height = 0.2
  ) +
  scale_x_continuous(
    name = expression("Effect on growth rate" ~ (day^{-1})),
    limits=c(-1.25, 1), expand=c(0, 0)
  ) +
  scale_y_discrete("Treatment", limits = rev) +
  ggtitle(expression(italic(L.~perpusilla))) +
  theme(
    panel.border = element_rect(linewidth = 1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

g3 <- ggplot(filter(data_reg_coef_r, species == "sp")) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(est, treatment), size = 2) +
  geom_errorbarh(
    aes(xmin = lwr, xmax = upr, y = treatment),
    height = 0.2
  ) +
  scale_x_continuous(
    name = expression("Effect on growth rate" ~ (day^{-1})),
    limits=c(-1.25, 1), expand=c(0, 0)
  ) +
  scale_y_discrete("Treatment", limits = rev) +
  ggtitle(expression(italic(S.~polyrhiza))) +
  theme(
    panel.border = element_rect(linewidth = 1),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

gcomb <- ggarrange(g1, g2, g3, g4, nrow=1, labels=c("A", "B", "C", "D"))

ggsave(
  here::here(
    "figure",
    "figure_regression_main.pdf"
  ),
  gcomb,
  width=12,
  height=4
)

ggsave(
  here::here(
    "figure",
    "figure_regression_main.png"
  ),
  gcomb,
  width=12,
  height=4
)

high_nutrient_reference <- function(model) {
  model_data <- model.frame(model)
  
  model_data$nutrient <- relevel(
    model_data$nutrient,
    ref = "high-nutrient"
  )
  
  lm(formula(model), data = model_data)
}

lfit_all_r_lm_high <- high_nutrient_reference(lfit_all_r_lm)
lfit_all_r_lp_high <- high_nutrient_reference(lfit_all_r_lp)
lfit_all_r_sp_high <- high_nutrient_reference(lfit_all_r_sp)

extract_eta <- function(model, species_code, condition_code) {
  ci <- confint(model, parm = "ETATRUE")
  
  data.frame(
    species = species_code,
    condition = condition_code,
    est = unname(coef(model)["ETATRUE"]),
    lwr = unname(ci[1, 1]),
    upr = unname(ci[1, 2])
  )
}

data_reg_coef_r_eta <- bind_rows(
  extract_eta(lfit_all_r_lm, "lm", "mesotrophic"),
  extract_eta(lfit_all_r_lm_high, "lm", "eutrophic"),
  extract_eta(lfit_all_r_lp, "lp", "mesotrophic"),
  extract_eta(lfit_all_r_lp_high, "lp", "eutrophic"),
  extract_eta(lfit_all_r_sp, "sp", "mesotrophic"),
  extract_eta(lfit_all_r_sp_high, "sp", "eutrophic")
)

data_processed_eta_pref <- read.csv(
  here::here("data_processed", "data_processed_eta_pref.csv")
) %>%
  rename(
    species = key,
    condition = treatment
  ) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("meso", "hyper"),
      labels = c("mesotrophic", "eutrophic")
    )
  )

data_reg_coef_r_eta_pref <- data_processed_eta_pref %>%
  left_join(
    data_reg_coef_r_eta,
    by = c("species", "condition")
  ) %>%
  mutate(
    species = factor(
      species,
      levels = c("lm", "lp", "sp"),
      labels = c("L. minor", "L. perpusilla", "S. polyrhiza")
    ),
    plot_scale = ifelse(condition == "eutrophic", 30, 10),
    effect_plot = -est * plot_scale,
    effect_lwr = -upr * plot_scale,
    effect_upr = -lwr * plot_scale
  )


g5 <- ggplot(data_reg_coef_r_eta_pref %>% filter(type=="nest", condition=="eutrophic")) +
  geom_bar(aes(species, preference), stat="identity", color="black", fill=NA) +
  geom_point(aes(species, -est*30), col="red") +
  geom_errorbar(aes(species, ymin=-upr*30, ymax=-lwr*30), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Nesting preference", expand=c(0, 0),
                     sec.axis = sec_axis(transform=~-./30,
                                         "Mean effect on growth rate")) +
  ggtitle("Nesting preference, high-nutrient") +
  coord_cartesian(ylim=c(0, 34)) +
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
  geom_errorbar(aes(species, ymin=-upr*10, ymax=-lwr*10), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Nesting preference", expand=c(0, 0),
                     sec.axis = sec_axis(transform=~-./10,
                                         "Mean effect on growth rate")) +
  coord_cartesian(ylim=c(0, 34)) +
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
  geom_errorbar(aes(species, ymin=-upr*30, ymax=-lwr*30), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Feeding preference", expand=c(0, 0),
                     sec.axis = sec_axis(transform=~-./30,
                                         "Mean effect on growth rate")) +
  ggtitle("Feeding preference, high-nutrient") +
  coord_cartesian(ylim=c(0, 34)) +
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
  geom_errorbar(aes(species, ymin=-upr*10, ymax=-lwr*10), width=0, col="red") +
  scale_x_discrete("Species") +
  scale_y_continuous("Feeding preference", expand=c(0, 0),
                     sec.axis = sec_axis(transform=~-./10,
                                         "Mean effect on growth rate")) +
  ggtitle("Feeding preference, low-nutrient") +
  coord_cartesian(ylim=c(0, 34)) +
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

ggsave(
  here::here(
    "figure",
    "figure_regression_main_pref.pdf"
  ),
  gpref,
  width=6*1.1,
  height=4*1.1
)

ggsave(
  here::here(
    "figure",
    "figure_regression_main_pref.png"
  ),
  gpref,
  width=6*1.1,
  height=4*1.1
)
