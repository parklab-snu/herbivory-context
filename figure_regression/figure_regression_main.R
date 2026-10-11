library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(lme4)
library(ggtext)
library(flextable)
library(here)
here::i_am("figure_regression/figure_regression_main.R")

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
    "figure_regression",
    "figure_regression_main.pdf"
  ),
  gcomb,
  width=12,
  height=4
)

ggsave(
  here::here(
    "figure_regression",
    "figure_regression_main.png"
  ),
  gcomb,
  width=12,
  height=4
)

plot_diagnostics <- function(model, model_name) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
  
  plot(model, which = 1:4, ask = FALSE, sub.caption = "", id.n = 3)
  
  mtext(model_name, side = 3, outer = TRUE, line = 1, cex = 1.1)
}

png(file = here::here("figure_regression", "figure_regression_diag_lm.png"), width = 8,height = 7,units = "in", res = 300)
plot_diagnostics(lfit_all_r_lm, "L. minor")
dev.off()

png(file = here::here("figure_regression", "figure_regression_diag_lp.png"), width = 8,height = 7,units = "in", res = 300)
plot_diagnostics(lfit_all_r_lp, "L. perpusilla")
dev.off()

png(file = here::here("figure_regression", "figure_regression_diag_sp.png"), width = 8,height = 7,units = "in", res = 300)
plot_diagnostics(lfit_all_r_sp, "S. polyrhiza")
dev.off()

png(file = here::here("figure_regression", "figure_regression_diag_div.png"), width = 8,height = 7,units = "in", res = 300)
plot_diagnostics(lfit_all, "Shannon diversity")
dev.off()
