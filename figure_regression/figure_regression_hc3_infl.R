library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(lme4)
library(ggtext)
library(car)
library(flextable)
library(here)
here::i_am("figure_regression/figure_regression_infl.R")

source(
  here::here(
    "script",
    "script_diversity.R"
  )
)

hc3_results <- function(model) {
  b <- coef(model)
  V <- car::hccm(model, type = "hc3")
  se <- sqrt(diag(V))[names(b)]
  df <- df.residual(model)
  critical <- qt(0.975, df)
  
  data.frame(
    est = unname(b),
    lwr = unname(b - critical * se),
    upr = unname(b + critical * se),
    p = unname(
      2 * pt(abs(b / se), df = df, lower.tail = FALSE)
    ),
    row.names = names(b)
  )
}

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
               data=data_diversity_all_reg[-108,])

hc3_div <- hc3_results(lfit_all)

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
    est = hc3_div[coef_order, "est"],
    lwr = hc3_div[coef_order, "lwr"],
    upr = hc3_div[coef_order, "upr"]
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
  data = filter(fitdata_comb, species == "lm")[-99,]
)

lfit_all_r_sp <- lm(
  growth_formula,
  data = filter(fitdata_comb, species == "sp")[-95,]
)

lfit_all_r_lp <- lm(
  growth_formula,
  data = filter(fitdata_comb, species == "lp")[-90,]
)

growth_labels <- levels(data_reg_coef$treatment)

extract_growth <- function(model, species_code) {
  hc <- hc3_results(model)[coef_order, , drop = FALSE]
  
  data.frame(
    treatment = factor(growth_labels, levels = growth_labels),
    est = hc$est,
    lwr = hc$lwr,
    upr = hc$upr,
    p = hc$p,
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
    "figure_regression_hc3_infl.pdf"
  ),
  gcomb,
  width=12,
  height=4
)

ggsave(
  here::here(
    "figure_regression",
    "figure_regression_hc3_infl.png"
  ),
  gcomb,
  width=12,
  height=4
)
