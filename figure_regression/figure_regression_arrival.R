library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(lme4)
library(ggtext)
library(flextable)
library(here)
here::i_am("figure_regression/figure_regression_arrival.R")

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
      RNE=grepl("RNE", treatment)
    ) %>%
    filter(ETA==FALSE),
  data_diversity_arrival %>%
    filter(day==5) %>%
    mutate(
      ETA=grepl("ETA", treatment),
      TLE=grepl("TLE", treatment),
      RNE=grepl("RNE", treatment)
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

fitdata_comb <- bind_rows(
  data_processed_all %>%
    group_by(condition, treatment, replication, species) %>%
    summarise(
      area3 = area[day == 3],
      area5 = area[day == 5]
    ) %>%
    ungroup %>%
    mutate(
      r = (log(area5) - log(area3)) / 2,
      ETA = grepl("ETA", treatment),
      TLE = grepl("TLE", treatment),
      RNE = grepl("RNE", treatment),
      nutrient = factor(
        condition,
        levels = c("meso", "hyper"),
        labels = c("low-nutrient", "high-nutrient")
      )
    ) %>%
    filter(ETA == FALSE),
  data_processed_arrival %>%
    group_by(condition, treatment, replication, species) %>%
    summarise(
      area3 = area[day == 3],
      area5 = area[day == 5]
    ) %>%
    ungroup %>%
    mutate(
      r = (log(area5) - log(area3)) / 2,
      ETA = grepl("ETA", treatment),
      TLE = grepl("TLE", treatment),
      RNE = grepl("RNE", treatment),
      nutrient = factor(
        condition,
        levels = c("meso", "hyper"),
        labels = c("low-nutrient", "high-nutrient")
      )
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

g1 <- ggplot(filter(data_reg_coef_r, species == "lm")) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(est, treatment), size = 2) +
  geom_errorbarh(
    aes(xmin = lwr, xmax = upr, y = treatment),
    height = 0.2
  ) +
  scale_x_continuous(
    name = expression("Effect on growth rate" ~ (day^{-1})),
    limits=c(-0.6, 0.6), expand=c(0, 0)
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
    limits=c(-0.6, 0.6), expand=c(0, 0)
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
    limits=c(-0.6, 0.6), expand=c(0, 0)
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
    "figure_regression_arrival.pdf"
  ),
  gcomb,
  width=12,
  height=4
)

ggsave(
  here::here(
    "figure_regression",
    "figure_regression_arrival.png"
  ),
  gcomb,
  width=12,
  height=4
)
