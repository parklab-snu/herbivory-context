library(rstan)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(here)
here::i_am("figure_prior/figure_prior.R")

nsim <- 20000

reslist <- vector('list', nsim)
for (i in 1:nsim) {
  lambda <- abs(rnorm(2, 1, 5)-1)+1
  alphavec <- abs(rnorm(4, 0, 5))
  
  alphamat <- matrix(alphavec, 2, 2)
  
  rho <- sqrt(alphamat[1,2]*alphamat[2,1]/(alphamat[1,1]*alphamat[2,2]))
  fdiff <- (lambda[2]-1)/(lambda[1]-1) * sqrt(alphamat[1,2] * alphamat[1,1]/(alphamat[2,1]*alphamat[2,2]))
  fdiff <- pmax(fdiff, 1/fdiff)
  
  reslist[[i]] <- data.frame(
    rho=rho,
    fdiff=fdiff
  )
}

resdata <- reslist %>%
  bind_rows

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

g1 <- ggplot(resdata) +
  geom_ribbon(data=plotdata1, aes(x, ymin=ymin, ymax=ymax), fill="gray70", alpha=0.5) +
  geom_ribbon(data=plotdata2, aes(x, ymin=ymin, ymax=ymax), fill="gray90", alpha=0.5) +
  geom_hdr(aes(log(1/rho), log(fdiff), fill = after_stat(probs)), alpha=1) +
  geom_abline(intercept=0, slope=1) +
  geom_abline(intercept=0, slope=-1) +
  scale_x_continuous("Niche difference") +
  scale_fill_viridis_d("Probability") +
  scale_y_continuous("Fitness difference", expand=c(0, 0), limits=c(0, NA)) +
  coord_cartesian(xlim=c(-8, 8)) +
  scale_color_viridis_d(end=0.8) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "right"
  )

ggsave(here::here("figure_prior", "figure_prior.png"), g1, width=3, height=2)
