"
Visualize the internal parameters of a fitted Kalman filter from performed tests
(from `test.R`).
"
library(EpiEstim)
library(ggplot2)
library(cowplot)

source("sample_test_simulations.R")

# Read the results of a random simulation
simulation_results = sample_test_simulations(use_latest_test = TRUE)
# Construct dataframe with series to be visualized
R = test$parameters$R
plot_df = data.frame(R)
plot_df$t = 1:length(R)
plot_df$R_hat = unlist(simulation_results$ekf$R_hat)
plot_df$`exp(a_1)` = exp(unlist(lapply(simulation_results$ekf$a, function(a_t) a_t[1])))
plot_df$M_hat = exp(unlist(lapply(simulation_results$ekf$alpha, function(a_t) a_t[2])))
plot_df$`exp(a_2)` = exp(unlist(lapply(simulation_results$ekf$a, function(a_t) a_t[2])))
plot_df$`exp(alpha_3)` = exp(unlist(lapply(simulation_results$ekf$alpha, function(a_t) a_t[3])))
M = overall_infectivity(simulation_results$XI, c(0, test$discrete_si))
plot_df$M = tail(M, length(R))
plot_df$I = simulation_results$I
plot_df$XI = tail(simulation_results$XI, length(R))
# Plot 1: Incidence series
I_plot = ggplot(plot_df, aes(x=t)) +
  geom_line(aes(y=I, color="I")) +
  geom_line(aes(y=XI, color="XI")) +
  scale_colour_manual(
    values = c(
      "I" = "black",
      "XI" = "blue"
    )
  )
# Plot 2: R_hat = exp(alpha_1)
R_hat_plot = ggplot(plot_df, aes(x=t)) +
  geom_line(aes(y=R, linetype="Ground Truth")) +
  geom_line(aes(y=R_hat, linetype="Smoothed")) +
  geom_line(aes(y=`exp(a_1)`, linetype="Unsmoothed")) +
  scale_linetype_manual(
    values = c(
      "Ground Truth" = 1,
      "Smoothed" = 2,
      "Unsmoothed" = 3
    )
  )
# Plot 3: M_hat = exp(alpha_2)
M_hat_plot = ggplot(plot_df, aes(x=t)) +
  geom_line(aes(y=M, linetype="Ground Truth"), color="blue") +
  geom_line(aes(y=M_hat, linetype="M_hat"), color="blue") +
  geom_line(aes(y=`exp(a_2)`, linetype="exp(a_2)"), color="blue") +
  scale_linetype_manual(
    values = c(
      "Ground Truth" = 1,
      "M_hat" = 2,
      "exp(a_2)" = 3
    )
  )
# Plot 4: R_hat derivative = exp(alpha_3)
alpha_3_plot = ggplot(plot_df, aes(x=t)) +
  geom_line(aes(y=`exp(alpha_3)`), color="purple", linetype=2)
# Draw all plots in a single frame
plot_grid(I_plot, R_hat_plot, M_hat_plot, alpha_3_plot, ncol=1, axis="lr", align="v")