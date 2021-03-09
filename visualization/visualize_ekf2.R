"
Visualize the internal parameters of a fitted Kalman filter from performed tests
(from `test.R`).
"
library(RJSONIO)
library(EpiEstim)
library(ggplot2)
library(cowplot)

test = fromJSON("./results/202103081305.json", nullValue = NaN, null = NaN)
# Read the results of a random simulation
simulation_results = sample(x=test$results$simulations, size=1)[[1]]
# Construct dataframe with series to be visualized
R = c(test$parameters$R)
plot_df = data.frame(R)
plot_df$t = 1:nrow(plot_df)
plot_df$R_hat = exp(unlist(lapply(simulation_results$ekf2$alpha_hat, function(x) x[1])))
plot_df$`exp(a_1)` = exp(unlist(lapply(simulation_results$ekf2$a, function(x) x[1])))
plot_df$M_hat = exp(unlist(lapply(simulation_results$ekf2$alpha_hat, function(x) x[2])))
plot_df$`exp(a_2)` = exp(unlist(lapply(simulation_results$ekf2$a, function(x) x[2])))
plot_df$`exp(alpha_3)` = exp(unlist(lapply(simulation_results$ekf2$alpha_hat, function(x) x[3])))
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
  geom_line(aes(y=M_hat, linetype="Smoothed"), color="blue") +
  geom_line(aes(y=`exp(a_2)`, linetype="Unsmoothed"), color="blue") +
  scale_linetype_manual(
    values = c(
      "Ground Truth" = 1,
      "Smoothed" = 2,
      "Unsmoothed" = 3
    )
  )
# Plot 4: R_hat derivative = exp(alpha_3)
alpha_3_plot = ggplot(plot_df, aes(x=t)) +
  geom_line(aes(y=`exp(alpha_3)`), color="purple", linetype=2)
# Draw all plots in a single frame
plot_grid(I_plot, R_hat_plot, M_hat_plot, alpha_3_plot, ncol=1, axis="lr", align="v")