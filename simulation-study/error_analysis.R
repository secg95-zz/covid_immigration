"
Visualize the internal parameters of a fitted Kalman filter from performed tests
(from `simulation-study/run.R`).
"
library(RJSONIO)
library(EpiEstim)
library(ggplot2)
library(cowplot)

RESULT_ID = readline("Result set ID (looks like `202201050426`): ")
CONFIG_FILE = "./simulation-study/config.json"
RESULTS_FILE = file.path("./simulation-study/results", paste0(RESULT_ID, ".json"))
OUTPUT_DIR = file.path("./simulation-study/error-analysis", RESULT_ID)
config = fromJSON(CONFIG_FILE)
results = fromJSON(RESULTS_FILE)
dir.create(OUTPUT_DIR, recursive=TRUE)

# Read the results of a random simulation
esimation_results = results[[1]]$epidemic_results[[1]]$raw_output
epidemic_results = results[[1]]$epidemic_results[[1]]
M = unlist(epidemic_results$M)

# Construct dataframe with series to be visualized
R = config$scenarios[[1]]$R
plot_df = data.frame(R)
plot_df$t = 1:nrow(plot_df)
# R_hat
R_hat = esimation_results$ekf$R_hat
R_hat[sapply(R_hat, is.null)] = NA
plot_df$R_hat = unlist(R_hat)
# exp(alpha_1)
alpha_1 = lapply(esimation_results$ekf$alpha, function(x) x[1])
alpha_1[sapply(alpha_1, is.null)] = NA
plot_df$`exp(alpha_1)` = exp(unlist(alpha_1))
# M_hat
M_hat = lapply(esimation_results$ekf$alpha, function(x) x[2])
M_hat[sapply(M_hat, is.null)] = NA
plot_df$`M_hat` = exp(unlist(M_hat))
# exp(alpha_3)
alpha_3 = lapply(esimation_results$ekf$alpha, function(x) x[3])
alpha_3[sapply(alpha_3, is.null)] = NA
plot_df$`exp(alpha_3)` = exp(unlist(alpha_3))
# Ground truth
plot_df$M = tail(M, length(R))
plot_df$I = epidemic_results$I
plot_df$XI = tail(epidemic_results$XI, length(R))
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
  ylim(0, max(plot_df$R * 2))
  scale_linetype_manual(
    values = c(
      "Ground Truth" = 1,
      "Smoothed" = 2
    )
  )
# Plot 3: M_hat = exp(alpha_2)
M_hat_plot = ggplot(plot_df, aes(x=t)) +
  geom_line(aes(y=M, linetype="Ground Truth"), color="blue") +
  geom_line(aes(y=M_hat, linetype="Smoothed"), color="blue") +
  scale_linetype_manual(
    values = c(
      "Ground Truth" = 1,
      "Smoothed" = 2
    )
  )
# Plot 4: R_hat derivative = exp(alpha_3)
alpha_3_plot = ggplot(plot_df, aes(x=t)) +
  geom_line(aes(y=`exp(alpha_3)`), color="purple", linetype=2)
# Draw all plots in a single frame
plot_grid(I_plot, R_hat_plot, M_hat_plot, alpha_3_plot, ncol=1, axis="lr", align="v")