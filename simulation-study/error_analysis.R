# Visualize the internal parameters of a fitted Kalman filter from performed tests
# (from `simulation-study/run.R`).

library(RJSONIO)
library(EpiEstim)
library(comprehenr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(argparse)

CONFIG_FILE = "./simulation-study/config.json"
COLORS = brewer.pal(name = "Dark2", n=8)
ALGORITHM_COLOR_MAP = list(
  "ekf"=COLORS[1],
  "mukf"=COLORS[2]
)

parser = ArgumentParser(description="Plot error analysis for a result set")
parser$add_argument("result_id", help="Result set ID (looks like `202201050426`)")
cl_args = parser$parse_args()
result_id = cl_args$result_id
results_file = file.path("./simulation-study/results", paste0(result_id, ".json"))
output_dir = file.path("./simulation-study/error-analysis", result_id)
config = fromJSON(CONFIG_FILE)
results = fromJSON(results_file)

for (scenario_idx in 1:length(config$scenarios)) {
  scenario_results = results[[scenario_idx]]
  scenario_config = config$scenarios[[scenario_idx]]
  scenario_output_dir = file.path(output_dir, paste0("scenario", scenario_idx))
  dir.create(scenario_output_dir, recursive=TRUE)
  for (algorithm in c("ekf", "mukf")) {
    # Construct dataframe with series to be visualized
    plot_df = data.frame(
      "ID" = "Ground Truth",
      "t" = 1:config$T,
      "I" = rep(NA, config$T),
      "XI" = rep(NA, config$T),
      "M" = rep(NA, config$T),
      "R_hat" = scenario_config$R,
      "M_hat" = rep(NA, config$T)
    )
    epidemic_idx = 1
    for (epidemic_result in scenario_results$epidemic_results) {
      # Format estimation result
      estimation_result = epidemic_result$raw_output[[algorithm]]
      R_hat = estimation_result$R_hat
      R_hat[sapply(R_hat, is.null)] = NA
      R_hat = unlist(R_hat)
      alpha_2 = lapply(estimation_result$alpha, function(x) x[2])
      alpha_2[sapply(alpha_2, is.null)] = NA
      M_hat = exp(unlist(alpha_2))
      M = epidemic_result$M
      M[sapply(M, is.null)] = NA
      M = unlist(M)
      epidemic_df = data.frame(
        "ID" = paste("Epidemic", epidemic_idx),
        "t" = 1:config$T,
        "I" = epidemic_result$I,
        "XI" = tail(epidemic_result$XI, config$T),
        "M" = tail(M, config$T),
        "R_hat" = R_hat,
        "M_hat" = M_hat
      )
      # Append epidemic results to full results
      plot_df = rbind(plot_df, epidemic_df)
      epidemic_idx = epidemic_idx + 1
    }
    
    # Plot 1: Incidence series
    I_plot = ggplot(plot_df) +
      geom_line(data=plot_df, aes(x=t, y=I, color="I", alpha="I", group=ID)) +
      geom_line(data=plot_df, aes(x=t, y=XI, color="XI", alpha="XI", group=ID)) +
      ylab("Incidence") +
      scale_colour_manual(
        name="Label",
        values=c("I"="black", "XI"="blue"),
        labels=c("I"="Local", "XI"="Immigrant")
      ) +
      scale_alpha_manual(
        values=c("I"=0.3, "XI"=0.3),
        guide="none"
      )
    
    # Plot 2: R_hat = exp(alpha_1)
    plot_df$is_ground_truth = plot_df$ID == "Ground Truth"
    plot_df$is_ground_truth = sapply(plot_df$is_ground_truth, toString)
    R_plot = ggplot(plot_df) +
      geom_line(data=plot_df, aes(x=t, y=R_hat, color=is_ground_truth, alpha=is_ground_truth, group=ID)) +
      ylab("R") +
      ylim(min(scenario_config$R) / 4, max(scenario_config$R) * 2) +
      scale_colour_manual(
        name="Label",
        values=c("TRUE"="black", "FALSE"=ALGORITHM_COLOR_MAP[[algorithm]]),
        labels=c("TRUE"="Ground Truth", "FALSE"=algorithm)
      ) +
      scale_alpha_manual(
        values=c("TRUE"=1, "FALSE"=0.3),
        guide="none"
      )
    
    # Plot 3: M_hat = exp(alpha_2)
    M_plot = ggplot(plot_df) +
      geom_line(data=plot_df, aes(x=t, y=M, color="Ground Truth", alpha="Ground Truth", group=ID)) +
      geom_line(data=plot_df, aes(x=t, y=M_hat, color="Estimation", alpha="Estimation", group=ID)) +
      ylab("M") +
      ylim(min(plot_df$M, na.rm=TRUE) / 4, max(plot_df$M, na.rm=TRUE) * 1.3) +
      scale_colour_manual(
        name="Label",
        values=c("Ground Truth"="black", "Estimation"=ALGORITHM_COLOR_MAP[[algorithm]]),
        labels=c("Ground Truth"="Ground Truth", "Estimation"=algorithm)
      ) +
      scale_alpha_manual(
        values=c("Ground Truth"=0.3, "Estimation"=0.3),
        guide="none"
      )
    
    # Draw all plots in a single frame
    stacked_plot_path = file.path(
      scenario_output_dir,
      paste0(scenario_config$label, "-", algorithm, ".png")
    )
    stacked_plot = plot_grid(I_plot, R_plot, M_plot, ncol=1, axis="lr", align="v")
    png(stacked_plot_path, pointsize=15, width=500, height=600)
    print(stacked_plot)
    dev.off()
  }
}