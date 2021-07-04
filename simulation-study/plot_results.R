library(RJSONIO)
library(ggplot2)

CONFIG_FILE = "./simulation-study/config.json"
RESULTS_FILE = "./simulation-study/results/202106271257.json"
MAPE_PLOT_PATH = "./simulation-study/plots/mean-MAPE.png"
R_PLOT_PATH = "./simulation-study/plots/scenario-%s-%s.png"
ALGORITHM_COLOR_MAP = list(
  "ekf"="blue",
  "epiestim"="deeppink1",
  "epiestim_immigration"='darkgreen'
)
config = fromJSON(CONFIG_FILE)
results = fromJSON(RESULTS_FILE)

# Plot mean MAPE across simulations for each scenario
full_mape_table = data.frame(matrix(ncol=3, nrow=0))
names(full_mape_table) = c("Scenario", "Algorithm", "MAPE")
for (scenario_results in results) {
  scenario_mape_table = data.frame(
    "Scenario" = scenario_results$scenario_id,
    "Algorithm" = names(scenario_results$mean_mape),
    "MAPE" = scenario_results$mean_mape,
    row.names=NULL
  )
  full_mape_table = rbind(full_mape_table, scenario_mape_table)
}
mape_plot = ggplot(full_mape_table) +
  geom_col(position="dodge", aes(x=Scenario, y=MAPE, fill=Algorithm)) +
  scale_fill_manual(values = ALGORITHM_COLOR_MAP) +
  xlab("Scenario") +
  ylab("Mean MAPE(R)") +
  theme(legend.title = element_blank(), legend.text=element_text(size=15),
        axis.title=element_text(size=15), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12), legend.position="bottom"
        )
png(MAPE_PLOT_PATH, pointsize=15, width=630, height=350)
print(mape_plot)
dev.off()

# Plot R estimations for all R simulations for each scenario
for (scenario_results in results) {
  # Locate scenario configuration
  for (scenario_config in config$scenarios) {
    if (scenario_config$id == scenario_results$scenario_id) break
  }
  # Plot each algorithm's R estimations
  for (algorithm in names(ALGORITHM_COLOR_MAP)) {
    # Build table with R estimations for this (scenario, algorithm)
    algorithm_R_table = data.frame(
      t=1:config$T,
      label="R",
      value=scenario_config$R
    )
    R_hat_idx = 1
    for (epidemic_results in scenario_results$epidemic_results) {
      R_hat = epidemic_results$R_hat[[algorithm]]
      R_hat[sapply(R_hat, is.null)] = NA
      R_hat = data.frame(
        t=1:config$T,
        label=paste0("R_hat", R_hat_idx),
        value=unlist(R_hat)
      )
      R_hat_idx = R_hat_idx + 1
      algorithm_R_table = rbind(algorithm_R_table, R_hat)
    }
    # Plot R estimations for this (scenario, algorithm)
    R_hat_plot = ggplot(algorithm_R_table) +
      geom_line(data=algorithm_R_table, aes(x=t, y=value, color=label, alpha=label)) +
      ylab("R(t)") +
      scale_colour_manual(
        name="",
        values=c("R"="black", "R_hat"=ALGORITHM_COLOR_MAP[[algorithm]]),
        labels=c("R"="Truth", "R_hat"=algorithm),
        na.value=ALGORITHM_COLOR_MAP[[algorithm]]
      ) +
      scale_alpha_manual(
        name="",
        values=c("R"=1, "R_hat"=0.3),
        guide="none",
        na.value=0.3
      ) +
      theme(legend.position="bottom")
    png(
      sprintf(R_PLOT_PATH, scenario_results$scenario_id, algorithm),
      pointsize=15, width=630, height=350
    )
    print(R_hat_plot)
    dev.off()
  }
}