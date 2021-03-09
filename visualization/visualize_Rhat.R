"
Visualize the resulting estimation of the R series of performed tests (from
`test.R`).
"
library(RJSONIO)
library(ggplot2)

test = fromJSON("./results/202103081847.json", nullValue = NaN, null = NaN)
# Read the results of a random simulation
simulation_results = sample(test$results$simulations, 1)[[1]]
# Build dataframe for ggplot
R = test$parameters$R
plot_df = data.frame(R)
plot_df$t = 1:length(R)
plot_df$`Epiestim R` = unlist(simulation_results$epiestim$R_hat)
plot_df$`Epiestim+Imm R` = unlist(simulation_results$epiestim_immigration$R_hat)
plot_df$`EKF` = unlist(simulation_results$ekf$R_hat)
plot_df$`EKF2` = unlist(simulation_results$ekf2$R_hat[1:50])
ggplot(plot_df[30:50,], aes(x=t)) +
  geom_line(aes(y=R, color="Ground Truth"), linetype="dashed") +
  geom_line(aes(y=`Epiestim R`, color="Epiestim")) +
  geom_line(aes(y=`Epiestim+Imm R`, color="Epiestim+Imm")) +
  geom_line(aes(y=`EKF`, color="EKF")) +
  geom_line(aes(y=`EKF2`, color="EKF2")) +
  scale_colour_manual(
    values = c(
      "Ground Truth" = "black",
      "Epiestim" = "blue",
      "Epiestim+Imm" = "darkgreen",
      "EKF" = "deeppink1",
      "EKF2" = "red"
      ),
    )