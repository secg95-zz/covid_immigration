library(RJSONIO)

sample_test_simulations = function(use_latest_test = TRUE) {
  "
  Loads the results of a single simulation, randomly drawn from the selected
  test results file.
  
  Parameters
  ----------
  use_latest_test : logical
    If TRUE, the simulation is drawn from the  latest generated test results
    file. Else, the user is prompted for the file name over the command line.
  
  Returns
  -------
  simulation_results : list
    List with keys `I`, `XI`, `epiestim`, `ekf`, .. (one for every estimation
    algorithm).
  "
  # Load results file (contains multiple simulations)
  result_filenames = list.files("results/")
  result_filenames = result_filenames[grep("*.json$", result_filenames)]
  if (use_latest_test) {
    selected_filename = result_filenames[length(result_filenames)]
  } else {
    print("Available result files:")
    print(result_filenames)
    selected_filename = readline(prompt="Enter file name (e.g. 202103170954.json): ")
    selected_filename = gsub("\"", "", selected_filename)
  }
  test = fromJSON(paste0("results/", selected_filename), nullValue = NaN, null=NA)
  # Read the results of a random simulation
  simulation_results = sample(test$results$simulations, 1)[[1]]
  return(simulation_results)
}
