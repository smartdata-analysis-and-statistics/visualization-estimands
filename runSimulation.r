source("R/visualisation.r")
source("R/sim.r")
source("R/ps.r")

################################################################################
#### Simulation Study
################################################################################
# Run the simulation study
run_sim(seed = 944, n = 120000, caliper_width = 0.001,  dir = "Datasets/")

# Generate all figures from the simulation study. Figures will be placed in the /Figures/ directory.
generate_figures(scenario_num = 1)
generate_figures(scenario_num = 2)

