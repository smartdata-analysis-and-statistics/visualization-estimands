library(mvtnorm) # generate correlated normal data

run_scenario <- function(pars) {

  dat <- gendat(pars)
  
  # Derive propensity score
  dat$psAB <- dat$lpAB <- dat$psAC <- dat$lpAC <- NA
  dat$psAB[dat$trtC == 0] <- derive_ps(subset(dat, trtC == 0), trt_var = "trtA")
  dat$lpAB[dat$trtC == 0] <- logit(dat$psAB[dat$trtC == 0])
  dat$psAC[dat$trtB == 0] <- derive_ps(subset(dat, trtB == 0), trt_var = "trtA")
  dat$lpAC[dat$trtB == 0] <- logit(dat$psAC[dat$trtB == 0])
  dat$psBC[dat$trtA == 0] <- derive_ps(subset(dat, trtA == 0), trt_var = "trtB")
  dat$lpBC[dat$trtA == 0] <- logit(dat$psBC[dat$trtA == 0])
  
  # Derive caliper
  caliper_AB <- pars$caliper * sd(logit(subset(dat, trtC == 0)$psAB)) # Caliper used for AvsB comparisons
  caliper_AC <- pars$caliper * sd(logit(subset(dat, trtB == 0)$psAC)) # Caliper used for AvsC comparisons
  caliper_BC <- pars$caliper * sd(logit(subset(dat, trtA == 0)$psBC)) # Caliper used for BvsC comparisons
  
  # Save output
  sim <- list(pars = pars,
              caliper_AB = caliper_AB, # Caliper for AB comparisons
              caliper_AC = caliper_AC,
              caliper_BC = caliper_BC,
              data = dat
  )
  
  if (!is.null(pars$save.dir)) {
    if (!dir.exists(pars$save.dir)) {
      dir.create(pars$save.dir)
    }
    save(sim, file = paste(pars$save.dir, "simdat.rda", sep = ""))
    
    # Remove existing simulation results
    fname_results <- paste(pars$save.dir, "results.rda", sep = "")
    if (file.exists(fname_results)) {
      file.remove(fname_results)
    }
  }
  
  ATT <- derive_ATT("A", "B", replace = TRUE, data = dat, pars = pars)
  ATT <- derive_ATT("B", "A", replace = TRUE, data = dat, pars = pars)
  ATT <- derive_ATT("A", "C", replace = TRUE, data = dat, pars = pars)
  ATT <- derive_ATT("C", "A", replace = TRUE, data = dat, pars = pars)
  ATT <- derive_ATT("B", "C", replace = TRUE, data = dat, pars = pars)
  ATT <- derive_ATT("C", "B", replace = TRUE, data = dat, pars = pars)
  
  # Sampling without replacement
  ATT <- derive_ATT("A", "B", replace = FALSE, data = dat, pars = pars)
  ATT <- derive_ATT("B", "A", replace = FALSE, data = dat, pars = pars)
  ATT <- derive_ATT("A", "C", replace = FALSE, data = dat, pars = pars)
  ATT <- derive_ATT("C", "A", replace = FALSE, data = dat, pars = pars)
  ATT <- derive_ATT("B", "C", replace = FALSE, data = dat, pars = pars)
  ATT <- derive_ATT("C", "B", replace = FALSE, data = dat, pars = pars)
}

run_sim <- function(seed, # Random seed
                    n, # Total sample size
                    caliper_width, # caliper width, specified as multiplication factor of the SD of the LP 
                    dir = NULL  # Directory to save all results
                    ) 
  {
  
  # Check if the folder exists
  if (!missing(dir) && !dir.exists(dir)) {
    dir.create(dir)
  }
  if (caliper_width <= 0) {
    stop("Invalid caliper width")
  }
  if (n <= 1) {
    stop("Not enough patients to simulate")
  }
  
  # Scenario 1
  pars_1 <- retrieve_configuration(HTE = FALSE, 
                                   caliper = caliper_width,
                                   n = n,
                                   seed = seed,
                                   dir = paste(dir, "Scenario 1/", sep = ""))
  run_scenario(pars = pars_1) 
  
  # Scenario 2
  pars_2 <- retrieve_configuration(HTE = TRUE, 
                                   caliper = caliper_width,
                                   n = n,
                                   seed = seed,
                                   dir = paste(dir, "Scenario 2/", sep = ""))
  run_scenario(pars = pars_2) 
 
}


# This function retrieves the ATT objects for a particular estimand
retrieve_ATT_files <- function(population, replace, fdir) {
  # Retrieve relevant ATT samples
  lbl_sampling <- ifelse(replace, "with_replacement", "without_replacement")
  fsuffix <- paste("_", population, "_", lbl_sampling, ".rda$", sep = "")
  list.files(path = fdir, pattern = fsuffix)
}

# Configure the super population
retrieve_configuration <- function(
  HTE = FALSE,
  caliper, #default caliper width
  n, # Sample size
  seed, # Random seed
  dir = NULL) # Directory to save the results 
{
  
  pars <- list()
  
  pars$HTE <- HTE
  pars$caliper <- caliper
  pars$n <- n # Sample size
  pars$seed <- seed
  pars$save.dir <- dir
  
  
  pars$mean.trtA <- c(x1 = 1, x2 = -1)
  pars$mean.trtB <- c(x1 = 0, x2 = 1) 
  pars$mean.trtC <- c(x1 = -1, x2 = -1)
  
  pars$sigma.trtA <-  matrix(c(2, 0, 0, 2), ncol = 2)
  pars$sigma.trtB <-  matrix(c(2, 0, 0, 2), ncol = 2)
  pars$sigma.trtC <-  matrix(c(2, 0, 0, 2), ncol = 2)
  colnames(pars$sigma.trtA) <- colnames(pars$sigma.trtB) <- colnames(pars$sigma.trtC) <- c("x1", "x2")
  rownames(pars$sigma.trtA) <- rownames(pars$sigma.trtB) <- rownames(pars$sigma.trtC) <- c("x1", "x2")
  
  pars$beta0 <- 0 # Baseline outcome risk
  pars$beta1 <- 1 # Confounder effect of x1
  pars$beta2 <- 0.25 # Confounder effect of x2
  pars$beta3 <- -0.5 # Treatment effect of B versus A
  pars$beta4 <- -0.25 # Treatment effect of C versus A
  pars$beta5 <- 0 # Effect modification between B and x1
  pars$beta6 <- 0 # Effect modification between C and x1
  pars$sigma <- 1 # SD of the residual error
  
  if (HTE) {
    pars$beta5 <- 0.25
    pars$beta6 <- 0.125
  } 
   
  return(pars)
}

# Function to generate data with 2 continuous outcomes, 3 treatments and 2 outcomes (with and without HTE)
gendat <- function(pars){
  
  set.seed(pars$seed)
  
  data <- data.frame(ID = seq(pars$n),
                     x1 = NA, 
                     x2 = NA, 
                     treatment = NA)
  
  if (pars$n %% 3 != 0) {
    warning("Sample size cannot be divided by 3")
  }
  
  # Generate 3 treatments with equal sample size
  data$treatment <- factor(rep(c("A", "B", "C"), each = pars$n/3))
  
  
  # Draw x1 x2 from bivariate Normal distribution
  data[data$treatment == "A", c("x1", "x2")] <- rmvnorm(sum(data$treatment == "A"),
                                                        mean = pars$mean.trtA,
                                                        sigma = pars$sigma.trtA)
  
  data[data$treatment == "B", c("x1", "x2")] <- rmvnorm(sum(data$treatment == "B"),
                                                        mean = pars$mean.trtB,
                                                        sigma = pars$sigma.trtB)
  
  data[data$treatment == "C", c("x1", "x2")] <- rmvnorm(sum(data$treatment == "C"),
                                                        mean = pars$mean.trtC,
                                                        sigma = pars$sigma.trtC)

  lp <- pars$beta0 + 
    pars$beta1 * data$x1 + 
    pars$beta2 * data$x2 + 
    pars$beta3 * (data$treatment == "B") + 
    pars$beta4 * (data$treatment == "C") +
    pars$beta5 * (data$treatment == "B") * data$x1 +
    pars$beta6 * (data$treatment == "C") * data$x1
  
  data$y <- rnorm(n = pars$n, mean = lp, sd = pars$sigma)
  
  data <- data %>% mutate(trtA = ifelse(treatment == "A", 1, 0),
                          trtB = ifelse(treatment == "B", 1, 0),
                          trtC = ifelse(treatment == "C", 1, 0))
  
  return(data)
}
