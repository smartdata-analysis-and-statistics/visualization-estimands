library(sandwich)
library(lmtest)
library(MatchIt)

derive_ps <- function(data, trt_var) {
  fit <- glm(paste(trt_var, " ~ x1 + x2"), data = data, family = binomial)
  predict(fit, type = "response")
}


apply_matching <- function(target_population, # A, B, or C
                           comparison, # AB, AC, or BC
                           caliper, # caliper width specified in SD units
                           replace, # TRUE or FALSE 
                           data) {

  if (comparison == "AB") {
    popdat <- subset(data, trtC == 0)
  } else if (comparison == "AC") {
    popdat <- subset(data, trtB == 0)
  } else if (comparison == "BC") {
    popdat <- subset(data, trtA == 0)
  } else {
    stop("Invalid treatment comparison")
  }

  fmla <- as.formula(paste("trt", target_population, "~ x1 + x2", sep = ""))
  
  fit <- matchit(fmla, 
          data = popdat, 
          family = binomial, 
          method = "nearest", 
          caliper = caliper, 
          std.caliper = TRUE,
          distance = "logit", 
          replace = replace)
  
  # Add information on estimand
  fit$target_comparison <- comparison
  fit$target_population <- target_population
  
  return(fit)
}



################################################################################
# Derive ATT for a particular comparison using a specific method
################################################################################
derive_ATT <- function(treatment_target, # Treatment used for target population ("A", "B", or "C")
                       treatment_comparison, # Comparative treatment ("A", "B", or "C")
                       replace, # Sampling with or without replacement?
                       data, 
                       pars) {
  
  # Set random seed
  set.seed(pars$seed)
  
  if (!(treatment_target %in% c("A", "B", "C"))) {
    stop("The active treatment should be A, B or C!")
  }
  if (!(treatment_comparison %in% c("A", "B", "C"))) {
    stop("The comparative treatment should be A, B or C!")
  }
  if (treatment_target == treatment_comparison) {
    stop("The analysis should involve two different treatments!")
  }
  
  treatments <- sort(c(treatment_target, treatment_comparison))
  comparison <- paste(treatments, collapse = "")
  
  
  message(paste("Estimating ATT", treatment_target, " in population ", comparison, 
                " using 1-NN matching ", ifelse(replace, "with", "without"), 
                " replacement", sep = "" ))
  
  # Derive PSM structure and reference values
  if (comparison == "AB") {
    fmla <- as.formula("y ~ trtA")
    lbl_trt <- "trtA"
    ref_beta_main <- -pars$beta3
    ref_beta_hte <- -pars$beta5
  } else if (comparison == "AC") {
    fmla <- as.formula("y ~ trtA")
    lbl_trt <- "trtA"
    ref_beta_main <- -pars$beta4
    ref_beta_hte <- -pars$beta6
  } else if (comparison == "BC") {
    fmla <- as.formula("y ~ trtB")
    lbl_trt <- "trtB"
    ref_beta_main <- pars$beta3 - pars$beta4
    ref_beta_hte <- pars$beta5 - pars$beta6
  } else {
    stop("Invalid treatment comparison!")
  }
  
  ##############################################################################
  # Derive matched sample targeting individuals treated by treatment_target
  ##############################################################################
  fit <- apply_matching(treatment_target, 
                        comparison = comparison, 
                        caliper = pars$caliper, 
                        replace = replace, 
                        data = data)
  
  if (!is.null(pars$save.dir)) {
    # Save the matched sample
    save(fit, file = paste(pars$save.dir, "ATT", 
                           treatment_target, "_", 
                           comparison, "_", 
                           ifelse(replace, "with", "without"), 
                           "_replacement.rda", sep = ""))
  }
  
  # Summarize the matched sample
  xsum <- summary(fit)
  
  ##############################################################################
  # Estimate the ATT in individuals treated by treatment_target
  ##############################################################################
  mdat <- match.data(fit)
  mdat$ID <- rownames(mdat)
  
  if (!fit$info$replace) {
    #Derive cluster-robust standard error
    fitEst <- lm(fmla, data = mdat, weights = weights)
    fitVar <- vcovCL(fitEst, cluster = ~ subclass) 
  } else {
    mdat <- get_matches(fit)
    
    # We need to adjust for ID since some individuals may be present multiple times
    fitEst <- lm(fmla, data = mdat, weights = weights) 
    fitVar <- vcovCL(fitEst, cluster = ~ subclass + id) 
  }
  
  result <- data.frame(
    Comparison = paste(treatments[1], "versus", treatments[2]),
    Estimand = paste("ATT", treatment_target, sep = ""),
    Replace = replace, 
    Reference = ref_beta_main + ref_beta_hte * (pars[[paste("mean.trt", treatment_target, sep = "")]])["x1"],
    Estimate = coef(fitEst)[lbl_trt], 
    SE = sqrt(diag(fitVar)[lbl_trt]),
    SS = rowSums(xsum$nn)["Matched"], # Sample size
    ESS = rowSums(xsum$nn)["Matched (ESS)"] # Effective sample size
  )
  
  if (!is.null(pars$save.dir)) {
    fname <- paste(pars$save.dir, "results.rda", sep = "")
    results <- data.frame()
    if (file.exists(fname)) {
      load(paste(pars$save.dir, "results.rda", sep = ""))
      results <- rbind(results, result)
    } else {
      results <- result
    }
    save(results, file = fname)
  }

  return(result)
}


################################################################################
# PROPENSITY SCORE MATCHING FOR MSPATHS
################################################################################


apply_matching_MS <- function(target_population, comparison,
                              treatments, 
                              psvars, # Covariates to use for matching
                              caliper, # caliper width specified in SD units
                              replace, # TRUE or FALSE 
                              data) {
  
  popdat <- subset(data, dmt_nm %in% treatments & visit == 1)
  
  # Create an indicator for the target population
  popdat$target <- ifelse(popdat$dmt_nm == target_population, 1,0)
  
  fmla <- as.formula(paste("target ~ ", paste(psvars, collapse = "+"), sep = ""))
  
  fit <- matchit(fmla, 
                 data = popdat, 
                 family = binomial, 
                 method = "nearest", 
                 caliper = caliper, 
                 std.caliper = TRUE,
                 distance = "logit", 
                 replace = replace)
  
  # Add information on estimand
  fit$target_comparison <- comparison
  fit$target_population <- target_population
  
  return(fit)
}



################################################################################
# Derive ATT for a particular comparison using a specific method
################################################################################
derive_ATT_MS <- function(treatment_target, # Treatment used for target population 
                          treatment_comparison, # Comparative treatment 
                          psvars,
                          replace, # Sampling with or without replacement?
                          data,
                          seed = 201,
                          caliper = 0.2,
                          save.dir = NULL) {
  
  # Set random seed
  set.seed(seed)
  
  if (!(treatment_target %in% c("TERI", "DMF", "NAT"))) {
    stop("The active treatment should be a DMT!")
  }
  if (!(treatment_comparison %in% c("TERI", "DMF", "NAT"))) {
    stop("The comparative treatment should be DMT!")
  }
  if (treatment_target == treatment_comparison) {
    stop("The analysis should involve two different treatments!")
  }
  
  x <- c("TERI", "DMF", "NAT")
  tr_vector <- factor(c(treatment_target, treatment_comparison), levels = x, ordered = T)
  
  
  treatments <- sort(tr_vector)
  comparison <- paste(treatments, collapse = "_vs_")
  
  
  message(paste("Estimating the ATT of ", treatment_target, " in population patients receiving ", 
                paste(treatments, collapse = " or "), 
                " using 1-NN matching ", ifelse(replace, "with", "without"), 
                " replacement", sep = "" ))
  
  ##############################################################################
  # Derive matched sample targeting individuals treated by treatment_target
  ##############################################################################
  fit <- apply_matching_MS(treatment_target, comparison = comparison,
                        treatments = treatments,
                        psvars = psvars,
                        caliper = caliper, 
                        replace = replace, 
                        data = data)
  
if (!is.null(save.dir)) {
  # Save the matched sample
  save(fit, file = paste(save.dir, "ATT", 
                         treatment_target, "_", 
                         comparison, "_", 
                         ifelse(replace, "with", "without"), 
                         "_replacement.rda", sep = ""))
}

  # Summarize the matched sample
  xsum <- summary(fit)
  
  ##############################################################################
  # Estimate the ATT in individuals treated by treatment_target
  ##############################################################################
  
  mdat <- get_matches(fit)
  mdat$ID <- rownames(mdat)
  
  # Derive PSM structure
  if (comparison == "TERI_vs_DMF") {
    mdat$target <- ifelse(mdat$dmt_nm == "TERI", 1, 0)
    fmla <- as.formula("z_mdtdh.365 ~ target")
    lbl_trt <- "target"
  } else if (comparison == "TERI_vs_NAT") {
    mdat$target <- ifelse(mdat$dmt_nm == "TERI", 1, 0)
    fmla <- as.formula("z_mdtdh.365 ~ target")
    lbl_trt <- "target"
  } else if (comparison == "DMF_vs_NAT") {
    mdat$target <- ifelse(mdat$dmt_nm == "DMF", 1, 0)
    fmla <- as.formula("z_mdtdh.365 ~ target")
    lbl_trt <- "target"
  } else {
    stop("Invalid treatment comparison!")
  }
  
  
  if (!fit$info$replace) {
    #Derive cluster-robust standard error
    fitEst <- lm(fmla, data = mdat, weights = weights)
    fitVar <- vcovCL(fitEst, cluster = ~ subclass) 
  } else {
    # We need to adjust for ID since some individuals may be present multiple times
    fitEst <- lm(fmla, data = mdat, weights = weights) 
    fitVar <- vcovCL(fitEst, cluster = ~ subclass + id) 
    fitConf <- coefci(fitEst, vcov = vcovCL, cluster = ~subclass + id)
    fitlow <- fitConf[2,1]
    fithigh <- fitConf[2,2]
  }
  
  result <- data.frame(
    Comparison = paste(treatments[1], "versus", treatments[2]),
    Estimand = paste("ATT", treatment_target, sep = " "),
    Replace = replace, 
    Estimate = coef(fitEst)[lbl_trt], 
    SE = sqrt(diag(fitVar)[lbl_trt]),
    low = fitlow,
    high = fithigh,
    SS = rowSums(xsum$nn)["Matched"], # Sample size
    ESS = rowSums(xsum$nn)["Matched (ESS)"] # Effective sample size
  )
  
  
  return(result)
}
