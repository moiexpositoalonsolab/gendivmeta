DEBUG=TRUE
DEBUG=FALSE
if(DEBUG) message("RUNNING IN DEBUG MODE!")

#### Functions ####

#####*********************************************************************######
##### Utility functions ####

# getsummary<-function(type="A"){
#   d %>% 
#     dplyr::filter(Diversity.Type==type) %>% 
#     mutate(diff=(Recent.Diversity-Early.Diversity)/Early.Diversity) ->
#     # dplyr::select(diff) ->
#     mysummary
#   return(mysummary)
# }
# getplot<-function(type="A"){
#   d %>% 
#     dplyr::filter(Diversity.Type==type) %>% 
#     ggplot(.)+
#     geom_histogram(aes(x=(Recent.Diversity-Early.Diversity)/Early.Diversity))->
#     myplot
#   return(myplot)
# }
compute_p_values <- function(g_values, g_se) {
  # Calculate the t-statistic
  t_values <- g_values / g_se
  
  # Compute two-tailed p-values using the normal distribution
  p_values <- 2 * (1 - pnorm(abs(t_values)))
  
  # Return results as a data frame
  results <- data.frame(Hedges_g = g_values, SE = g_se, t_value = t_values, p_value = p_values)
  
  return(results)
}
compute_p_values_conservative <- function(g_values, g_se, n1, n2) {
  n1<-as.numeric(n1)  
  n2<-as.numeric(n2)  
  # Degrees of freedom
  df = n1 + n2 - 2
  
  # Calculate the t-statistic
  t_values <- g_values / g_se
  
  # Compute two-tailed p-values using the normal distribution
  # p_values <- 2 * (1 - pnorm(abs(t_values)))
  p_values <- 2 * (1 - pt(abs(t_values), df))
  
  # Return results as a data frame
  results <- data.frame(Hedges_g = g_values, SE = g_se, t_value = t_values, p_value = p_values)
  
  return(results)
}
# compute_p_values_lowtest <- function(g_values, g_se) {
#   # Calculate the t-statistic
#   t_values <- g_values / g_se
#   
#   # Compute one-tailed p-values (testing for declines, i.e., g < 0)
#   p_values <- pnorm(t_values)  # One-tailed p-value for g < 0
#   
#   # Return results as a data frame
#   results <- data.frame(Hedges_g = g_values, SE = g_se, t_value = t_values, p_value = p_values)
#   
#   return(results)
# }
# combine_p_values_fisher <- function(p_values) {
#   test_statistic <- -2 * sum(log(p_values))
#   df <- 2 * length(p_values)  # Degrees of freedom
#   p_combined <- 1 - pchisq(test_statistic, df)
#   
#   return(list(statistic = test_statistic, p_value = p_combined))
# }
# compute_diversity_change_old <- function(df) {
#   # Ensure numerical
#   df$Early.Error<-as.numeric(df$Early.Error)
#   df$Recent.Error<-as.numeric(df$Recent.Error)
#   df$Early.Sample.N<-as.numeric(df$Early.Sample.N)
#   df$Recent.Sample.N<-as.numeric(df$Recent.Sample.N)
# 
#   # Ensure proper handling of SD vs. SE
#   convert_to_se <- function(error, sample_size, error_type) {
#     ifelse(error_type == "SD", error / sqrt(sample_size), error)
#   }
# 
#   # Convert Early and Recent Errors to SE if needed
#   df$Early.SE <- convert_to_se(df$Early.Error, df$Early.Sample.N, df$Error.Type)
#   df$Recent.SE <- convert_to_se(df$Recent.Error, df$Recent.Sample.N, df$Error.Type)
# 
#   # Compute percent change
#   df$Percent.Change <- ((df$Recent.Diversity - df$Early.Diversity) / df$Early.Diversity) * 100
# 
#   # Compute SE of the difference
#   df$SE.Diff <- sqrt(df$Early.SE^2 + df$Recent.SE^2)
# 
#   # Compute SE for percent change
#   df$SE.Percent.Change <- (df$SE.Diff / df$Early.Diversity) * 100
# 
#   # Compute 95% Confidence Interval
#   df$Lower.CI <- df$Percent.Change - 1.96 * df$SE.Percent.Change
#   df$Upper.CI <- df$Percent.Change + 1.96 * df$SE.Percent.Change
# 
#   # Compute SD
#   df$SD.Percent.Change <- df$SE.Percent.Change * sqrt(df$Early.Sample.N)  # Assuming df$Early.Sample.N represents the sample size
# 
#   # Compute CV
#   df$CV.Percent.Change <- 100*df$SD.Percent.Change / df$SD.Percent.Change # Assuming df$Early.Sample.N represents the sample size
# 
# 
#   # Return only relevant columns
#   return(df[, c("Recent.Diversity", "Early.Diversity", "Percent.Change", 
#                 "Lower.CI", "Upper.CI","SD.Percent.Change","CV.Percent.Change")])
# }

# updated
compute_diversity_change <- function(df) {
  # Ensure numeric conversion
  df$Early.Error <- as.numeric(df$Early.Error)
  df$Recent.Error <- as.numeric(df$Recent.Error)
  df$Early.Sample.N <- as.numeric(df$Early.Sample.N)
  df$Recent.Sample.N <- as.numeric(df$Recent.Sample.N)
  
  # Convert any error type to SE
  convert_to_se <- function(error, sample_size, error_type) {
    if (error_type == "SD") {
      return(error / sqrt(sample_size))
    } else if (error_type == "SE") {
      return(error)
    } else if (error_type == "95CL") {
      return((error / 2) / 1.96)  # Assuming full width of the 95% CI
    } else {
      return(NA)
    }
  }
  
  # Apply conversion to SE
  df$Early.SE <- mapply(convert_to_se, df$Early.Error, df$Early.Sample.N, 
                        df$Error.Type)
  df$Recent.SE <- mapply(convert_to_se, df$Recent.Error, df$Recent.Sample.N, 
                         df$Error.Type)
  
  # Recover SDs from SEs
  df$Early.SD <- df$Early.SE * sqrt(df$Early.Sample.N)
  df$Recent.SD <- df$Recent.SE * sqrt(df$Recent.Sample.N)
  
  # Compute pooled SD
  df$SD.Pooled <- sqrt(
    ((df$Early.Sample.N - 1) * df$Early.SD^2 +
       (df$Recent.Sample.N - 1) * df$Recent.SD^2) /
      (df$Early.Sample.N + df$Recent.Sample.N - 2)
  )
  
  # Effective sample is useful as harmonic of both
  df$Effective.N <- 2 * df$Early.Sample.N * df$Recent.Sample.N / (df$Early.Sample.N + df$Recent.Sample.N)
  
  # Compute percent change
  df$Percent.Change <- ((df$Recent.Diversity - df$Early.Diversity) / df$Early.Diversity) * 100
  
  # SE of the difference
  df$SE.Diff <- sqrt(df$Early.SE^2 + df$Recent.SE^2)
  
  # SE for percent change
  df$SE.Percent.Change <- (df$SE.Diff / df$Early.Diversity) * 100
  
  # 95% Confidence Intervals
  df$Lower.CI <- df$Percent.Change - 1.96 * df$SE.Percent.Change
  df$Upper.CI <- df$Percent.Change + 1.96 * df$SE.Percent.Change
  
  # Convert pooled SD to percent change scale
  df$SD.Percent.Change <- (df$SD.Pooled / df$Early.Diversity) * 100
  
  # Coefficient of variation of percent change
  df$CV.Percent.Change <- 100 * df$SD.Percent.Change / abs(df$Percent.Change)
  
  # Degrees of freedom
  df$df <- df$Early.Sample.N + df$Recent.Sample.N - 2
  
  # Bias correction factor J
  df$J <- 1 - (3 / (4 * df$df - 1))
  
  # Hedges' g
  df$Hedges.g.recalc <- df$J * (df$Recent.Diversity - df$Early.Diversity) / df$SD.Pooled
  df$Hedges.g.SE.recalc <- sqrt(
    (df$Early.Sample.N + df$Recent.Sample.N) / 
      (df$Early.Sample.N * df$Recent.Sample.N) + 
      (df$Hedges.g.recalc^2) / (2 * (df$Early.Sample.N + df$Recent.Sample.N))
  )
  df$Hedges.g.CI.lower.recalc <- df$Hedges.g.recalc - 1.96 * df$Hedges.g.SE.recalc
  df$Hedges.g.CI.upper.recalc <- df$Hedges.g.recalc + 1.96 * df$Hedges.g.SE.recalc
  
  # Return selected columns
  return(df[, c("Recent.Diversity", "Early.Diversity", "Percent.Change",
                "Lower.CI", "Upper.CI", "SE.Percent.Change", "Effective.N",
                "Hedges.g.recalc","Hedges.g.SE.recalc",
                "Hedges.g.CI.lower.recalc","Hedges.g.CI.upper.recalc",
                "SD.Percent.Change", "CV.Percent.Change")])
}



printHPI<-function(x){
  x=na.omit(x)
  IQ1=quantile(x,prob=c(0.025),na.rm=T) %>% format(digits = 2)
  IQ3=quantile(x,prob=c(0.975),na.rm=T) %>% format(digits = 2)
  paste0(IQ1, " - ", IQ3)
}
printIQR<-function(x){
  x=na.omit(x)
  IQ1=quantile(x,prob=c(0.25),na.rm=T) %>% format(digits = 2)
  IQ3=quantile(x,prob=c(0.75),na.rm=T) %>% format(digits = 2)
  paste0(IQ1, " - ", IQ3)
}
printranges<-function(x){
  x=na.omit(x)
  IQ1=min(x,na.rm=T) %>% format(digits = 2)
  IQ3=max(x,na.rm=T) %>% format(digits = 2)
  paste0(IQ1, " - ", IQ3)
}
printmedian<-function(x){
  x=na.omit(x)
  MED=quantile(x,prob=c(0.5),na.rm=T) %>% format(digits = 2)
  MED
}
printmean<-function(x){
  x=na.omit(x)
  MEAN=mean(x,na.rm=T) %>% format(digits = 2)
  MEAN
}
printmeanSE <- function(x) {
  x <- na.omit(x)
  SE <- sd(x) / sqrt(length(x))
  format(SE, digits = 2)
}
printCI <- function(x) {
  x <- na.omit(x)
  m <- mean(x)
  se <- sd(x) / sqrt(length(x))
  ci_half_width <- 1.96 * se
  sprintf("%.2f%% ± %.2f%% (95%% CI)", m, ci_half_width)
  return(ci_half_width)
}
clean<-function(x){
  format(x,digits=2)
}
#####*********************************************************************######
#### Meta statistics ####
metaztest <- function(g, se) {
  d <- data.frame(g, se)
  d <- na.omit(d)
  
  g <- d$g
  se <- d$se
  
  # inverse-variance weights
  w <- 1 / se^2
  
  # weighted mean
  g_mean <- sum(w * g) / sum(w)
  
  # standard error of the weighted mean
  se_mean <- sqrt(1 / sum(w))
  
  # z-test statistic
  z <- g_mean / se_mean
  
  # p-value (two-tailed)
  p <- 2 * pnorm(-abs(z))
  
  # 95% confidence interval
  ci_half_width <- 1.96 * se_mean
  ci_lower <- g_mean - ci_half_width
  ci_upper <- g_mean + ci_half_width
  
  return(
    list(
      Weighted_mean = g_mean,
      SE_mean = se_mean,
      CI_half_width = ci_half_width,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      z = z,
      p = p,
      summary_string = sprintf("%.2f ± %.2f (95%% CI)", g_mean, ci_half_width)
    )
  )
}
# metamcmc <- function(df, 
#                      variable = "Percent.Change", 
#                      se = "SE.Percent.Change",
#                      random_effects = c("PaperID"),  # vector of random effects
#                      ...) {
#   # Select and clean data
#   df$y <- df[[variable]]
#   df$se <- df[[se]]
#   df <- df[!is.na(df$y) & !is.na(df$se), ]
#   
#   # Compute sampling variances
#   df$var <- df$se^2
#   
#   # Construct random formula
#   random_formula <- as.formula(paste("~", paste(random_effects, collapse = " + ")))
#   
#   # Construct corresponding prior
#   prior <- list(
#     R = list(V = 1, nu = 0.002),
#     G = lapply(seq_along(random_effects), function(i) list(V = 1, nu = 0.002))
#   )
#   names(prior$G) <- paste0("G", seq_along(random_effects))
#   
#   # Run the model
#   model <- MCMCglmm(
#     y ~ 1,
#     random = random_formula,
#     data = df,
#     mev = df$var,
#     verbose = FALSE,
#     prior = prior,
#     ...
#   )
#   
#   return(model)
# }
metamcmc <- function(df, 
                     variable = "Percent.Change", 
                     se = "SE.Percent.Change",
                     fixed_effects = NULL,     # optional fixed effects
                     random_effects = NULL,    # optional random effects
                     prior = NULL,             # optional prior
                     family="gaussian",
                     ...) {
  # Select and clean data
  df$y <- df[[variable]]
  df$se <- df[[se]]
  df <- df[!is.na(df$y) & !is.na(df$se), ]
  
  # Compute sampling variances
  df$var <- df$se^2
  
  # Build fixed effect formula
  fixed_formula <- if (is.null(fixed_effects)) {
    y ~ 1
  } else {
    as.formula(paste("y ~", fixed_effects))
  }
  
  # If no prior provided, construct default
  if (is.null(prior)) {
    prior <- list(R = list(V = 1, nu = 0.002))
    
    if (!is.null(random_effects) && length(random_effects) > 0) {
      prior$G <- lapply(seq_along(random_effects), function(i) list(V = 1, nu = 0.002))
      names(prior$G) <- paste0("G", seq_along(random_effects))
    }
  }
  
  # Run model with or without random effects
  if (!is.null(random_effects) && length(random_effects) > 0) {
    random_formula <- as.formula(paste("~", paste(random_effects, collapse = " + ")))
    
    model <- MCMCglmm(
      fixed = fixed_formula,
      random = random_formula,
      data = df,
      mev = df$var,
      verbose = FALSE,
      prior = prior,
      # family = family,
      ...
    )
  } else {
    model <- MCMCglmm(
      fixed = fixed_formula,
      data = df,
      mev = df$var,
      verbose = FALSE,
      prior = prior,
      # family = family,
      ...
    )
  }
  
  return(model)
}

metamcmclog <- function(df, 
                     variable = "Percent.Change", 
                     se = "SE.Percent.Change",
                     fixed_effects = NULL,
                     random_effects = NULL,
                     prior = NULL,
                     family = "gaussian",
                     log_transform = TRUE,
                     ...) {
  
  # Clean and prepare data
  df$raw_y <- df[[variable]]
  df$raw_se <- df[[se]]
  df <- df[!is.na(df$raw_y) & !is.na(df$raw_se), ]
  
  # Transform: log(1 + percent change / 100) and delta method SE
  if (log_transform) {
    df$y <- log(df$raw_y / 100 + 1)
    df$se <- df$raw_se / (100 + df$raw_y)
  } else {
    df$y <- df$raw_y
    df$se <- df$raw_se
  }
  
  df$var <- df$se^2
  
  # Build fixed effect formula
  fixed_formula <- if (is.null(fixed_effects)) {
    y ~ 1
  } else {
    as.formula(paste("y ~", fixed_effects))
  }
  
  # Build prior if not provided
  if (is.null(prior)) {
    prior <- list(R = list(V = 1, nu = 0.002))
    
    if (!is.null(random_effects) && length(random_effects) > 0) {
      prior$G <- lapply(seq_along(random_effects), function(i) list(V = 1, nu = 0.002))
      names(prior$G) <- paste0("G", seq_along(random_effects))
    }
  }
  
  # Build random effects formula
  if (!is.null(random_effects) && length(random_effects) > 0) {
    random_formula <- as.formula(paste("~", paste(random_effects, collapse = " + ")))
    
    model <- MCMCglmm(
      fixed = fixed_formula,
      random = random_formula,
      data = df,
      mev = df$var,
      verbose = FALSE,
      prior = prior,
      family = family,
      ...
    )
  } else {
    model <- MCMCglmm(
      fixed = fixed_formula,
      data = df,
      mev = df$var,
      verbose = FALSE,
      prior = prior,
      family = family,
      ...
    )
  }
  
  return(model)
}

# meanmcmc<-function(mod){
#   mean(mod$Sol[, "(Intercept)"])
# }

posteriormcmc <- function(mod, use_median = FALSE, hpi_prob = 0.95) {
  require(coda)
  posterior=mod$Sol[, "(Intercept)"]
  
  # DIC
  DIC<-mod$DIC
  
  # Ensure input is an MCMC object
  posterior_mcmc <- as.mcmc(posterior)
  
  # Central estimate
  center <- if (use_median) median(posterior) else mean(posterior)
  
  # Standard error (posterior SD)
  se <- sd(posterior)
  
  # Highest Posterior Interval
  hpi <- HPDinterval(posterior_mcmc, prob = hpi_prob)
  hpi_lower <- hpi[1, "lower"]
  hpi_upper <- hpi[1, "upper"]
  hpi_width <- hpi_upper - hpi_lower
  hpi_half_width <- hpi_width / 2
  
  # Return as list or data.frame
  return(list(
    center = center,
    se = se,
    hpi_lower = hpi_lower,
    hpi_upper = hpi_upper,
    hpi_half_width = hpi_half_width,
    DIC= DIC,
    # summary_string = sprintf("%.2f ± %.2f (%.0f%% HPI) - DIC %.2f", center, hpi_half_width, hpi_prob * 100, DIC)
    summary_string = sprintf("%.2f ± %.2f (%.0f%% HPI)", center, hpi_half_width, hpi_prob * 100)
  ))
}

extract_random_effects_hpi <- function(model, factor_name) {
  # Extract random effect matrix
  rand_effects <- model$Sol[, grep(paste0("^", factor_name), colnames(model$Sol))]
  
  # Extract level names
  level_names <- sub(paste0("^", factor_name), "", colnames(rand_effects))
  level_names <- sub("^\\.", "", level_names)  # clean up leading dot
  
  # Build summary table
  summary_df <- data.frame(
    Level = level_names,
    Mean = apply(rand_effects, 2, mean),
    Lower = apply(rand_effects, 2, function(x) quantile(x, 0.025)),
    Upper = apply(rand_effects, 2, function(x) quantile(x, 0.975)),
    Significant = apply(rand_effects, 2, function(x) {
      q <- quantile(x, c(0.025, 0.975))
      !(q[1] < 0 & q[2] > 0)  # TRUE if 95% HPI excludes 0
    }),
    row.names = NULL
  )
  
  return(summary_df)
}


#####*********************************************************************######
#### Popgen relationships ####

GDAR_pred<-function(dat, 
                    column_interest="Percent.Change",
                    z = 0.05
){
  stopifnot(column_interest %in% colnames(dat))
  Mx <- dat[[column_interest]] / 100  # Convert percentage to fraction
  
  # Ensure Mxzero stays in a valid range [0,1]
  Mxzero <- pmin(pmax(Mx, 0), 1)
  
  # Create new dataset
  myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
  
  # Compute Ax only for valid cases where 1 - Mxzero is non-negative
  valid_indices <- (1 - myd$Mxzero) >= 0
  myd$Ax <- NA  # Initialize with NA
  myd$Ax <- 1 - (1 - myd$Mx)^(1/z)
  # myd$Ax[valid_indices] <- 
  # 1 - (1 - myd$Mx[valid_indices])^(1/z)
  # 
  # Ensure Ax remains in a valid range
  myd$Ax <- pmax(pmin(myd$Ax, 1), 0)  # Clamp Ax between 0 and 1
  return(myd)
}

plot_GDAR<-function(){
  ggplot(myd) +
    geom_point(aes(x = Ax*100, y = -(Mxzero)*100 )) +
    labs(
      # title = "Power-Law Relationship of Area Loss vs Genetic Diversity Loss",
      y = "Genetic diversity loss (π %)",
      x = "Area or population loss (%)") +
    scale_x_continuous(breaks = seq(0,100,by=25),
                       labels = paste("-",seq(0,100,by=25)))+
    theme_minimal()
}


generate_pvalues_withpower<-function(n_sim=1000,power_target=0.33){
  set.seed(123)
  
  # Parameters
  # n_sim <- 100000  # Number of simulated tests
  # power_target <- 0.33  # Desired power
  alpha <- 0.05  # Significance level
  
  # Calculate effect size corresponding to 33% power
  # For simplicity, assuming t-test with df = large, so approx normal
  z_alpha <- qnorm(1 - alpha / 2)
  z_power <- qnorm(power_target + (1 - power_target) * alpha)
  
  effect_size <- z_power - z_alpha
  
  # Simulate test statistics under alternative
  z_values <- rnorm(n_sim, mean = effect_size, sd = 1)
  
  # Convert to two-sided p-values
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  # # Plot histogram
  # hist(p_values, breaks = 50, main = "P-value Distribution under 33% Power",
  #      xlab = "P-value", col = "lightblue", freq = FALSE)
  # 
  # # Add uniform null for comparison
  # curve(dunif(x, 0, 1), col = "red", lwd = 2, add = TRUE)
  # legend("topright", legend = c("33% Power", "Null Uniform"), col = c("lightblue", "red"), lwd = 2)
}

## FST quick effect

# Function to compute fraction of diversity lost
fraction_diversity_lost_FST <- function(k, x, Fst) {
  if (x >= k) {
    stop("Number of demes lost (x) must be less than total demes (k).")
  }
  
  term1 <- (1 - Fst) / (k - x)
  term2 <- ((k - x - 1) * (k - (1 - Fst))) / ((k - x) * (k - 1))
  
  frac_lost <- 1 - (term1 + term2)
  
  return(frac_lost)
}
if(DEBUG){
  fraction_diversity_lost(k = 4, x = 2, Fst = 0.2)
  fraction_diversity_lost(k = 10, x = 5, Fst = 0.2)
  fraction_diversity_lost(k = 100, x = 50, Fst = 0.2)
  
  for(i in c(4,5,10,20,25,30,50,100,1000)){
    fraction_diversity_lost(k = i, x = round(0.5*i), Fst = 0.2) %>% print
  }  
}

# M_fraction_loss_panmictic_drift(){
# # Parameters
# N0 <- 1000  # initial population size
# N1 <- 100   # bottleneck population size
# 
# # Proportions based on neutral SFS
# log_2N <- log(2 * N0)
# 
# P_rare <- log(0.01 * 2 * N0) / log_2N
# P_uncommon <- (log(0.1) - log(0.01)) / log_2N
# P_common <- (log(0.5) - log(0.1)) / log_2N
# 
# # Representative allele frequencies for each class
# p_hat_rare <- 0.005  # 0.5%
# p_hat_uncommon <- 0.05  # 5%
# p_hat_common <- 0.2  # 20%
# 
# # Loss probabilities per class (per generation)
# loss_rare <- (1 - p_hat_rare)^(2 * N1)
# loss_uncommon <- (1 - p_hat_uncommon)^(2 * N1)
# loss_common <- (1 - p_hat_common)^(2 * N1)
# 
# # Weighted sum: fraction of alleles lost
# X_frac_lost <- P_rare * loss_rare + P_uncommon * loss_uncommon + P_common * loss_common
# 
# cat("Fraction of alleles lost (X):", round(X_frac_lost, 4), "\n")
# }

#####*********************************************************************######
#### Inverse modeling functions #### 

##### immediate  ######
#' Title
#'
#' @return
#' @export
#'
#' @examples
M_inv_immediate_population_loss_panmictic<-function(dat,N0=10000){
  
  # Estimate function
  M_estim_subsample<-function(X,N0){ Nx = 1-N0^(-X); return(Nx) }
  
  # Convert to fraction lost (without negative value)
  Mx <- -dat$Percent.Change / 100  # Convert percentage to fraction
  
  # Ensure Mxzero stays in a valid range [0,1]
  Mxzero <- pmin(pmax(Mx, 0), 1)
  
  myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
  
  # Compute Ax only for valid cases where 1 - Mxzero is non-negative
  valid_indices <- (1 - Mxzero) >= 0
  
  myd$Nx <- NA  # Initialize with NA
  
  # Computer estimate of what would be the population bottleneck needed
  myd$Nx[valid_indices] <- M_estim_subsample(Mxzero,N0) 
  
  myd$Nx
  
  RES<-myd %>% 
    mutate(Pop.Bottleneck=Nx*100,
           Pop.start=N0,
           Percent.Change=-(Mxzero)*100
    )
  return(RES)
}
if(DEBUG){  
  # Test
  dat=data.frame(Percent.Change=seq(0,-100,by=-1))
  dat=M_inv_immediate_population_loss_panmictic(dat,10000)
  dat2=M_inv_immediate_population_loss_panmictic(dat,1000)
  dat3=M_inv_immediate_population_loss_panmictic(dat,100)
  dat4=M_inv_immediate_population_loss_panmictic(dat,1e6)
  qplot()+
    geom_line(data=dat,aes(y=Percent.Change,x=Pop.Bottleneck))+
    geom_line(data=dat2,aes(y=Percent.Change,x=Pop.Bottleneck))+
    geom_line(data=dat3,aes(y=Percent.Change,x=Pop.Bottleneck))+
    geom_line(data=dat4,aes(y=Percent.Change,x=Pop.Bottleneck))+
    xlab("Population bottleneck (%)")+
    ylab("Percent M diversity change (%)")+
    xlim(0,100)+ylim(-100,0)
}

#' Title
#'
#' @param dat 
#'
#' @return
#' @export
#'
#' @examples
G_inv_immediate_population_loss_panmictic<-function(dat,minvalue=1e-4){
  # Estimate function
  # Equation is N1= - 1/X, but if we do the negative of percentage change not necessary
  G_estim_subsample<-function(X,minvalue){ N = (1/(X+minvalue)); return(N) } 
  
  # Convert to fraction lost (without negative value)
  Mx <-  -dat$Percent.Change / 100  # Convert percentage to fraction
  
  # Ensure Mxzero stays in a valid range [0,1]
  Mxzero <- pmin(pmax(Mx, 0), 1)
  Mxzero[is.na(Mxzero)]<-0
  
  myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
  
  # Compute Ax only for valid cases where 1 - Mxzero is non-negative
  valid_indices <- (1 - Mxzero) >= 0
  
  myd$Nx <- NA  # Initialize with NA
  
  # Computer estimate of what would be the population bottleneck needed
  myd$Nx[valid_indices] <- G_estim_subsample(Mxzero,minvalue) %>% round
  
  myd$Nx
  
  
  RES<-myd %>% 
    mutate(Pop.Bottleneck=Nx,
           Percent.Change=-(Mxzero)*100
    )
  return(RES)
}
if(DEBUG){  
  # Test
  dat=data.frame(Percent.Change=seq(0,-100,by=-1))
  dat=G_inv_immediate_population_loss_panmictic(dat)
  qplot(y=dat$Percent.Change ,x= dat$Pop.Bottleneck, geom="point")+scale_x_log10()
}

G_testdouble<-function(dat,minvalue=1e-4,myn){
  # Estimate function
  # Equation is N1= - 1/X, but if we do the negative of percentage change not necessary
  G_estim_subsample<-function(X,minvalue,n){ N = -1/(2*(1-((X+1+minvalue)/(1-(1/myn))))) ; return(N) } 
  
  # Convert to fraction lost (without negative value)
  Mx <-  -dat$Percent.Change / 100  # Convert percentage to fraction
  
  # Ensure Mxzero stays in a valid range [0,1]
  Mxzero <- pmin(pmax(Mx, 0), 1)
  Mxzero[is.na(Mxzero)]<-0
  
  myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
  
  # Compute Ax only for valid cases where 1 - Mxzero is non-negative
  valid_indices <- (1 - Mxzero) >= 0
  
  myd$Nx <- NA  # Initialize with NA
  
  # Computer estimate of what would be the population bottleneck needed
  myd$Nx[valid_indices] <- G_estim_subsample(Mxzero,minvalue,myn) %>% round
  
  
  RES<-myd %>% 
    mutate(Pop.Bottleneck=Nx,
           Percent.Change=-(Mxzero)*100
    )
  return(RES)
}
if(DEBUG){  
  # Test
  dat=data.frame(Percent.Change=seq(0,-100,by=-1), n=10)
  dat=G_testdouble(dat, myn=dat$n)
  qplot(y=dat$Percent.Change ,x= dat$Pop.Bottleneck, geom="point")
}



#' Title
#'
#' @param dat 
#' @param generations 
#'
#' @return
#' @export
#'
#' @examples
M_inv_drift_population_loss<- function(dat,N0){
    # check columns
    stopifnot("Percent.Change" %in% colnames(dat))
    stopifnot("Gen.interval" %in% colnames(dat))
  
    # function(dat,N0, Xm_vec, t_vec) {
    ## Function to estimate
    M_estim_drift<-function(N0, Xm_vec, t_vec){
      # Log base
      log_2N <- log(2 * N0)
      
      # Proportions from SFS
      p_hat_rare = 0.005
      p_hat_uncommon = 0.05
      p_hat_common = 0.2
      P_rare <- log(0.01 * 2 * N0) / log_2N
      P_uncommon <- (log(0.1) - log(0.01)) / log_2N
      P_common <- (log(0.5) - log(0.1)) / log_2N
      
      # Ensure vectors are of equal length
      if (length(Xm_vec) != length(t_vec)) {
        stop("Xm_vec and t_vec must be the same length.")
      }
      
      # Vectorized apply function
      N1_results <- mapply(function(Xm, t) {
        loss_function <- function(N1) {
          drift_factor <- 2 * N1 * t
          loss_rare <- (1 - p_hat_rare)^drift_factor
          loss_uncommon <- (1 - p_hat_uncommon)^drift_factor
          loss_common <- (1 - p_hat_common)^drift_factor
          
          Xm_predicted <- P_rare * loss_rare +
            P_uncommon * loss_uncommon +
            P_common * loss_common
          
          return(Xm_predicted - Xm)
        }
        
        # Solve for N1
        result <- tryCatch({
          round(uniroot(loss_function, lower = 1, upper = 100000)$root)
        }, error = function(e) {
          NA  # return NA if root finding fails
        })
        
        return(result)
      }, Xm_vec, t_vec)
      
      return(N1_results)
    }
    ## Extract generations
    generations<-dat$Gen.interval
    
    ## Convert to fraction lost (without negative value)
    Mx <- -dat$Percent.Change / 100  # Convert percentage to fraction
    
    # Ensure Mxzero stays in a valid range [0,1]
    Mxzero <- pmin(pmax(Mx, 0), 1)
    Mxzero[is.na(Mxzero)]<-0
    
    # Start dataset
    myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
    
    # Compute Ax only for valid cases where 1 - Mxzero is non-negative
    valid_indices <- (1 - Mxzero) >= 0
    myd$Nx <- NA  # Initialize with NA
    
    # Computer estimate of what would be the population bottleneck needed
    myd$Nx[valid_indices] <- M_estim_drift(N0,Mxzero, generations)
    
    # Return
    RES<-myd %>% 
      mutate(Pop.Bottleneck=myd$Nx,
             Percent.Change=-(Mxzero)*100
      )
    return(RES)
  }
if(DEBUG){
  # all one generation
  dat=data.frame(Percent.Change=seq(1,-99,by=-1),
                 Gen.interval=rep(1,length(seq(1,-99,by=-1)))                     )
  dat=M_inv_drift_population_loss(dat,1000,dat$Gen.interval)
  qplot(y=dat$Percent.Change ,x= dat$Pop.Bottleneck, geom="point")+
    scale_x_log10()
  # all poisson generations
  dat=data.frame(Percent.Change=seq(1,-99,by=-1),
                 Gen.interval=rpois(length(seq(1,-99,by=-1)),3))
  dat=M_inv_drift_population_loss(dat,1000,dat$Gen.interval)
  qplot(y=dat$Percent.Change ,x= dat$Pop.Bottleneck, geom="point")+
    scale_x_log10()
}


#' Title
#'
#' @param dat 
#' @param generations 
#'
#' @return
#' @export
#'
#' @examples
G_inv_drift_population_loss<-function(dat,generations=1){
  # transform funciton
  G_estim_drift <- function(X, t) {
    # X is the % loss of diversity
    retained <- 1 - X / 100
    N <- 1 / (2 * (1 - retained^(1 / t)))
    
    return(N)
  }  
  ## Convert to fraction lost (without negative value)
  Mx <- -dat$Percent.Change / 100  # Convert percentage to fraction
  
  # Ensure Mxzero stays in a valid range [0,1]
  Mxzero <- pmin(pmax(Mx, 0), 1)
  Mxzero[is.na(Mxzero)]<-0
  
  myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
  
  # Compute Ax only for valid cases where 1 - Mxzero is non-negative
  valid_indices <- (1 - Mxzero) >= 0
  
  myd$Nx <- NA  # Initialize with NA
  
  # Computer estimate of what would be the population bottleneck needed
  myd$Nx[valid_indices] <- G_estim_drift(Mxzero, generations)
  
  
  RES<-myd %>% 
    mutate(Pop.Bottleneck=myd$Nx,
           Percent.Change=-(Mxzero)*100
    )
  return(RES)
}
if(DEBUG){  
  # Test
  dat=data.frame(Percent.Change=seq(-1,-99,by=-1))
  dat=G_inv_drift_population_loss(dat)
  qplot(y=dat$Percent.Change ,x= dat$Pop.Bottleneck, geom="line")
  
  dat=data.frame(Percent.Change=seq(-1,-99,by=-1))
  dat=G_inv_drift_population_loss(dat, rep(10,nrow(dat)))
  qplot(y=dat$Percent.Change ,x= dat$Pop.Bottleneck, geom="line")
}
# Test
#   dat=data.frame(Percent.Change=c(10,50,90), Gen.interval=c(2,3,2))
#   dat=data.frame(Percent.Change=c(1,2,5), Gen.interval=c(2,3,2))
#   tmp=G_inv_drift_population_loss(dat,dat$Gen.interval)
#   plot(tmp$Percent.Change ~ tmp$Nx)
# 
#   G_estim_drift(c(10,50,90), c(1,1,1))
# 1/(2*Mxzero^(1/generations))
# 
# t=2
# X=0.05
# 
# 
# t=3
# X=0.05
# 1 / ( 2* (1- (1- (X/100)^(1/t) )) )


# π diversity 

#' Infer the immediate area loss based on area relationships
#'
#' @param dat Contains a column Percent.Change used to do inverse modeling
#' @param z The parameter of the power law G=A^z (default 0.05)
#'
#' @return
#' @export
#'
#' @examples
G_inv_immediate_area_loss_GDAR<-function(dat, z=0.05){
  
  # Convert to fraction lost (without negative value)
  Mx <- -dat$Percent.Change / 100  # Convert percentage to fraction
  
  # Ensure Mxzero stays in a valid range [0,1]
  Mxzero <- pmin(pmax(Mx, 0), 1)
  Mxzero[is.na(Mxzero)]<-0
  
  #Initialize results  
  myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
  
  # Compute Ax only for valid cases where 1 - Mxzero is non-negative
  valid_indices <- (1 - Mxzero) >= 0
  myd$Ax <- NA  # Initialize with NA
  myd$Ax[valid_indices] <- 1 - (1 - Mxzero[valid_indices])^(1/z)
  
  # Ensure Ax remains in a valid range
  myd$Ax <- pmax(pmin(myd$Ax, 1), 0)  # Clamp Ax between 0 and 1
  
  RES<-myd %>% 
    mutate(Area.Loss=Ax*100,
           Percent.Change=-(Mxzero)*100
    )
  return(RES)
}
if(DEBUG){  
  # Test
  dat=data.frame(Percent.Change=seq(1,-99,by=-1))
  dat=G_inv_immediate_area_loss_GDAR(dat)
  qplot(y=dat$Percent.Change ,x= dat$Area.Loss, geom="line") %>% print
}

#' Infer the immediate area loss based on area relationships
#'
#' @param dat Contains a column Percent.Change used to do inverse modeling
#' @param z The parameter of the power law M=A^z (default 0.3)
#'
#' @return
#' @export
#'
#' @examples
M_inv_immediate_area_loss_MAR<-function(dat, z=0.2){
  
  Mx <- -dat$Percent.Change / 100  # Convert percentage to fraction
  
  # Ensure Mxzero stays in a valid range [0,1]
  Mxzero <- pmin(pmax(Mx, 0), 1)
  Mxzero[is.na(Mxzero)]<-0
  
  #Initialize results  
  myd <- data.frame(Mx = Mx, Mxzero = Mxzero)
  
  # Compute Ax only for valid cases where 1 - Mxzero is non-negative
  valid_indices <- (1 - Mxzero) >= 0
  myd$Ax <- NA  # Initialize with NA
  myd$Ax[valid_indices] <- 1 - (1 - Mxzero[valid_indices])^(1/z)
  
  # Ensure Ax remains in a valid range
  myd$Ax <- pmax(pmin(myd$Ax, 1), 0)  # Clamp Ax between 0 and 1
  
  RES<-myd %>% 
    mutate(Area.Loss=Ax*100,
           Percent.Change=-(Mxzero)*100
    )
  return(RES)
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,-99,by=-1))
  dat=M_inv_immediate_area_loss_MAR(dat)
  qplot(y=dat$Percent.Change ,x= dat$Area.Loss, geom="line") %>% print
}


#' Title
#'
#' @param dat 
#' @param wfmoments 
#' @param weighting_function 
#'
#' @return
#' @export
#'
#' @examples
# Possible weight functions
weighting_function_Gaussian <- function(d) exp(-d^2) # Gaussian Kernel
weighting_function_Gaussian_Quantile <- function(d, myquantile=0.5){
  probs<-exp(-d^2) # Gaussian Kernel
  probs[probs<quantile(probs,probs = myquantile)]<-0
  return(probs)
}
weighting_function_Exp <- function(d) d==min(d) # Exponential Kernel
weighting_function_Minimum <- function(d) abs(d)==min(abs(d)) # 
weighting_function_Linear <- function(d) {(max(abs(d))-abs(d))/(max(abs(d))-min(abs(d))) } # Minimum Kernel
# weighting_function_Minimum <- function(d) abs(d)) # Linear Kernel
weighting_function_Quantile <- function(d, myquantile=0.1){
  as.numeric(d<quantile(d,myquantile)) # Minimum Kernel
}

# Main function
G_inv_drift_area_loss_WFmoments<-function(dat, wfmoments, weighting_function=function(d) exp(-d^2)){
  
  # Initialize
  Mx<-dat$Percent.Change # Empirical gendiv change from study
  MxSD<-dat$SD.Percent.Change # Noise in gendiv change from study
  Mxsamples<-c() # Initialized vector
  Mxindex<-c() # Initialized vector
  
  Mxsim<-wfmoments$Percent.Change # Simulated gendiv change under area contraction
  Axsim<-wfmoments$Area.Loss # Simulated area contraction
  Axhat<-c()
  
  # Scale values that make sense -100% to +100%, the others are outliers
  Mx <- pmin(pmax(Mx, -100), +0)
  # Mx <- pmin(pmax(Mx, min(Mxsim)), max(Mxsim)) # POSSIBLE BUG! bound it to the simulated data
  
  # Iterate over all empirical values and get expectation based on simulation
  i=1
  for(i in 1:length(Mx)){
    Mxsample= Mx[i]
    D<- Mxsample-Mxsim
    probs<-weighting_function(D)
    Ax<- weighted.mean(Axsim,probs,na.rm = T)
    # print(Ax)
    Axhat<-c(Axhat,Ax)
  }
  # Save in dataset
  res<-
    data.frame(Percent.Change=Mx, Area.Loss=Axhat)
  return(res)
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,-99,by=-1))
  dat=G_inv_drift_area_loss_WFmoments(dat, wfedge)
  qplot(y=dat$Percent.Change ,x= dat$Area.Loss, geom="point") %>% print
}

#' Title
#'
#' @param dat 
#' @param wfmoments 
#' @param weighting_function 
#'
#' @return
#' @export
#'
#' @examples
G_inv_drift_area_loss_WFmoments_sample<-function(dat, wfmoments, weighting_function=function(d) exp(-d^2), debug=F){
  
  # Initialize
  Mx<-dat$Percent.Change # Empirical gendiv change from study
  MxSD<-dat$SD.Percent.Change # Noise in gendiv change from study
  Mxsamples<-c() # Initialized vector
  Mxindex<-c() # Initialized vector
  
  Mxsim<-wfmoments$Percent.Change # Simulated gendiv change under area contraction
  Axsim<-wfmoments$Area.Loss # Simulated area contraction
  Axhat<-c()
  
  # Scale values that make sense -100% to +100%, the others are outliers
  # Mx <- pmin(pmax(Mx, -100), +0) #if prefer typical reange
  Mx <- pmin(pmax(Mx, min(Mxsim)), max(Mxsim)) # POSSIBLE BUG! bound it to the simulated data
  Mx[is.na(Mx)]<-0
  
  # Iterate over all empirical values and get expectation based on simulation
  i=1
  for(i in 1:length(Mx)){
    if(debug) print(i)
    Mxsample= Mx[i]
    if(debug)print(Mxsample)
    D<- Mxsample-Mxsim
    probs<-weighting_function(D)
    if(debug)hist(probs)
    Ax<- sample(Axsim,1,prob = probs,replace=T) 
    # Ax<- sample(Axsim,10,prob = probs,replace=T) %>% mean  Ax<- sample(Axsim,10,prob = probs,replace=T) %>% mean
    Axhat<-c(Axhat,Ax)
  }
  # Assign zero to the positive trend
  # Axhat[dat$Percent.Change>0]<-0
  # Save in dataset
  res<-
    data.frame(Percent.Change=Mx, Area.Loss=Axhat)
  return(res)
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,-99,by=-1))
  dat=G_inv_drift_area_loss_WFmoments_sample(dat, wfedge,weighting_function_Gaussian)
  qplot(y=dat$Percent.Change ,x= dat$Area.Loss, geom="point") %>% print
}

G_inv_drift_area_loss_WFmoments_sample_future<-function(dat, wfmoments, weighting_function=function(d) exp(-d^2), debug=F){
  
  # Initialize
  Mx<-dat$Percent.Change # Empirical gendiv change from study
  MxSD<-dat$SD.Percent.Change # Noise in gendiv change from study
  Mxsamples<-c() # Initialized vector
  Mxindex<-c() # Initialized vector
  
  Mxsim<-wfmoments$Percent.Change.Fut # Simulated gendiv change under area contraction
  Axsim<-wfmoments$Area.Loss # Simulated area contraction
  Axhat<-c()
  
  # Scale values that make sense -100% to +100%, the others are outliers
  # Mx <- pmin(pmax(Mx, -100), +0) #if prefer typical reange
  Mx <- pmin(pmax(Mx, min(Mxsim)), max(Mxsim)) # POSSIBLE BUG! bound it to the simulated data
  Mx[is.na(Mx)]<-0
  
  # Iterate over all empirical values and get expectation based on simulation
  i=1
  for(i in 1:length(Mx)){
    if(debug) print(i)
    Mxsample= Mx[i]
    if(debug)print(Mxsample)
    D<- Mxsample-Mxsim
    probs<-weighting_function(D)
    if(debug)hist(probs)
    Ax<- sample(Axsim,1,prob = probs,replace=T) 
    # Ax<- sample(Axsim,10,prob = probs,replace=T) %>% mean  Ax<- sample(Axsim,10,prob = probs,replace=T) %>% mean
    Axhat<-c(Axhat,Ax)
  }
  # Assign zero to the positive trend
  # Axhat[dat$Percent.Change>0]<-0
  # Save in dataset
  res<-
    data.frame(Percent.Change=Mx, Area.Loss=Axhat)
  return(res)
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,-99,by=-1))
  dat=G_inv_drift_area_loss_WFmoments_sample(dat, wfedge,weighting_function_Gaussian)
  qplot(y=dat$Percent.Change ,x= dat$Area.Loss, geom="point") %>% print
}

exponential_interpolation_t <- function(D0, D1, N, g) {
  p <- g / (4 * N)
  D0 * (D1 / D0)^p
}
exponential_interpolation <- function(D0, D1, p) {
  D0 * (D1 / D0)^p
}
#' Title
#'
#' @param dat 
#' @param wfmoments 
#' @param weighting_function 
#'
#' @return
#' @export
#'
#' @examples
G_inv_drift_area_loss_WFmoments_drift_sample<-function(dat, Nt, wfmoments, weighting_function=function(d) exp(-d^2), debug=F){
  stopifnot("Percent.Change" %in% colnames(dat))
  stopifnot("Gen.interval" %in% colnames(dat))
  
  # Initialize
  t=generations=dat$Gen.interval
  Mx<-dat$Percent.Change # Empirical gendiv change from study
 
  Mxsamples<-c() # Initialized vector
  
  Mxsim<-wfmoments$Percent.Change # Simulated gendiv change under area contraction
  Mxsimfut<-wfmoments$Percent.Change.Fut # Simulated gendiv change under area contraction future equilibrium
  Axsim<-wfmoments$Area.Loss # Simulated area contraction
  Axhat<-c()
  
  # Scale values that make sense -100% to +100%, the others are outliers
  Mx <- pmin(pmax(Mx, min(Mxsim,Mxsimfut)), max(Mxsim,Mxsimfut)) # POSSIBLE BUG! this bound may necessarily be the future one?
  Mx[is.na(Mx)]<-0
  
  # Iterate over all empirical values and get expectation based on simulation
  i=1
  for(i in 1:length(Mx)){
    # Get the observed value and generation interval
    Mxsample= Mx[i]
    tsample= t[i]
    # Transform simulated data adding time
    timeMxsim<-exponential_interpolation_t(Mxsim,Mxsimfut, Nt, tsample)
      timeMxsim[is.na(timeMxsim)]<-0
    # print(i)
    # print(Mxsample)
    #Distance
    D<- abs(Mxsample-timeMxsim)
    # scale distance to avoid error in probabilities
    D=D-min(D)
    # Compute probabilty based on weighting or kernel function
    probs<-weighting_function(D)
    # Sample area of matched simulated values and og
    Ax<- sample(Axsim,1,prob = probs,replace=T) 
    Axhat<-c(Axhat,Ax)
  }
  # Assign zero to the positive trend
  # Axhat[dat$Percent.Change>0]<-0 # we do not want this for generality
  # Save in dataset
  res<-
    data.frame(Percent.Change=Mx, Area.Loss=Axhat)
  return(res)
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,-99,by=-1), generations=rpois(101,lambda = 3))
  dat=G_inv_drift_area_loss_WFmoments_drift_sample(dat,generations = dat$generations,
                                                  Nt = 100, wfedge,weighting_function_Gaussian)
  qplot(y=dat$Percent.Change ,x= dat$Area.Loss, geom="point") %>% print
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,99,by=1), Gen.interval=rpois(99,lambda = 3))
  dat=G_inv_drift_area_loss_WFmoments_drift_sample(dat,
                                                  Nt = 10000, 
                                                  wffrag,
                                                  weighting_function_Gaussian)
  qplot(y=dat$Percent.Change ,x= dat$Area.Loss, geom="point") %>% print
}



#' Title
#'
#' @param dat 
#' @param wfmoments 
#' @param weighting_function 
#'
#' @return
#' @export
#'
#' @examples

G_inv_drift_area_loss_WFmoments_drift<-function(dat, wfmoments, weighting_function=function(d) exp(-d^2)){
  
  # Initialize
  Mx<-dat$Percent.Change # Empirical gendiv change from study
  MxSD<-dat$SD.Percent.Change # Noise in gendiv change from study
  Mxsamples<-c() # Initialized vector
  Mxindex<-c() # Initialized vector
  
  Mxsim<-wfmoments$Percent.Change # Simulated gendiv change under area contraction
  Axsim<-wfmoments$Area.Loss # Simulated area contraction
  Axhat<-c()
  
  # Scale values that make sense -100% to +100%, the others are outliers
  Mx <- pmin(pmax(Mx, -100), +100)
  Mx <- pmin(pmax(Mx, min(Mxsim)), max(Mxsim)) # POSSIBLE BUG! bound it to the simulated data
  Mx[is.na(Mx)]<-0
  
  # Iterate over all empirical values and get expectation based on simulation
  i=1
  for(i in 1:length(Mx)){
    Mxsample= Mx[i]
    D<- Mxsample-Mxsim
    probs<-weighting_function(D)
    Ax<- weighted.mean(Axsim,probs,na.rm = T)
    # print(Ax)
    Axhat<-c(Axhat,Ax)
  }
  # Save in dataset
  res<-
    data.frame(Percent.Change=Mx, Area.Loss=Axhat)
  return(res)
}


#####*********************************************************************######
#### Re-sample noise ####

#' Resample dataset with noise
#'
#' @param dat 
#' @param reps 
#'
#' @return
#' @export
#'
#' @examples
resample<-function(dat,reps=1){
  sim<-dat
  sim$Percent.Change<- 
    apply(dat,1,FUN = function(i){ rnorm(n=1 , mean=as.numeric(i["Percent.Change"]), 
                                         sd=as.numeric(i["SD.Percent.Change"])) 
    })
  # Iterate
  if(reps>1){
    for( i in 2:reps){
      tmp<-dat
      # Sample following norm
      tmp$Percent.Change<- 
        apply(dat,1,FUN = function(i){ rnorm(n=1 , mean=as.numeric(i["Percent.Change"]), 
                                             sd=as.numeric(i["SD.Percent.Change"])) 
        })
      sim<-rbind(sim, tmp)
    } # for
  } # if
  
  # Return
  return(sim)
}



wrapper_model<-function(datori, replicates, relevantvariable ,
                        inferencefunction,...){
  
  # Reample dataset based on noise
  datsim<-resample(datori,replicates)
  datsim$Percent.Change[is.na(datsim$Percent.Change)]<-0
  
  # Inverse modeling  
  
  resori<-inferencefunction(datori,...)
  res=inferencefunction(datsim,...)
  
  # Put together dataset with extra info on trends
  datori[[relevantvariable]]<-resori[[relevantvariable]]
  datori$pvalue_color<-datori$pvalue_color
  datsim[[relevantvariable]]<-res[[relevantvariable]]
  datsim$pvalue_color<-rep(datori$pvalue_color,replicates) # nolte this will keep filling
  
  # Add the original to the simulated
  datsim[[paste0(relevantvariable,".original")]]<-rep(resori[[relevantvariable]],replicates)
  datsim$Percent.Change.original<-rep(resori$Percent.Change,replicates)
  
  return(datsim)
}

#####*********************************************************************######
#### Diversity increase ####


# Relative diversity at time t after expansion
pi_growth <- function(N0, N1, t) {
  pi_ratio <- 1 - (1 - N0 / N1) * (1 - 1 / (2 * N1))^t
  return(pi_ratio)  # π_t / π_infinity
}

# Example usage
pi_growth(N0 = 100, N1 = 1000, t = 500)  # Fraction of equilibrium π reached



#' Title
#'
#' @param dat 
#' @param generations 
#'
#' @return
#' @export
#'
#' @examples
# Estimate population size change (Xn) from genetic diversity change (Xt), given N1 and t
G_inv_growth <- function(dat, N1, t) {
  # Mx: proportional change in pi (e.g. 0.2 for +20%)
  # N1: current population size
  # t: generations since population size change
  
  Mx<-dat$Percent.Change/100 # Empirical gendiv change from study
  Mx <- pmin(pmax(Mx, 0), +1)
  
  # Avoid division by 0
  decay_factor <- 1 - (1 - 1 / (2 * N1))^t
  
  # If decay_factor is too small (e.g., t = 0), return NA
  Xn <- ifelse(decay_factor > 0, Mx / decay_factor, NA)
  
  res<-data.frame(Percent.Change=(Mx)*100,
                  Pop.Growth=Xn*100)
  return(res)
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,99,by=1), 
                 generations=rpois(99,lambda = 3))
  dat=G_inv_growth(dat,
                   N1=1000,
                   t = dat$generations
                   )
  qplot(y=dat$Percent.Change ,
        x= dat$Pop.Growth, geom="point") %>% print
}


G_inv_Fst<-function(dat){
  # Mx: proportional change in pi (e.g. 0.2 for +20%)
  
  Mx<-dat$Percent.Change/100 # Empirical gendiv change from study
  Mx <- pmin(pmax(Mx, 0), +1)
  
  Fst<-Mx/(1+Mx)
  
  res<-data.frame(Percent.Change=(Mx)*100,
                  Fst=Fst)
  return(res) 
}
if(DEBUG){
  # Test
  dat=data.frame(Percent.Change=seq(1,99,by=1))
  dat=G_inv_Fst(dat)
  qplot(y=dat$Percent.Change ,x= dat$Fst, geom="point") %>% print
}

