library(dplyr)  # Data manipulation
library(purrr)  # Functional programming tools
library(tibble) # Data frames
library(lme4)   # Generalized linear mixed models
library(numDeriv) # Numerical derivatives

# Helper function: obtain gradients and Hessians of the GLMM model
# Parameters:
# - formula: Formula specifying the model structure.
# - data: Data to fit the model.
# - beta_prior: Prior values for fixed effects.
# - theta: Prior values for variance component parameters.
# - ...: Additional arguments for the GLMM.
get_glmer_derivative <- function(formula, data, beta_prior, theta, ...){
  # Get deviance function for GLMM
  devfun <- glmer(formula, data, devFunOnly = T, ...)
  
  # Calculate gradient and Hessian using numerical differentiation
  gradient <- grad(devfun, c(theta, beta_prior)); Hessian <- hessian(devfun, c(theta, beta_prior))
  
  return(list(gradient = gradient, Hessian = Hessian))
}

# Helper function: Rearranges fixed effects for each site
# Parameters:
# - deriv: Derivative vector or matrix (gradient or Hessian).
# - FE_site_loc: Indices of site-specific fixed effects. NULL if all fixed effects are shared.
# - n_theta: Number of variance component parameters.
get_FE_site_rearranged <- function(deriv, FE_site_loc, n_theta){
  if(is.null(FE_site_loc)){
    return(deriv)
  } else {
    # Adjust indices of site-varying fixed effects by n_theta slots 
    FE_site_loc_full <- FE_site_loc + n_theta
    
    # Rearrange gradient
    if(is.vector(deriv)){
      return(deriv[-FE_site_loc_full] %>% c(deriv[FE_site_loc_full]))
    } else {
      # Rearrange Hessian
      deriv_tranCol <- deriv[,-FE_site_loc_full] %>% cbind(deriv[,FE_site_loc_full]) 
      return(deriv_tranCol[-FE_site_loc_full,] %>% rbind(deriv_tranCol[FE_site_loc_full,])) 
    }  
  }
}

# Helper function: One iteration of Fed-GLMM estimation
# Parameters:
# - formula: Formula specifying the model structure.
# - data: Data to fit the model.
# - site_var: Column name for site variable.
# - FE_site_loc: Indices of site-specific fixed effects. NULL if all fixed effects are shared.
# - estimate_prior: Initial estimates for variance components and fixed effects.
# - ...: Additional arguments for the GLMM.
get_one_iteration <- function(formula, data, site_var, FE_site_loc, estimate_prior, common_theta = F, ...){
  # Identify the number of sites (K), variance component parameters (n_theta), and site-varying fixed effects (n_alpha)
  site <- data[[site_var]]
  K <- length(unique(site))
  n_alpha <- length(FE_site_loc)
  n_theta <- glFormula(formula, data, ...)$reTrms$theta %>% length
  
  # Split prior estimates for variance component parameters and site-sharing fixed effects
  if(common_theta == T){
    K_theta <- 1
    theta_list <- rep(list(estimate_prior[1: n_theta]), K)
  } else {
    K_theta <- K
    theta_list <- estimate_prior[(1): (K*n_theta)] %>% split(ceiling(seq_along(.)/n_theta))  
  }
  
  n_beta <- length(estimate_prior) - K_theta*n_theta - K*n_alpha 
  beta_shared <- estimate_prior[(K_theta*n_theta + 1): (K_theta*n_theta + n_beta)]
  
  # Split data by site
  site_data <- map(unique(site), ~data[site == .,])
  
  # Create a list of full fixed effect estimates for each site
  beta_list <- list(beta_shared) %>% rep(K)
  if(!is.null(FE_site_loc)){
    beta_site_list <- estimate_prior[(K_theta*n_theta + n_beta + 1): (K_theta*n_theta + n_beta + K*n_alpha)] %>% split(ceiling(seq_along(.)/n_alpha))
    for(i in 1:length(FE_site_loc)){
      beta_list <- map2(beta_list, beta_site_list, ~append(.x, .y[i], after = FE_site_loc[i]-1))
    }  
  }
  
  # Compute gradients and Hessians for each site
  derivative_list_full <- derivative_list <- glmer_derivative <- pmap(list(site_data, beta_list, theta_list), possibly(function(A, B, C){get_glmer_derivative(formula, A, B, C, ...)}, otherwise = NULL))
  
  # Handle missing gradient/Hessian by averaging existing ones if needed
  if(length(compact(glmer_derivative)) < length(unique(site))){
    dev_gradient_avg <- map(derivative_list %>% compact, ~.$gradient) %>% Reduce("+", .) %>% `/`(length(derivative_list %>% compact))  
    dev_Hessian_avg <- map(derivative_list %>% compact, ~.$Hessian) %>% Reduce("+", .) %>% `/`(length(derivative_list %>% compact)) 
    derivative_list_full <- derivative_list %>% map(function(x){if(is.null(x)){list(gradient = dev_gradient_avg, Hessian = dev_Hessian_avg)} else {x}})
  }
  
  # Rearrange the gradients and Hessians for fixed effects
  gradient_list_rearranged <- derivative_list_full %>% map(~.$gradient/(-2)) %>% map(~get_FE_site_rearranged(.,FE_site_loc, n_theta))
  Hessian_list_rearranged <- derivative_list_full %>% map(~.$Hessian/(-2)) %>% map(~get_FE_site_rearranged(.,FE_site_loc, n_theta))
  
  # Extract and combine the gradients for variance component parameters, site-sharing fixed effects, and site-varying fixed effects
  if(common_theta == T){
    gradient_theta <- gradient_list_rearranged %>% map(~.[1:n_theta]) %>% Reduce("+", .)
  } else {
    gradient_theta <- gradient_list_rearranged %>% map(~.[1:n_theta]) %>% unlist  
  }
  gradient_beta <- gradient_list_rearranged %>% map(~.[(n_theta + 1):(n_theta + n_beta)]) %>% Reduce("+", .)
  
  if(!is.null(FE_site_loc)){
    gradient_alpha <- gradient_list_rearranged %>% map(~.[(n_theta + n_beta + 1):(n_theta + n_beta + n_alpha)]) %>% unlist  
  } else {
    gradient_alpha <- NULL
  }
  gradient_t <- c(gradient_theta, gradient_beta, gradient_alpha)
  
  # Initialize Hessian matrix
  Hessian_t <- matrix(0, K_theta*n_theta + n_beta + K*n_alpha, K_theta*n_theta + n_beta + K*n_alpha) %>% 
    as(Class = "sparseMatrix")
  
  # Fill the Hessian matrix by combining Hessians from all sites
  for(k in 1:K){
    
    Hessian_site <- Hessian_list_rearranged[[k]]
    
    if(common_theta == T){
      Hessian_t[(1):(n_beta+n_theta), (1):(n_beta+n_theta)] <- Hessian_t[(1):(n_beta+n_theta), (1):(n_beta+n_theta)] + Hessian_site[(1):(n_beta+n_theta), (1):(n_beta+n_theta)]
      if(!is.null(FE_site_loc)){
        Hessian_t[(1):(n_beta+n_theta), (n_beta+n_theta+(k-1)*n_alpha+1):(n_beta+n_theta+k*n_alpha)] <- Hessian_site[(1):(n_beta+n_theta), (n_beta+n_theta+1):(n_beta+n_theta+n_alpha)]
        Hessian_t[(n_beta+n_theta+(k-1)*n_alpha+1):(n_beta+n_theta+k*n_alpha), (1):(n_beta+n_theta)] <- Hessian_site[(n_beta+n_theta+1):(n_beta+n_theta+n_alpha), (1):(n_beta+n_theta)]
        Hessian_t[(n_beta+n_theta+(k-1)*n_alpha+1):(n_beta+n_theta+k*n_alpha), (n_beta+n_theta+(k-1)*n_alpha+1):(n_beta+n_theta+k*n_alpha)] <- Hessian_site[(n_beta+n_theta+1):(n_beta+n_theta+n_alpha), (n_beta+n_theta+1):(n_beta+n_theta+n_alpha)]
      }
      
    } else {
      
      Hessian_t[((k-1)*n_theta+1):(k*n_theta), ((k-1)*n_theta+1):(k*n_theta)] <- Hessian_site[1:n_theta, 1:n_theta]
      Hessian_t[((k-1)*n_theta+1):(k*n_theta), (K*n_theta+1):(K*n_theta+n_beta)] <- Hessian_site[1:n_theta, (n_theta+1):(n_theta+n_beta)]
      
      if(!is.null(FE_site_loc)){
        Hessian_t[((k-1)*n_theta+1):(k*n_theta), (K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha)] <- Hessian_site[1:n_theta, (n_theta+ n_beta+1):(n_theta+ n_beta + n_alpha)]  
      }
      
      Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), ((k-1)*n_theta+1):(k*n_theta)] <- Hessian_site[(n_theta+1):(n_theta+ n_beta), 1:n_theta]
      Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), (K*n_theta+1):(K*n_theta+n_beta)] <- Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), (K*n_theta+1):(K*n_theta+n_beta)] + Hessian_site[(n_theta+1):(n_theta+ n_beta), (n_theta+1):(n_theta+ n_beta)]
      
      if(!is.null(FE_site_loc)){
        Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), (K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha)] <- Hessian_site[(n_theta+1):(n_theta+ n_beta), (n_theta+ n_beta+1):(n_theta+ n_beta + n_alpha)]
        Hessian_t[(K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha), ((k-1)*n_theta+1):(k*n_theta)] <- Hessian_site[(n_theta+ n_beta+1):(n_theta+ n_beta+ n_alpha), 1:n_theta]
        Hessian_t[(K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha), (K*n_theta+1):(K*n_theta+n_beta)] <- Hessian_site[(n_theta+ n_beta+1):(n_theta+ n_beta+ n_alpha), (n_theta+1):(n_theta+ n_beta)]
        Hessian_t[(K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha), (K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha)] <- Hessian_site[(n_theta+ n_beta+1):(n_theta+ n_beta+ n_alpha), (n_theta+ n_beta+1):(n_theta+ n_beta + n_alpha)]
      }  
    }
  }
  
  # Update the estimates using the Newton-Raphson method
  estimate_bar <- as.vector(estimate_prior - solve(Hessian_t) %*% (gradient_t))  # Convert to vector for naming
  
  # Generate parameter names dynamically
  if(common_theta == T){
    theta_names <- paste0("theta", 1:n_theta)
  } else {
    theta_names <- paste0("theta", rep(1:n_theta, times = K), "_site", rep(unique(site), each = n_theta))  
  }
  beta_names <- paste0("beta", 1:n_beta)
  if(n_alpha == 0){alpha_names <- NULL} else {alpha_names <- paste0("alpha", rep(1:n_alpha, times = K), "_site", rep(unique(site), each = n_alpha))}
  param_names <- c(theta_names, beta_names, alpha_names)
  names(estimate_bar) <- param_names
  
  return(list(estimate_bar = estimate_bar, gradient_t = gradient_t, Hessian_t = Hessian_t))
}

# Fed-GLMM function: Iteratively compute GLMM estimates across sites
# Parameters:
# - formula: Formula specifying the model structure.
# - data: Data to fit the model.
# - site_var: Column name for site variable.
# - FE_site_loc: Indices of site-specific fixed effects. NULL if all fixed effects are shared.
# - initial_values: Initial estimates for variance components and fixed effects.
# - n_iter: Number of iterations to run.
# - ...: Additional arguments for the GLMM.
get_FedGLMM_estimates <- function(formula, data, site_var, FE_site_loc, initial_values, n_iter = 5, common_theta = F, ...){
  # Initialize variables to store results
  estimate_history <- list() 
  diff_estimate_bar_history <- list()
  gradient_history <- list()
  Hessian_history <- list()
  
  estimate_history[[1]] <- final_result <- estimate_prior <- initial_values 
  iter_count <- 0
  
  # Iterate to obtain estimates
  for (x in 1:n_iter){
    converge_mark <- F
    converge_iter <- NaN
    diverge_mark <- F
    
    writeLines(paste0("#### iter ", x, " ####\n"))
    
    # Run one iteration to update estimates
    iteration <- get_one_iteration(formula, data, site_var, FE_site_loc, estimate_prior, common_theta, ...)
    
    # Calculate difference from the previous estimate
    diff_estimate_bar <- diff_estimate_bar_history[[x]] <- sqrt(sum((estimate_prior - iteration$estimate_bar)^2))
    writeLines(paste0("Beta Difference: ", diff_estimate_bar))
    
    # Store results from this iteration
    final_result <- iteration$estimate_bar
    estimate_history[[x+1]] <- iteration$estimate_bar 
    gradient_history[[x]] <- iteration$gradient_t
    Hessian_history[[x]] <- iteration$Hessian_t
    iter_count <- x
    
    # Check for divergence
    if (diff_estimate_bar > 200) {
      writeLines("The values have diverged.\n")
      diverge_mark <- T
      break
    }
    
    # Check for convergence
    if(diff_estimate_bar < 1e-6){
      writeLines("The values have converged.\n")
      converge_mark <- T
      converge_iter <- x
      break
    } else {
      estimate_prior <- iteration$estimate_bar
    }
  }
  
  # Assign parameter names to the final result
  final_result <- setNames(as.vector(final_result), names(iteration$estimate_bar))
  
  return(list(final_result = final_result, estimate_history = estimate_history, gradient_history = gradient_history, Hessian_history = Hessian_history, var_cov = solve(-Hessian_history[[length(Hessian_history)]]), diff_estimate_bar_history = diff_estimate_bar_history, iter_count = iter_count, converge_mark = converge_mark, converge_iter = converge_iter, diverge_mark = diverge_mark))
}



# Data simulation: Generates synthetic data for testing Fed-GLMM

# Helper function: generate covariance matrix from variance component parameters
get_covmatrix_block <- function(theta){
  # Generate block covariance matrix based on variance component parameters
  nc <- (sqrt(length(theta)*8 + 1) - 1)/2 
  rowIndices <- rep(1:nc, 1:nc)
  colIndices <- sequence(1:nc)
  templatet <- sparseMatrix(rowIndices, colIndices, x = 1 * (rowIndices == colIndices)) %>% t 
  LindTemplate <- rowIndices + nc * (colIndices - 1) - choose(colIndices, 2) 
  
  templatet@x <- theta[LindTemplate]
  D_block <- t(templatet) %*% templatet 
  return(D_block)
}

# Helper sigmoid function
get_sigmoid <- function(x){
  1/(1+exp(-x))
}

# Parameters:
# - Xall: Matrix of fixed effect covariates for all sites.
# - patient: Vector of patient identifiers.
# - site: Vector of site identifiers.
# - RE_col_ind: Indices for random effects covariates.
# - beta_true_list: List of true fixed effect coefficients for each site.
# - theta_list: List of variance component parameters for each site.
get_simulated_df <- function(Xall, patient, site, RE_col_ind, beta_true_list, theta_list){
  # Set up random effect and fixed effect design matrices for each site
  J_k <- map(unique(site), ~t(as(factor(patient)[site == .], Class = "sparseMatrix"))) 
  X_k <- map(unique(site), ~Xall[site == .,])
  X_RE_k <- map(unique(site), ~Xall[site == .,RE_col_ind]) 
  Z_k <- map2(J_k, X_RE_k, ~t(KhatriRao(t(.x), t(.y))))
  
  # Compute the covariance matrix for variance component parameters for each site
  D_k <- theta_list %>% map(get_covmatrix_block) %>% map2(unique(site), ~.bdiag(rep(list(.x), length(unique(patient[site == .y])))))
  
  # Sample random effects for each site
  b_k <- map(D_k, ~t(chol(.)) %*% rnorm(nrow(.)))
  
  # Compute linear predictor and outcome
  ita_k <- pmap(list(X = X_k, Z = Z_k, b = b_k, beta_true = beta_true_list), ~with(list(...), X %*% beta_true + Z %*% b))
  ita <- do.call(rbind, ita_k) %>% get_sigmoid %>% as.numeric
  y <- rbinom(length(site), 1, ita)
  
  # Return generated dataset as a data frame
  df <- cbind(y, Xall, patient, site, ita) %>% as.matrix %>% as.data.frame %>% rename(Y = 1, X0 = 2)
  return(df)
}

# Function to simulate test data (fixed effects and variance components across sites)
# Meta-analysis to obtain initial parameter estimates
# Perform Fed-GLMM to estimate model parameters iteratively across sites

# Simulate test data
K <- 8 # number of sites
n <- 4 # number of visits per patient
m <- seq(from = 100, to = 50*(K+1), by = 50) # number of patients per site
N <- sum(m)*n # number of total observations
M <- sum(m) # number of total patients
site <- rep(1:K, m*n) # site index for each observation
patient <- rep(1:M, each = n) # patient index for each observation

# Generate site-specific fixed effects
set.seed(100); beta_site <- map(rnorm(K, 0.5, 0.5) %>% round(1), ~c(1.5, 1, ., 0.5, 0.5)) # Fixed effects for each site, including the intercept; the coefficient for X2 varies across sites, while others are shared

# Generate site-specific variance component parameters
set.seed(200); theta_site <- rnorm(K, 1, 0.5) %>% round(1) %>% as.list # Variance component parameter for each site

# Generate fixed effect covariates
X <- data.frame(X0 = 1, X1 = rbinom(N, 1, 0.3)) %>% mutate(X2 = rnorm(N), X3 = runif(N), X4 = rbinom(N, 1, 0.5)) %>% as.matrix

# Generate synthetic dataset for testing
test_data <- get_simulated_df(X, patient, site, 1, beta_site, theta_site)

# Perform meta-analysis to obtain initial estimates for the parameters
library(metafor)

# Define GLMM formula
formula <- Y ~ 1 + X1 + X2 + X3 + X4 + (1 | patient)

# Parameter names
par_name <- c("theta1", "Intercept", paste0("X", 1:4))

# Fit GLMM for each site separately
glmm_local <- group_split(test_data, site) %>% map(~glmer(formula, data = ., family = binomial, nAGQ = 1, control = glmerControl("bobyqa", optCtrl = list(maxfun = 1e5))))

# Extract parameter estimates and variances for each site
estimate_local <- glmm_local %>% map(~data.frame(Parameter = par_name, Estimate = c(.@theta, .@beta), Variance = solve(.@optinfo$derivs$Hessian/2) %>% diag) %>% column_to_rownames("Parameter"))

# Combine site-specific estimates for meta-analysis
meta_stat_df <- map(par_name, ~lapply(estimate_local, function(x){x[.,]})) %>% set_names(par_name) %>% 
  map(~do.call(rbind, .)) %>% 
  map(~mutate(., site = unique(site)))

# Perform random-effects meta-analysis for X2 and fixed-effects meta-analysis for other parameters
RE_meta <- map_dbl("X2", ~rma(Estimate, Variance, data = meta_stat_df[[.]], method = "REML", control = list(maxiter = 1e9, stepadj = 0.1))$beta) %>% set_names("X2")
FE_meta <- map_dbl(par_name %>% setdiff("X2"), ~rma(Estimate, Variance, data = meta_stat_df[[.]], method = "FE")$beta) %>% set_names(par_name %>% setdiff("X2"))

# Combine meta-analysis results as initial estimates
meta_estimate <- c(RE_meta, FE_meta) %>% .[par_name]

# Split the initial estimates into variance component, shared fixed effects, and site-specific fixed effects
n_theta <- 1; FE_site_loc <- 3
alpha_prior <- meta_stat_df[n_theta + FE_site_loc] %>% map("Estimate") %>% transpose %>% map(unlist)
theta_prior <- meta_stat_df[1:n_theta] %>% map("Estimate") %>% transpose %>% map(unlist)
beta_prior <- meta_estimate[-(1:n_theta)][-FE_site_loc]

# Create a combined initial value vector for Fed-GLMM
initial_value <- c(unlist(theta_prior), beta_prior, unlist(alpha_prior))


# A few testing examples
initial_value <- c(unlist(theta_prior), beta_prior, unlist(alpha_prior))
test_FedGLMM <- get_FedGLMM_estimates(formula, test_data, "site", 3, initial_value, n_iter = 5, family = binomial(), nAGQ = 1, control = glmerControl("bobyqa", optCtrl = list(maxfun = 1e5)))

initial_value <- c(unlist(theta_prior), meta_estimate[-(1:n_theta)])
test_FedGLMM <- get_FedGLMM_estimates(formula, test_data, "site", NULL, initial_value, n_iter = 5, family = binomial(), nAGQ = 1, control = glmerControl("bobyqa", optCtrl = list(maxfun = 1e5)))

initial_value <- c(meta_estimate[1:n_theta], beta_prior, unlist(alpha_prior))
test_FedGLMM <- get_FedGLMM_estimates(formula, test_data, "site", 3, initial_value, n_iter = 5, common_theta = T, family = binomial(), nAGQ = 1, control = glmerControl("bobyqa", optCtrl = list(maxfun = 1e5)))

initial_value <- meta_estimate
test_FedGLMM <- get_FedGLMM_estimates(formula, test_data, "site", NULL, initial_value, n_iter = 5, common_theta = T, family = binomial(), nAGQ = 1, control = glmerControl("bobyqa", optCtrl = list(maxfun = 1e5)))


# Display final parameter estimates
test_FedGLMM$final_result %>% round(4) # estimates
test_FedGLMM$var_cov %>% diag %>% sqrt %>% round(4) # standard error

# Compare with the pooled analysis
pooled_model <- glmer(formula, data = test_data, family = binomial(), nAGQ = 1, control = glmerControl("bobyqa", optCtrl = list(maxfun = 1e5)))
summary(pooled_model)

