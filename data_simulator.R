library(dplyr)  # Data manipulation
library(purrr)  # Functional programming tools
library(tibble) # Data frames

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


