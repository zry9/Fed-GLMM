library(dplyr)
library(purrr)
library(tibble)
library(lme4)
library(numDeriv)

#Helper function 1
get_glmer_derivative <- function(formula, data, beta_prior, theta, ...){
  
  dglmm_subset_t1 <- Sys.time()
  
  devfun <- glmer(formula, data, devFunOnly = T, ...)
  
  gradient <- grad(devfun, c(theta, beta_prior)); Hessian <- hessian(devfun, c(theta, beta_prior))
  
  return(list(gradient = gradient, Hessian = Hessian))
  
}

#Helper function 2
get_FE_site_rearranged <- function(deriv, FE_site_loc, n_theta){
  
  FE_site_loc_full <- FE_site_loc + n_theta
  
  if(is.vector(deriv)){
    deriv[-FE_site_loc_full] %>% c(deriv[FE_site_loc_full])
  } else {
    deriv_tranCol <- deriv[,-FE_site_loc_full] %>% cbind(deriv[,FE_site_loc_full]) 
    deriv_tranCol[-FE_site_loc_full,] %>% rbind(deriv_tranCol[FE_site_loc_full,]) 
  }
}


#Helper function 3
get_one_iteration <- function(formula, data, site_var, FE_site_loc, estimate_prior, ...){
  
  
  site <- data[[site_var]]
  K <- length(unique(site))
  n_alpha <- length(FE_site_loc)
  n_theta <- glFormula(formula, data, ...)$reTrms$theta %>% length
  n_beta <- length(estimate_prior) - K*n_theta - K*n_alpha 
  
  theta_list <- estimate_prior[(1): (K*n_theta)] %>% split(ceiling(seq_along(.)/n_theta))
  beta_shared <- estimate_prior[(K*n_theta + 1): (K*n_theta + n_beta)]
  beta_site_list <- estimate_prior[(K*n_theta + n_beta + 1): (K*n_theta + n_beta + K*n_alpha)] %>% split(ceiling(seq_along(.)/n_alpha))
  
  
  site_data <- map(unique(site), ~data[site == .,])
  
  
  beta_list <- list(beta_shared) %>% rep(length(beta_site_list))
  for(i in 1:length(FE_site_loc)){
    beta_list <- map2(beta_list, beta_site_list, ~append(.x, .y[i], after = FE_site_loc[i]-1))
  }
  
  
  derivative_list_full <- derivative_list <- glmer_derivative <- pmap(list(site_data, beta_list, theta_list), possibly(function(A, B, C){get_glmer_derivative(formula, A, B, C, ...)}, otherwise = NULL))
  
  
  
  if(length(compact(glmer_derivative)) < length(unique(site))){
    dev_gradient_avg <- map(derivative_list %>% compact, ~.$gradient) %>% Reduce("+", .) %>% `/`(length(derivative_list %>% compact))  
    dev_Hessian_avg <- map(derivative_list %>% compact, ~.$Hessian) %>% Reduce("+", .) %>% `/`(length(derivative_list %>% compact)) 
    derivative_list_full <- derivative_list %>% map(function(x){if(is.null(x)){list(gradient = dev_gradient_avg, Hessian = dev_Hessian_avg)} else {x}})
  }
  
  
  gradient_list_rearranged <- derivative_list_full %>% map(~.$gradient/(-2)) %>% map(~get_FE_site_rearranged(.,FE_site_loc, n_theta))
  
  Hessian_list_rearranged <- derivative_list_full %>% map(~.$Hessian/(-2)) %>% map(~get_FE_site_rearranged(.,FE_site_loc, n_theta))
  
  
  gradient_theta <- gradient_list_rearranged %>% map(~.[1:n_theta]) %>% unlist
  gradient_beta <- gradient_list_rearranged %>% map(~.[(n_theta + 1):(n_theta + n_beta)]) %>% Reduce("+", .)
  gradient_alpha <- gradient_list_rearranged %>% map(~.[(n_theta + n_beta + 1):(n_theta + n_beta + n_alpha)]) %>% unlist
  gradient_t <- c(gradient_theta, gradient_beta, gradient_alpha)
  
  gradient_t <- c(gradient_theta, gradient_beta, gradient_alpha)
  
  
  Hessian_t <- matrix(0, K*n_theta + n_beta + K*n_alpha, K*n_theta + n_beta + K*n_alpha) %>% 
    as(Class = "sparseMatrix")
  
  for(k in 1:K){
    Hessian_site <- Hessian_list_rearranged[[k]]
    
    Hessian_t[((k-1)*n_theta+1):(k*n_theta), ((k-1)*n_theta+1):(k*n_theta)] <- Hessian_site[1:n_theta, 1:n_theta]
    Hessian_t[((k-1)*n_theta+1):(k*n_theta), (K*n_theta+1):(K*n_theta+n_beta)] <- Hessian_site[1:n_theta, (n_theta+1):(n_theta+n_beta)]
    Hessian_t[((k-1)*n_theta+1):(k*n_theta), (K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha)] <- Hessian_site[1:n_theta, (n_theta+ n_beta+1):(n_theta+ n_beta + n_alpha)]
    Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), ((k-1)*n_theta+1):(k*n_theta)] <- Hessian_site[(n_theta+1):(n_theta+ n_beta), 1:n_theta]
    Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), (K*n_theta+1):(K*n_theta+n_beta)] <- Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), (K*n_theta+1):(K*n_theta+n_beta)] + Hessian_site[(n_theta+1):(n_theta+ n_beta), (n_theta+1):(n_theta+ n_beta)]
    Hessian_t[(K*n_theta+1):(K*n_theta+n_beta), (K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha)] <- Hessian_site[(n_theta+1):(n_theta+ n_beta), (n_theta+ n_beta+1):(n_theta+ n_beta + n_alpha)]
    Hessian_t[(K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha), ((k-1)*n_theta+1):(k*n_theta)] <- Hessian_site[(n_theta+ n_beta+1):(n_theta+ n_beta+ n_alpha), 1:n_theta]
    Hessian_t[(K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha), (K*n_theta+1):(K*n_theta+n_beta)] <- Hessian_site[(n_theta+ n_beta+1):(n_theta+ n_beta+ n_alpha), (n_theta+1):(n_theta+ n_beta)]
    Hessian_t[(K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha), (K*n_theta+n_beta+(k-1)*n_alpha+1):(K*n_theta+n_beta+k*n_alpha)] <- Hessian_site[(n_theta+ n_beta+1):(n_theta+ n_beta+ n_alpha), (n_theta+ n_beta+1):(n_theta+ n_beta + n_alpha)]
    
  }
  
  estimate_bar <- estimate_prior - solve(Hessian_t) %*% (gradient_t)
  
  
  return(list(estimate_bar = estimate_bar, gradient_t = gradient_t, Hessian_t = Hessian_t))
  
}


#Fed-GLMM
get_FedGLMM_estimates <- function(formula, data, site_var, FE_site_loc, initial_values, n_iter = 5, ...){
  
  estimate_history <- list() 
  diff_estimate_bar_history <- list()
  gradient_history <- list()
  Hessian_history <- list()
  
  
  
  estimate_history[[1]] <- final_result <- estimate_prior <- initial_values 
  
  
  iter_count <- 0
  
  for (x in 1:n_iter){
    
    converge_mark <- F
    converge_iter <- NaN
    diverge_mark <- F
    
    writeLines(paste0("#### iter ", x, " ####\n"))
    
    iteration <- get_one_iteration(formula, data, site_var, FE_site_loc, estimate_prior,...)
    
    diff_estimate_bar <- diff_estimate_bar_history[[x]] <- sqrt(sum((estimate_prior - iteration$estimate_bar)^2))
    
    writeLines(paste0("Beta Difference: ", diff_estimate_bar))
    
    
    final_result <- iteration$estimate_bar
    estimate_history[[x+1]] <- iteration$estimate_bar 
    gradient_history[[x]] <- iteration$gradient_t
    Hessian_history[[x]] <- iteration$Hessian_t
    
    iter_count <- x
    
    if (diff_estimate_bar > 200) {
      writeLines("The values have diverged.\n")
      diverge_mark <- T
      break
    }
    
    if(diff_estimate_bar < 1e-6){
      writeLines("The values have converged.\n")
      converge_mark <- T
      converge_iter <- x
      break
    } else {
      estimate_prior <- iteration$estimate_bar
    }
    
    
  }
  return(list(final_result = final_result, estimate_history = estimate_history, gradient_history = gradient_history, Hessian_history = Hessian_history, diff_estimate_bar_history = diff_estimate_bar_history, iter_count = iter_count, converge_mark = converge_mark, converge_iter = converge_iter, diverge_mark = diverge_mark))

}


#Data Simulator
get_covmatrix_block <- function(theta){
  
  nc <- (sqrt(length(theta)*8 + 1) - 1)/2 
  rowIndices <- rep(1:nc, 1:nc)
  colIndices <- sequence(1:nc)
  templatet <- sparseMatrix(rowIndices, colIndices, x = 1 * (rowIndices == colIndices)) %>% t 
  LindTemplate <- rowIndices + nc * (colIndices - 1) - choose(colIndices, 2) 
  
  templatet@x <- theta[LindTemplate]
  D_block <- t(templatet) %*% templatet 
  return(D_block)
  
}


get_sigmoid <- function(x){
  1/(1+exp(-x))
}

get_simulated_df <- function(Xall, patient, site, RE_col_ind, beta_true_list, theta_list){
  
  J_k <- map(unique(site), ~t(as(factor(patient)[site == .], Class = "sparseMatrix"))) 
  X_k <- map(unique(site), ~Xall[site == .,])
  X_RE_k <- map(unique(site), ~Xall[site == .,RE_col_ind]) 
  Z_k <- map2(J_k, X_RE_k, ~t(KhatriRao(t(.x), t(.y))))
  
  
  D_k <- theta_list %>% map(get_covmatrix_block) %>% map2(unique(site), ~.bdiag(rep(list(.x), length(unique(patient[site == .y])))))
  
  b_k <- map(D_k, ~t(chol(.)) %*% rnorm(nrow(.)))
  
  ita_k <- pmap(list(X = X_k, Z = Z_k, b = b_k, beta_true = beta_true_list), ~with(list(...), X %*% beta_true + Z %*% b))
  ita <- do.call(rbind, ita_k) %>% get_sigmoid %>% as.numeric
  y <- rbinom(length(site), 1, ita)
  
  df <- cbind(y, Xall, patient, site, ita) %>% as.matrix %>% as.data.frame %>% rename(Y = 1, X0 = 2)
  
  return(df)
}


#Simulate test data
K <- 8 #site
n <- 4 #visit per patient
m <- seq(from = 100, to = 50*(K+1), by = 50) #patient per site
N <- sum(m)*n # total obs
M <- sum(m) #total patients
site <- rep(1:K, m*n) #site index for each obs
patient <- rep(1:M, each = n) #patient index for each obs

set.seed(100); beta_site <- map(rnorm(K, 0.5, 0.5) %>% round(1), ~c(1.5, 1, ., 0.5, 0.5))
set.seed(200); theta_site <- rnorm(K, 1, 0.5) %>% round(1) %>% as.list


X <- data.frame(X0 = 1, X1 = rbinom(N,1,0.3)) %>% mutate(X2 = rnorm(N), X3 = runif(N), X4 = rbinom(N,1,0.5)) %>% as.matrix

test_data <- get_simulated_df(X, patient, site, 1, beta_site, theta_site)


#Meta-analysis results to be used as initial value
library(metafor)

formula <- Y ~ 1 + X1 + X2 + X3 + X4 + (1 | patient)

par_name <- c("theta1", "Intercept", paste0("X", 1:4))

glmm_local <- group_split(test_data, site) %>% map(~glmer(formula, data = ., family = binomial, nAGQ = 1, control=glmerControl("bobyqa", optCtrl=list(maxfun=1e5))))

estimate_local <- glmm_local %>% map(~data.frame(Parameter = par_name, Estimate = c(.@theta, .@beta), Variance = solve(.@optinfo$derivs$Hessian/2) %>% diag) %>% column_to_rownames("Parameter"))


meta_stat_df <- map(par_name, ~lapply(estimate_local, function(x){x[.,]})) %>% set_names(par_name) %>% 
  map(~do.call(rbind, .)) %>% 
  map(~mutate(., site = unique(site))) 

RE_meta <- map_dbl("X2", ~rma(Estimate, Variance, data = meta_stat_df[[.]], method = "REML", control=list(maxiter=1e9, stepadj=0.1))$beta) %>% set_names("X2")
FE_meta <- map_dbl(par_name %>% setdiff("X2"), ~rma(Estimate, Variance, data = meta_stat_df[[.]], method = "FE")$beta) %>% set_names(par_name %>% setdiff("X2"))

meta_estimate <- c(RE_meta, FE_meta) %>% .[par_name] 

n_theta <- 1; FE_site_loc <- 3

alpha_prior <- meta_stat_df[n_theta+FE_site_loc] %>% map("Estimate") %>% transpose %>% map(unlist)
theta_prior <- meta_stat_df[1:n_theta] %>% map("Estimate") %>% transpose %>% map(unlist)
beta_prior <- meta_estimate[-(1:n_theta)][-FE_site_loc]

initial_value <- c(unlist(theta_prior), beta_prior, unlist(alpha_prior))


#Perform Fed-GLMM on simulated data
test_FedGLMM <- get_FedGLMM_estimates(formula, test_data, "site", 3, initial_value, n_iter = 5, family = binomial(), nAGQ = 1, control=glmerControl("bobyqa", optCtrl=list(maxfun=1e5)))

test_FedGLMM$final_result %>% round(1) %>% as.vector()
