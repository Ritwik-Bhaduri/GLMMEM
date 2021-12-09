setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
remove(list = ls())
source("function_simu.R")
library("SurvDisc")

lambda_arr = c(seq(0,1,by = 0.05), seq(1,1.5,by = 0.05))

n_simu = 7
result_mat = matrix(nrow = n_simu, ncol = 9)
colnames(result_mat) <- sapply(c("actual", "naive", "simex"), function(x) 
  sapply(c("intercept", "w", "z"), function(y) paste0(x,"_",y)))
result_all = list()

######################################################################################################
simu_apply <- function(seed, m = 100, n = 8, x.model = "heterogeneous", 
                       lambda_arr = lambda_arr, B = 100, time.check = FALSE){
  simu <- function(m = 50, n = 3, x.model = "homogeneous", seed = 1, 
                   lambda_arr = seq(0,1,by = 0.05), B = 100, time.check = FALSE, plot = FALSE, message = FALSE){
    set.seed(seed)
    
    #parameters
    mu_x = 0
    sigma_x_sq = 1
    sigma_x_u_sq = 1.5
    sigma_u_sq = 0.5
    sigma_y_sq = 0.5
    beta_0 = 0
    beta_x = 2
    beta_z = 1
    theta = 0.5
    
    cluster = rep(1:m, each = n)
    z = rnorm(m*n, 0, 1)
    A  = matrix( rep(diag(m), each = n ) , ncol = m , byrow = FALSE)
    
    if(x.model == "homogeneous"){
      x = mu_x + rnorm(m*n, 0, sigma_x_sq) # for homogeneous
    }else{
      a = rnorm(m, 0, sigma_x_u_sq) # for heterogeneous
      x = mu_x + rep(a, each = n) + rnorm(m*n, 0, sigma_x_sq) # for heterogeneous
    }
    
    w = x + rnorm(m*n, 0, sigma_u_sq)
    b_i = rnorm(m, 0, theta)
    beta_fixed = c("intercept" = beta_0, "beta_x" = beta_x, "beta_z" = beta_z)
    eta = beta_0 + beta_x * x + beta_z * z + rep(b_i, each = n) 
    mu = eta
    y = sapply(mu, function(p) rnorm(1, p, sigma_y_sq))
    # A = rep(as.factor(1:m),each = n)
    data_sim = data.frame(cluster,y, w, x, z)
    
    if(message == TRUE){
      t_one = 0.01
      message("estimated time: ", round(lubridate::seconds_to_period(t_one * length(lambda_arr)*B)))
      if(time.check == TRUE){
        message("Do you want to proceed?")
        if(askYesNo("Do you want to proceed?", default = TRUE, 
                    prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel")))) == FALSE) return(NULL)
      }
      message("Running...")
    }
    
    start_time = Sys.time() 
    fm2=nlme::lme.formula(fixed = y ~ w + z,
                          data = data_sim, random = ~ 1|cluster,
                          control=nlme::lmeControl(returnObject=TRUE))
    (s1 = simexlme(model=fm2, model.model=data_sim[,c("y","cluster","w","z")],
                   SIMEXvariable="w",respvar="y",grpvar="cluster",corform="~1 | cluster",
                   measurement.error=sigma_u_sq,measurement.error.resp=sigma_y_sq,
                   lambda = lambda_arr,B = B, fitting.method = "quadratic",
                   jackknife.estimation = FALSE))
    elapsed_time = Sys.time()-start_time
    if(message == TRUE) message("Actual time:", Sys.time()-start_time )
    
    coeffs <- list("actual" = c(beta_0,beta_x,beta_z), "naive" = fm2$coefficients$fixed, "simex" = s1$coefficients)
    if(plot == TRUE) plot(s1)
    
    return(list("coeffs" = coeffs, "elapsed.time" = elapsed_time, "simex.output" = s1, "naive.lme.output" = fm2))
  }
  result = simu(m = m, n = n, x.model = x.model, seed = seed, 
                lambda_arr = lambda_arr, B = B, time.check = time.check)
  result_mat[seed,] = unlist(result$coeffs)   
  result_all = c(result_all, result)
}

library(parallel)
ncores <- detectCores()
cl <- makeCluster(ncores-1)
clusterApply(cl, x = 1:n_simu, fun = simu_apply)

# start = Sys.time()
# 
# for(i in 1:n_simu){
#   result = simu(m = 100, n = 8, x.model = "heterogeneous", seed = i, 
#                 lambda_arr = lambda_arr, B = 100, time.check = FALSE)
#   result_mat[i,] = unlist(result$coeffs)   
#   result_all = c(result_all, result)
#   if(i == 1) message("total estimated time = ", round(lubridate::seconds_to_period(n_simu * result$elapsed.time)))
#   cat("Done", i/n_simu * 100, "%", " ... estimated time remaining = ", paste0(round(lubridate::seconds_to_period((n_simu-i) * result$elapsed.time))))
#   cat("\r")
# }
# message("Actual time:", Sys.time()-start)
# 
# saveRDS(result_all, file = "result_all.RDS")
# saveRDS(result_mat, file = "result_mat.RDS")

colMeans(result_mat[,4:9])
mse_simex = colMeans((result_mat[,7:9]-result_mat[,1:3])^2)
mse_naive = colMeans((result_mat[,4:6]-result_mat[,1:3])^2)
improvement = round(mse_simex/mse_naive * 100)
improvement
plot(result$simex.output)
