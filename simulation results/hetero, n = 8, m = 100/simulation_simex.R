setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
remove(list = ls())
source("function_simu.R")

lambda_arr = c(seq(0,1,by = 0.05), seq(1,1.5,by = 0.05))

n_simu = 20
result_mat = matrix(nrow = n_simu, ncol = 9)
colnames(result_mat) <- sapply(c("actual", "naive", "simex"), function(x) 
  sapply(c("intercept", "w", "z"), function(y) paste0(x,"_",y)))
result_all = list()

start = Sys.time()
for(i in 1:n_simu){
  result = simu(m = 100, n = 8, x.model = "heterogeneous", seed = i, 
                lambda_arr = lambda_arr, B = 100, time.check = FALSE)
  result_mat[i,] = unlist(result$coeffs)   
  result_all = c(result_all, result)
  if(i == 1) message("total estimated time = ", round(lubridate::seconds_to_period(n_simu * result$elapsed.time)))
  cat("Done", i/n_simu * 100, "%", " ... estimated time remaining = ", paste0(round(lubridate::seconds_to_period((n_simu-i) * result$elapsed.time))))
  cat("\r")
}
message("Actual time:", Sys.time()-start)

saveRDS(result_all, file = "result_all.RDS")
saveRDS(result_mat, file = "result_mat.RDS")

colMeans(result_mat[,4:9])
mse_simex = colMeans((result_mat[,7:9]-result_mat[,1:3])^2)
mse_naive = colMeans((result_mat[,4:6]-result_mat[,1:3])^2)
improvement = round(mse_simex/mse_naive * 100)
improvement
plot(result$simex.output)
