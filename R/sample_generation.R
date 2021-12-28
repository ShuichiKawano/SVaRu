### generating artificial data

source("True_Func.R")

sample_generation <- function(n_sample, start_data=-2, end_data=2, true_func_num=6, true_stan_dev=sqrt(0.05)) {
  exp_data  <- seq(start_data, end_data, length=n_sample)
  S_Gene <- True_Func(n_sample, exp_data, true_func_num, true_stan_dev)
  res_data <- S_Gene[[1]]
  true_curve <- S_Gene[[2]]
  return(list(exp_data=exp_data, true_curve=true_curve, res_data=res_data))
}
