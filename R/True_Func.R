### true fucntion

True_Func <- function(n_sample, exp_data, true_func_num, true_stan_dev){
  true_curve <- numeric(n_sample)
  for (i in 1:n_sample){
    if (true_func_num == 1){
      true_curve[i] <- sin(32*exp_data[i]^3)
    }
    if (true_func_num == 2){
      true_curve[i] <- (10/sqrt(2*pi)) * exp(-exp_data[i]^2/2)
    }
    # Doppler
    if (true_func_num == 3){
      ips = 0.1
      true_curve[i] <- (exp_data[i] * (1 - exp_data[i])) ^ (1/4) * sin((2 * pi * (1 + ips))/(exp_data[i] + ips))
    }
    if (true_func_num == 4){ 
      if ((0 <= exp_data[i]) && (exp_data[i] <= 0.5)){
        true_curve[i] <- - sin(2* pi*exp_data[i])
      }else if((0.5 < exp_data[i]) && (exp_data[i] <= 1)){
        true_curve[i] <-  sin(15*pi*exp_data[i])
      }
    }
    # Sinc
    if (true_func_num == 5){
      if(exp_data(i)==0){
        true_curve[i] <- 1
      }else{
        true_curve[i] <- sin(exp_data[i]) / exp_data[i]
      }
    }
    if (true_func_num == 6){
      true_curve[i] <- sin(exp_data[i]) + 2 * exp(-30 * exp_data[i]^2)
    }
    
  }
  
  res_data <- numeric(n_sample)
  for (i in 1:n_sample){
    res_data[i] <- true_curve[i] + rnorm(1, mean=0, sd = true_stan_dev)
  }
  return(list(res_data,true_curve))
}
