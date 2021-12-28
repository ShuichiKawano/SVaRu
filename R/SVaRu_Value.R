# This program is used to compute the value of GIC in the function "VSM_GIC". 

library(MASS)
source("detailed_calculate.R")

VSM_GIC_Value <- function(res_data, Basis_function, gam_1, gam_2, VSM_coef, reg_var){
  
  n_sample <- length(res_data)
  n_basis <- length(VSM_coef)
  
  R <- matrix(0, nrow=n_basis, ncol = n_basis)
  for (p in 4: (n_basis-4)){
    for (j in 4: (n_basis-4)){
      if (p == j){
        R[p, j] <- w_p_p(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p)
      }else if(p == j + 1){
        R[p, j] <- w_p_p_plus_1(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p)
      }else if(p == j - 1){
        R[p, j] <- w_p_p_minus_1(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p)
      }else if(p == j + 2){
        R[p, j] <- w_p_p_plus_2(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p)
      }else if(p == j - 2){
        R[p, j] <- w_p_p_minus_2(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p)
      }else if(p == j + 3){
        R[p, j] <- w_p_p_plus_3(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p)
      }else if(p == j - 3){
        R[p, j] <- w_p_p_minus_3(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p)
      } else {
        R[p, j] <- log_likelihood(reg_var, Basis_function, p, j)
      }
    }
  }
  
  Q_1 <- matrix(0, nrow=n_basis, ncol = 1)
  Q_2 <- matrix(0, nrow=n_basis, ncol = 1)
  
  for (p in 4: (n_basis-4)){
    Q_1[p] <- w_p(res_data, VSM_coef, reg_var, Basis_function, gam_1, gam_2, p);
  }
  
  reg_var <- sum(reg_var)
  Q_2 <- (1/reg_var) * ( t(Basis_function) %*% res_data -  t(Basis_function) %*%  Basis_function %*% VSM_coef)
  Q <- Q_1 %*% t(Q_2)
  RQ <- ginv(R) *  Q;

  return(RQ)
}
