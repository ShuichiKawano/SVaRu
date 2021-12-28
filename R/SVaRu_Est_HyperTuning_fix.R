# This program is used to estimate parameters and tuning parameters in VSM.

VSM_Est_HyperTuning_fix <- function(Basis_function, res_data, gam_1, gam_2, opt_GIC_smooth){
  
  n_sample <- dim(Basis_function)[1]
  n_basis <- dim(Basis_function)[2]
  VSM_Coef <- seq(1,1,length = n_basis)
  VSM_Lam <- opt_GIC_smooth *  seq(1,1,length = n_basis)
  
  s <- 1
  while (s > 10^(-5)){
    VSM_Coef = solve(t(Basis_function) %*% Basis_function + 2 * diag(VSM_Lam)) %*% t(Basis_function) %*% res_data
    VSM_Lam_old = VSM_Lam
    for (j in 1: n_basis){
      if (j == 1){
        b_1 <- ( VSM_Coef[j]^(2) -2 * gam_1 * VSM_Lam[j+1] )
        c_1 <-  b_1^(2) + 16 * gam_1 * gam_2;
        VSM_Lam[j] <-  (-b_1 + sqrt( c_1 )) / (8 * gam_1) 
      }else if ((1 < j) && (j < n_basis)){
        b_2 <- ( VSM_Coef[j]^(2) - 2 * gam_1 * VSM_Lam[j+1] - 2 * gam_1 * VSM_Lam[j-1] );
        c_2 <-  b_2^(2) + 16 * gam_1 * gam_2;
        VSM_Lam[j] <-  (-b_2 + sqrt( c_2 )) / (8 * gam_1) 
      }else if (j == n_basis){
        b_3 <- ( VSM_Coef[j]^(2) -2 * gam_1 * VSM_Lam[j-1])
        c_3 <-  b_3^(2) + 16 * gam_1 * gam_2;
        VSM_Lam[j] <-  (-b_3 + sqrt( c_3 )) / (8 * gam_1) 
      }
      s = t(VSM_Lam_old - VSM_Lam) %*% (VSM_Lam_old - VSM_Lam);
      if(s > 10^(3)){
        break
      }
    }
  }
  VSM_Coef <- as.vector(VSM_Coef)
  VSM_Lam <- as.vector(VSM_Lam)
  
  return(list(VSM_Coef, VSM_Lam))
}
