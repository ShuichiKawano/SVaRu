### SVaRu: This program performs smoothing varying regularization with generalized information criterion (GIC).

### INPUT
# exp_data: a vector of explanatory variable
# res_data: a vector of response variable
# basis_num_vec: candidates of the number of basis functions
# width_para: width of basis functions
# gam_1_vec: candidates of the hyper-tuning parameter gamma_1
# gam_2_vec: candidates of the hyper-tuning parameter gamma_2

### OUTPUT
# opt_n_basis: the optimal number of basis functions selected by GIC
# opt_gam_1: the optimal value of the hyper-tuning parameter gamma_1 selected by GIC
# opt_gam_2: the optimal value of the hyper-tuning parameter gamma_2 selected by GIC
# VSM_opt_coef: estimate of coefficients
# VSM_reg_var: estimate of error variance
# VSM_opt_lam: estimate of smoothing parameters

source("Gauss_Kernel.R")
source("SVaRu_Value.R")
source("SVaRu_Est_HyperTuning_fix.R")

SVaRu <- function(exp_data, res_data, basis_num_vec, width_para, gam_1_vec, gam_2_vec){
  
  n_sample <- length(exp_data)
  n_gamma <- length(gam_1_vec) 
  VSM_GIC_Mat <- matrix (0, nrow = n_gamma, ncol= (length(basis_num_vec) * n_gamma))
  for (k in 1:length(basis_num_vec)){
    n_basis <- basis_num_vec[k];
    width_para_vec <- width_para * seq(1,1,length = n_basis)
    basis_mean_point <- seq(min(exp_data), max(exp_data),length = n_basis)
    Basis_function <- Gauss_Kernel(exp_data, basis_mean_point, width_para_vec)
    opt_GIC_smooth <- 10^(-1)
    for (i in 1: n_gamma){
      for (j in 1: n_gamma){
        gam_1 = gam_1_vec[i]
        gam_2 = gam_2_vec[j]
        VSM_Est_paras <- VSM_Est_HyperTuning_fix(Basis_function, res_data, gam_1, gam_2, opt_GIC_smooth)
        VSM_coef <- VSM_Est_paras[[1]]
        VSM_coef <- as.vector(VSM_coef)
        VSM_Lam <- VSM_Est_paras[[2]]
        VSM_Lam <- as.vector(VSM_Lam)
        est_VSM <- Basis_function %*% VSM_coef
        reg_var <- (1/n_sample) * t(res_data - est_VSM) %*% (res_data - est_VSM)
        GIC_RQ <- VSM_GIC_Value(res_data, Basis_function, gam_1, gam_2, VSM_coef, reg_var)
        VSM_GIC_Mat[i, j + (k-1) * length(gam_1_vec)] <- n_sample * (log(2 * pi) + 1) + n_sample * log(reg_var) + 2 * sum(diag(GIC_RQ))
      }
    }
  }
  
  min_GIC_VSM = min(min(min(VSM_GIC_Mat)));
  for (k in 1: length(basis_num_vec)){
    for (i in 1:n_gamma){
      for (j in 1:n_gamma){
        if(VSM_GIC_Mat[i, j + (k-1) * n_gamma] == min_GIC_VSM){
          k_0 <- k
          i_0 <- i
          j_0 <- j
        }
      }
    }
  }
  opt_n_basis = basis_num_vec[k_0]
  opt_gam_1 = gam_1_vec[i_0]
  opt_gam_2 = gam_2_vec[j_0]
  
  VSM_basis_mean_point <- seq(min(exp_data), max(exp_data),length = opt_n_basis)
  VSM_opt_Basis_function <- Gauss_Kernel(exp_data, VSM_basis_mean_point, width_para_vec)
  VSM_opt_est_paras <- VSM_Est_HyperTuning_fix(VSM_opt_Basis_function, res_data, opt_gam_1, opt_gam_2, opt_GIC_smooth)
  VSM_opt_coef <- VSM_opt_est_paras[[1]]
  VSM_opt_lam <- VSM_opt_est_paras[[2]]
  VSM_reg_var <- (1/n_sample) * t(res_data - VSM_opt_Basis_function %*% VSM_opt_coef) %*% (res_data - VSM_opt_Basis_function %*% VSM_opt_coef)
  return(list(opt_n_basis=opt_n_basis, opt_gam_1=opt_gam_1, opt_gam_2=opt_gam_2, VSM_opt_coef=VSM_opt_coef, VSM_reg_var=VSM_reg_var, VSM_opt_lam=VSM_opt_lam))
}
