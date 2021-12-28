### Gaussian basis function
Gauss_Kernel <- function(exp_data, basis_mean_point, width_para_vec){
  n_basis <- length(basis_mean_point)
  Gaussian_kernel <- matrix(0, nrow = length(exp_data), ncol =  n_basis)
  for (j in 1: n_basis){
    for (i in 1: length(exp_data)){
      Gaussian_kernel[i,j] <- exp(-(exp_data[i]-basis_mean_point[j])^(2)/(2 * width_para_vec[j]))
    }
  }
  return(Gaussian_kernel)
}
