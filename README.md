# Smoothing varying regularization.
This is an R source code for performing smoothing varying regularization (SVaRu). The directory *R* includes many files that perform SVaRu. Among them, most important file is given as follows.
- **SVaRu.R** provides SVaRu with the values of tuning parameters selected by generalized information criterion. 

SVaRu is introduced in the paper:
Kim, D., Kawano, S. and Ninomiya, Y. (2021) Smoothly varying regularization. arXiv:2102.00136 [URL](https://arxiv.org/abs/2102.00136).

## Usage example
Read source files.
```
source("R/SVaRu.R")
```

Setting of simulation.
```
source("R/sample_generation.R")
n_sample <- 200 # sample size
sample <- sample_generation(n_sample)
exp_data <- sample$exp_data
res_data <- sample$res_data
basis_num_vec <- 40 # candidates of the number of basis functions
width_para <- 10^(-1.5) # width of basis functions
start_gamma_1 <- exp(-5) # minimum value of candidates of the hyper-tuning parameter gamma_1
end_gamma_1 <- exp(-1) # maximum value of candidates of the hyper-tuning parameter gamma_1
start_gamma_2 <- exp(-5)  # minimum value of candidates of the hyper-tuning parameter gamma_2
end_gamma_2 <- exp(-1) # maximum value of candidates of the hyper-tuning parameter gamma_2
n_gamma <- 10 # the number of candidates of the hyper-tuning parameters
gam_1_vec <- seq(start_gamma_1, end_gamma_1, length = n_gamma)
gam_2_vec <- seq(start_gamma_2, end_gamma_2, length = n_gamma)
```

Perform SVaRu.
```
# Perform SVaRu in the file SVaRu.R
SVaRu(exp_data, res_data, basis_num_vec, width_para, gam_1_vec, gam_2_vec)
```
