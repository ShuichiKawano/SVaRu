source("SVaRu.R")

#library(MASS)


source("sample_generation.R")
#source("True_Func.R")
#source("Gauss_Kernel.R")

######6
#true_func_num <- 6  #真の関数
#start_data=-2  	#データ_説明変数の発生開始
#end_data=2		#データ_説明変数の発生終了

###調整パラメータ###

#monte_num <- 1 		#Monte-Carlo回数
n_sample <- 200  	#Sample_size
#true_stan_dev <- sqrt(0.05)  #真の誤差標準偏差

basis_num_vec <- 40	#基底関数の数_候補ベクトル
width_para <- 10^(-1.5)	#基底関数のバンド幅

#start_smooth <- 10^(-9)	#ridge平滑化パラメータ
#end_smooth <- 10^(-1)	#
#n_smooth <- 7		#平滑化パラメータ数

start_gamma_1 <- exp(-5)#提案手法ハイパーチューニングパラメータ
end_gamma_1 <- exp(-1)
start_gamma_2 <- exp(-5)
end_gamma_2 <- exp(-1)
n_gamma <- 10		
#####調整パラメータ終了#####

#smooth_vec <- seq(start_smooth, end_smooth, length = n_smooth)
gam_1_vec <- seq(start_gamma_1, end_gamma_1, length = n_gamma)
gam_2_vec <- seq(start_gamma_2, end_gamma_2, length = n_gamma)

sample <- sample_generation(n_sample)
exp_data <- sample$exp_data
res_data <- sample$res_data
#VSM
SVaRu(exp_data, res_data, basis_num_vec, width_para, gam_1_vec, gam_2_vec)
