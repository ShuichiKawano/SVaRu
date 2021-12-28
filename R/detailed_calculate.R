# This program is used in the function "VSM_GIC_Value". 

####################################
lambda <- function(VSM_coef, p, gam_1, gam_2){
  Value <- gam_2 * (VSM_coef[p]^2 + 2 * gam_1 * gam_2 * (VSM_coef[p-1]^(-2) - 2 * VSM_coef[p]^(-2) + VSM_coef[p+1]^(-2)))^(-1)
  return(Value)
}

####################################
log_likelihood <- function(reg_var, Basis_function, p1, p2){
  Value <- {
    Value <- -(1/reg_var) * Basis_function[p1, ] %*% Basis_function[p2, ]
  }
  return(Value)
}

####################################
Func_eta <-function(VSM_coef, p, gam_1, gam_2){
  Value <- VSM_coef[p]^2 + 2 * gam_1 * gam_2 * (VSM_coef[p-1]^(-2) - 2 * VSM_coef[p]^(-2) + VSM_coef[p+1]^(-2))
  return(Value)
}

####################################
w_p <- function(res_data, VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  
  aa <- 4 * gam_1 * gam_2^2 * VSM_coef[p]^(3)*Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2)
  bb <- -2 * gam_2 *Func_eta(VSM_coef, p, gam_1, gam_2)^(-2) * (VSM_coef[p] + 4 * gam_1 *gam_2 * VSM_coef[p]^(-3))
  cc <- 4 * gam_1 * gam_2^2 * VSM_coef[p]^(3)*Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2)
  
  Value <- {
    (1/reg_var) * sum( (res_data - Basis_function[p,] %*% VSM_coef) %*% t(Basis_function[p,]) )
    -(1/2) * (aa * VSM_coef[p-1]^(2) + bb * VSM_coef[p]^(2) + 2 * lambda(VSM_coef, p, gam_1, gam_2) * VSM_coef[p] + cc *  VSM_coef[p+1]^(2))
    -(gam_1 / 2) * (-2 * cc * (lambda(VSM_coef, p+2, gam_1, gam_2) -lambda(VSM_coef, p+1, gam_1, gam_2))
                    -2 * (cc * lambda(VSM_coef, p, gam_1, gam_2) + lambda(VSM_coef, p+1, gam_1, gam_2) * bb)
                    +4 * lambda(VSM_coef, p, gam_1, gam_2) * bb 
                    -2 * (bb * lambda(VSM_coef, p-1, gam_1, gam_2) + lambda(VSM_coef, p, gam_1, gam_2) * aa)
                    +2 * lambda(VSM_coef, p-1, gam_1, gam_2) * aa)
    +(gam_2/2) * (aa * lambda(VSM_coef, p-1, gam_1, gam_2)^(-1) 
                  + bb * lambda(VSM_coef, p, gam_1, gam_2)^(-1) 
                  + cc * lambda(VSM_coef, p+1, gam_1, gam_2)^(-1))
  }
  return(Value)
}

####################################
w_p_p <- function(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  
  a4 <- 4*gam_1*gam_2^(2) * (-3 * VSM_coef[p]^(-4)*Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2) 
                             - 8 *gam_1 * gam_2 * VSM_coef[p]^(-6)*Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-5))
  b4 <- -2 * gam_2 * (4 *Func_eta(VSM_coef, p, gam_1, gam_2)^(-5) * (VSM_coef[p] + 4*gam_1 *gam_2 * VSM_coef[p]^(-3) )^2 
                      + Func_eta(VSM_coef, p, gam_1, gam_2)^(-2) * (1+12*gam_1*gam_2 *VSM_coef[p]^(-4)))
  c4 <- -2 * gam_2 * Func_eta(VSM_coef, p, gam_1, gam_2)^(-2) * (VSM_coef[p] - 4 *gam_1*gam_2 * VSM_coef[p]^(-3))
  d4 <- 4 * gam_1 *gam_2^2 * (-3 * VSM_coef[p]^(-4) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2) 
                              -8 * gam_1 * gam_2 *  VSM_coef[p]^(-6) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-5))
  e4 <- 4 * gam_1 *gam_2^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2)
  f4 <- a4
  g4 <- 4 *gam_1 *gam_2^2 * VSM_coef[p]^(-3) *Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2)
                      

  Value <- 
    { log_likelihood(reg_var, Basis_function, p, p)
      -(1/2)*(a4 * VSM_coef[p-1]^(2) + b4 * VSM_coef[p]^2 + 2*c4 * VSM_coef[p] 
            + 2 * c4 * VSM_coef[p] + 2 * lambda(VSM_coef, p, gam_1, gam_2) + d4 * lambda(VSM_coef, p+1, gam_1, gam_2)^2)
      -(gam_1/2)*(- 2 * ( d4 * (lambda(VSM_coef, p+2, gam_1, gam_2) - lambda(VSM_coef, p+1, gam_1, gam_2)) - e4^2 ) 
                - 2 * (f4 * (lambda(VSM_coef, p-2, gam_1, gam_2) - lambda(VSM_coef, p-1, gam_1, gam_2)) - g4 ^2)
                + 2 * (e4^2 + lambda(VSM_coef, p+1, gam_1, gam_2) * d4)
                - 2 * (d4 * lambda(VSM_coef, p, gam_1, gam_2) + 2 * e4 * c4 + lambda(VSM_coef, p+1, gam_1, gam_2) * b4)
                + 4 * (c4^2 + lambda(VSM_coef, p, gam_1, gam_2) * b4)
                - 2 * (b4 * lambda(VSM_coef, p-1, gam_1, gam_2) + 2 * c4 * g4 + lambda(VSM_coef, p, gam_1, gam_2) * a4) 
                + 2 * (g4^2 + lambda(VSM_coef, p-1, gam_1, gam_2) * a4))
      +(gam_2/2) * (a4*lambda(VSM_coef, p-1, gam_1, gam_2)^(-1) - g4^2 * lambda(VSM_coef, p-1, gam_1, gam_2)^(-2) + b4 * lambda(VSM_coef, p, gam_1, gam_2)^(-1) 
                  - c4^2 * lambda(VSM_coef, p, gam_1, gam_2)^(-2) + d4 *lambda(VSM_coef, p+1, gam_1, gam_2)^(-1) - e4^2 * lambda(VSM_coef, p+1, gam_1, gam_2)^(-2) )}
  
  return(Value)
}

####################################
w_p_p_minus_1 <- function(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  
  a2 <- -16*gam_1*gam_2^(2)*Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-3)*VSM_coef[p]^(-3)*(VSM_coef[p-1] - 4*gam_1*gam_2*VSM_coef[p-1]^(-3))
  b2 <- 4 * gam_1 * gam_2 ^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2)
  c2 <- -16*gam_1*gam_2*Func_eta(VSM_coef, p, gam_1, gam_2)^(-3)*VSM_coef[p-1]^(-3)*(VSM_coef[p] + 4*gam_1*gam_2*VSM_coef[p]^(-3))
  d2 <- 4 * gam_1 * gam_2^(2) * VSM_coef[p-1]^(-3) * Func_eta(VSM_coef, p, gam_1, gam_2)^(-2)
  e2 <- 4 * gam_1 * gam_2^2 * VSM_coef[p-1]^(-3) * Func_eta(VSM_coef, p-2, gam_1, gam_2)
  f2 <- -2 * gam_2 * Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2) * (VSM_coef[p-1] - 4 * gam_1 * gam_2 * VSM_coef[p-1]^(-3))
  g2 <- 4 * gam_1 * gam_2^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2)
  h2 <- -2 * gam_2 * Func_eta(VSM_coef, p, gam_1, gam_2)^(-2) * (VSM_coef[p] + 4 * gam_1 * gam_2 * VSM_coef[p]^(-3))

  Value <- 
    {log_likelihood(reg_var, Basis_function, p, p-1)
     -(1/2)*(a2 * VSM_coef[p-1]^(2) + 2 * b2 * VSM_coef[p-1] + c2 * VSM_coef[p]^(2) + 2 * d2 * VSM_coef[p])
     -(gam_1/2)*(- 2 * ( a2 * (lambda(VSM_coef, p-2, gam_1, gam_2) - lambda(VSM_coef, p-1, gam_1, gam_2)) + b2 * (e2-f2) ) 
                - 2 * (g2 * d2 + lambda(VSM_coef, p+1, gam_1, gam_2) * c2)
                + 4 * (d2 * h2 + lambda(VSM_coef, p, gam_1, gam_2) * c2)
                - 2 * (c2 * lambda(VSM_coef, p-1, gam_1, gam_2) + h2 * f2 + d2 * b2 + lambda(VSM_coef, p, gam_1, gam_2) * a2) 
                + 2 * (f2 * b2 + lambda(VSM_coef, p-1, gam_1, gam_2) * a2))
    +(gam_2/2) * (a2*lambda(VSM_coef, p-1, gam_1, gam_2)^(-1) - b2 * lambda(VSM_coef, p-1, gam_1, gam_2)^(-2) * f2
                  + c2 * lambda(VSM_coef, p, gam_1, gam_2)^(-1) 
                  - h2 * lambda(VSM_coef, p, gam_1, gam_2)^(-2) * d2)}
  
  return(Value)
}

####################################
w_p_p_minus_2 <- function(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  a <- 32*gam_1^(2)*gam_2^(3)*VSM_coef[p]^(-3)*VSM_coef[p-2]^(-3)*Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-3)
  b <- 4 * gam_1 * gam_2 ^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2)
  c <- -2 * gam_2 * Func_eta(VSM_coef, p-2, gam_1, gam_2)^(-2) * (VSM_coef[p-2] + 4 * gam_1 * gam_2 * VSM_coef[p-2]^(-3))
  d <- 4 * gam_1 * gam_2^2 * VSM_coef[p-2]^(-3) * Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2)
  e <- c <- -2 * gam_2 * Func_eta(VSM_coef, p, gam_1, gam_2)^(-2) * (VSM_coef[p] + 4 * gam_1 * gam_2 * VSM_coef[p]^(-3))
  
  Value <- {log_likelihood(reg_var, Basis_function, p, p-2)
      -(1/2)*(a * VSM_coef[p-1]^(2))
      -(gam_1/2)*(-2*(a*(lambda(VSM_coef, p-2, gam_1, gam_2)-lambda(VSM_coef, p-1, gam_1, gam_2)) + b*(c-d))
                -2*(e*d + lambda(VSM_coef, p, gam_1, gam_2) * a)+2*(d * b + lambda(VSM_coef, p-1, gam_1, gam_2) * a))
      +(gam_2/2) * (a * lambda(VSM_coef, p-1, gam_1, gam_2)^(-1) - b * lambda(VSM_coef, p-1, gam_1, gam_2)^(-2) * d)}
  
  return(Value)
}

####################################
w_p_p_minus_3 <- function(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  Value <- {log_likelihood(reg_var, Basis_function, p, p-3)
    +gam_1 * (4 * gam_1 * gam_2 ^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2))*( 4 * gam_1 * gam_2^2 * VSM_coef[p-3]^(-3) * Func_eta(VSM_coef, p-2, gam_1, gam_2)^(-2) )}

  return(Value)
}

####################################
w_p_p_plus_1 <- function(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  
  a3 <- -16*gam_1*gam_2^(2)*Func_eta(VSM_coef, p, gam_1, gam_2)^(-3)*VSM_coef[p+1]^(-3)*(VSM_coef[p] - 4*gam_1*gam_2*VSM_coef[p]^(-3))
  b3 <- 4 * gam_1 * gam_2 ^2 * VSM_coef[p+1]^(-3) * Func_eta(VSM_coef, p, gam_1, gam_2)^(-2)
  c3 <- -16*gam_1*gam_2^2*Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-3)*VSM_coef[p]^(-3)*(VSM_coef[p+1] + 4*gam_1*gam_2*VSM_coef[p+1]^(-3))
  d3 <- 4 * gam_1 * gam_2^(2)  * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2)
  e3 <- 4 * gam_1 * gam_2^2 * VSM_coef[p+1]^(-3) * Func_eta(VSM_coef, p+2, gam_1, gam_2)
  f3 <- -2 * gam_2 * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2) * (VSM_coef[p+1] - 4 * gam_1 * gam_2 * VSM_coef[p+1]^(-3))
  g3 <- -2 * gam_2 * Func_eta(VSM_coef, p, gam_1, gam_2)^(-2) * (VSM_coef[p] - 4 * gam_1 * gam_2 * VSM_coef[p]^(-3))
  h3 <- 4 * gam_1 * gam_2^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p-1, gam_1, gam_2)^(-2)
  
  
  Value <- 
    {log_likelihood(reg_var, Basis_function, p, p+1)
      -(1/2)*(a3 * VSM_coef[p]^(2) + 2 * b3 * VSM_coef[p] + c3 * VSM_coef[p+1]^(2) + 2 * d3 * VSM_coef[p])
      -(gam_1/2)*(- 2 * ( c3 * (lambda(VSM_coef, p+2, gam_1, gam_2) - lambda(VSM_coef, p+1, gam_1, gam_2)) + d3 * (e3-f3) ) 
              + 2 * (f3 * d3 + lambda(VSM_coef, p+1, gam_1, gam_2) * c3)
              - 2 * (c3 * lambda(VSM_coef, p, gam_1, gam_2) + d3 * b3 + f3 * g3 + lambda(VSM_coef, p+1, gam_1, gam_2) * a3)
              + 4 * (b3 * g3 + lambda(VSM_coef, p, gam_1, gam_2) * a3)
              - 2 * (a3 * lambda(VSM_coef, p-1, gam_1, gam_2) + h3 * h3)) 
    +(gam_2/2) * (a3*lambda(VSM_coef, p, gam_1, gam_2)^(-1) - g3 * lambda(VSM_coef, p, gam_1, gam_2)^(-2) * b3 
                + c3 * lambda(VSM_coef, p+1, gam_1, gam_2)^(-1) 
                - d3 * lambda(VSM_coef, p+1, gam_1, gam_2)^(-2) * f3)}

  return(Value)
}

####################################
w_p_p_plus_2 <- function(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  a1 <- 32*gam_1^(2)*gam_2^(3)*VSM_coef[p]^(-3)*VSM_coef[p+2]^(-3)*Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-3)
  b1 <- 4 * gam_1 * gam_2 ^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2)
  c1 <- -2 * gam_2 * Func_eta(VSM_coef, p+2, gam_1, gam_2)^(-2) * (VSM_coef[p+2] + 4 * gam_1 * gam_2 * VSM_coef[p+2]^(-3))
  d1 <- 4 * gam_1 * gam_2^2 * VSM_coef[p+2]^(-3) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2)
  e1 <- c <- -2 * gam_2 * Func_eta(VSM_coef, p, gam_1, gam_2)^(-2) * (VSM_coef[p] + 4 * gam_1 * gam_2 * VSM_coef[p]^(-3))
  
  Value <- 
    {log_likelihood(reg_var, Basis_function, p, p+2)
      -(1/2)*(a1 * VSM_coef[p+1]^(2))
      -(gam_1/2)*(-2*(a1*(lambda(VSM_coef, p+2, gam_1, gam_2)-lambda(VSM_coef, p+1, gam_1, gam_2)) + b1*(c1-d1))
      -2*(e1*d1 + lambda(VSM_coef, p, gam_1, gam_2) * a1) + 2*(d1 * b1 + lambda(VSM_coef, p+1, gam_1, gam_2) * a1))
      +(gam_2/2) * (a1 * lambda(VSM_coef, p+1, gam_1, gam_2)^(-1) - b1 * lambda(VSM_coef, p+1, gam_1, gam_2)^(-2) * d1)}
  
  return(Value)
}

####################################
w_p_p_plus_3 <- function(VSM_coef, reg_var, Basis_function, gam_1, gam_2, p){
  Value <- {log_likelihood(reg_var, Basis_function, p, p+3)
   + gam_1 * (4 * gam_1 * gam_2 ^2 * VSM_coef[p]^(-3) * Func_eta(VSM_coef, p+1, gam_1, gam_2)^(-2) )*( 4 * gam_1 * gam_2^2 * VSM_coef[p+3]^(-3) * Func_eta(VSM_coef, p+2, gam_1, gam_2)^(-2) )}

  return(Value)
}
