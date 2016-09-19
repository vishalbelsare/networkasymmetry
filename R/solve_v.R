########################################################################
# solve_v.R
# Solve model, given parameters.
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

#solve_v <- function(R,N,args) { 

# to debug:
rm(list=ls(all=TRUE))
gc()
source("R/initialize_functions.R")
source("R/helpers.R")
source("R/benchmark.R")
source("R/solve_lambda_gamma.R")
R <- 2
N <- 30
argsx <- initialize_fake_links(R,N)
args <- initialize_fakes(R,N,args=argsx)
# could try "attach" instead.
#initialize_fakes <- function(R,N,args) {
  #args <- benchmark(R,N,args=argsx)
  
Ax <- args$A
Gx <- args$G

beta <- args$beta
C <- args$C
ir <- args$ir
s <- args$s
Ti <- args$Ti
Tr <- args$Tr
z <- args$z

true_sigma <- 2

#print(system.time(
  lg <- solve_lambda_gamma(R=R,N=N,sigma=true_sigma,args=args)
  #))

lambda <- lg$lambda
gamma <- lg$gamma

sigma <- 2

# now what. solve some shib. want to solve for s, given lambda, gamma, beta, z, Ti, Tr, C, ir.
# first, price equations?
# then solve for A, G. and then s is more obvious.

#
tol <- 1e-5 / N^2
# region prices
p_r0 <- .sparseDiagonal(n=R,x=0)
p_r1 <- .sparseDiagonal(n=R,x=1)

# plant prices
p_i0 <- .sparseDiagonal(n=N,x=0)
p_i1 <- .sparseDiagonal(n=N,x=1)

lambda.df <- lambda %>% s_to_df()
gamma.df <- gamma %>% s_to_df()

# set objective: Frobenius norms of all matrices 
# (note: the `p`s are diagonal matrices for technical reasons)
obj = norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")

### Solve for \Lambda and \Gamma:   

# might depend on w.

# while the difference between iterations is greater than tolerance
while (obj > tol) {
#  print(obj)
  # save newest iteration of parameters
  p_i0 <- p_i1
  p_r0 <- p_r1

  # calculate new p_r1
  m1 <- Tr %*% 
    p_i0 %>% 
    s_to_df() %>% 
    mutate(x = x^(1-sigma)) %>% 
    df_to_s(dims=c(R,N))
  
  p_r1 <- rowSums(lambda * m1) %>% 
    as("matrix") %>% 
    s_to_df() %>% 
    mutate(x=x^(1/(1-sigma))) %>% 
    df_to_s(dims=c(R,R)) %>% 
    to_sdiag()
  
  # calculate new eta ( = unit intermediate cost)
  m2 <- (Ti %*% p_i0) %>% 
    s_to_df() %>% 
    mutate(x=x^(1-sigma)) %>% 
    df_to_s(dims=c(N,N))
  
  m2.0 <- rowSums(gamma * m2) %>% 
    as("matrix") %>% 
    s_to_df() %>%
    mutate(x=x^(1/(1-sigma)))
  
  m2.1 <- (m2.0 %>% 
             df_to_s(dims=c(N,1))
           )^(1-beta)
  
  eta <- m2.1 %>%
    to_sdiag() * 
    (C %>% to_sdiag())
  
  # calculate new p_i1    
  p_i1 <- eta / z
  
  # solve for w, normalize p and p?
  obj = norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
}

# calculate new observed region-plant demand shares
#  LAM1 <- (p_r1 %>% s_to_df() %>% mutate(x=x^(1-sigma)) %>% df_to_s(dims=c(R,R))) %*% A * (Tr %*% p_i1 %>% s_to_df() %>% mutate(x=x^(sigma-1)) %>% df_to_s(dims=c(R,N)))
A <- (p_r1 %>% 
        s_to_df() %>%
        mutate(x=x^(sigma-1)) %>%
        df_to_s(dims=c(R,R))) %*% 
  lambda * 
  (Tr %*% 
     p_i1 %>% 
     s_to_df() %>%
     mutate(x=x^(1-sigma)) %>%
     df_to_s(dims=c(R,N))
   )

# calculate new plant-plant demand shares
m3 <- ((1-beta) %>% to_sdiag()) %*% 
  (eta %>% 
     s_to_df() %>% 
     mutate(x=x^(sigma-1)) %>% 
     df_to_s(dims=c(N,N)))

r3 <- gamma * 
  (Ti %*% 
     p_i1 %>% 
     s_to_df() %>% 
     mutate(x=x^(1-sigma)) %>% 
     df_to_s(dims=c(N,N)))

G <- m3 %*% r3

# productivity isn't enough to justify output, only other way is aggregate gamma?

w0 <- 0
w1 <- 1
s0 <- rep_len(0,N)
s1 <- rep_len(1,N)
tol <- 1e-5 / N^2

# while ((w1-w0)^2 + sum((s1-s0)^2) > tol) {
#   print(sum((s1-s0)^2 + (w1-w0)^2))
while (sum(abs(s1-s0)) + (w1-w0)^2 > tol) {
  print(sum(abs(s1-s0))+ (w1-w0)^2)
  w0 <- w1
  s0 <- s1
  I <- ir %*% (t(t(beta))*s0)
#    rowSums(t(t(ir) * t(t(beta))*s0)) #sum(beta * s0)
  s1 <- w0 * t(A) %*% I + t(G) %*% s0
  w1 <- mean(s1/s0)
  s1 <- s1 / sum(s1)
}

plot(log(s1),log(s))

# A <- args$A
# G <- args$G


# the