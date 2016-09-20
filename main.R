########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

library(dplyr)
library(tibble)
library(Matrix)
library(stringr)

rm(list=ls(all=TRUE))
# Sometimes R or RStudio (?) doesn't automatically clean removed objects. I think.
# That's a problem with a NxN matrix when N really big. So clean just in case:
gc()

# Tip from http://jeromyanglim.tumblr.com/post/33418725712/how-to-source-all-r-files-in-a-directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir(paste0(getwd(),"/R"))

#set.seed(10) # set.seed(10), R=2, N=300 screws up.
R <- 3
N <- 50
sigma <- 2 # here.
print("fake edges")
argsx <- initialize_fake_links(R,N)
print("fake shares")
args <- initialize_fakes(R,N,args=c(sigma=sigma,argsx))

print("solve for lambda, gamma; use A and G as initial guesses?")
lg <- solve_lambda_gamma(R,N,args=c(sigma=sigma,args))

print("solve for s, A, G") # should try to use existing s, A, G, as intial arguments.
solved <- solve_v(R,N,args=c(sigma=sigma,args,lg))
# returns v, p_r, p_i, A, G.

# ok, ok. now?
#plot(log(args$s),log(solved$v))

# 0. guess volatility parameters.
# 1. take given parameters, z, gamma, lambda, draw random shocks to them.
# 2. solve model again
# 3. repeat 1-2, say 25 times.
# 4. calculate: g_{it} = log(s_{it}) - log(s_i), etc, t=1,...,T
# 5. calculate cov(g_{it},g_{jt})
# 6. see how close that is to data. then go to step 0.

# shocks:
v <- solved$v
z <- args$z
lambda <- lg$lambda
gamma <- lg$gamma

sigma_z <- 0.02
sigma_l <- 0.02
sigma_g <- 0.02

# lambda_prime <- lambda
# lambda_prime@x <- lambda_prime@x * (1 + rnorm(length(lambda_prime@x),0,sigma_l))
#
# gamma_prime <- gamma
# gamma_prime@x <- gamma_prime@x * (1 + rnorm(length(gamma_prime@x),0,sigma_g))
#
# z_prime <- z * (1 + rnorm(N,0,sigma_z))
#
# # do this:
# #solved_prime <-
# lg_prime <- list(gamma=gamma_prime,lambda=lambda_prime)
# args_prime <- args
# args_prime$z <- z_prime
#
# #initial_p <- list(p_r=solved$p_r,p_i=solved$p_i)
# source("R/solve_v.R")
# solved_prime <- solve_v(R,N,args=c(sigma=sigma,args_prime,lg_prime,solved))
# str(solved_prime)
#
# plot(log(v),log(solved_prime$v))

# beta, C, ir, lambda, gamma, Ti, Tr, v, z
invariant <- list(
  beta=args$beta,
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  sigma=sigma,
  Ti=args$Ti,
  Tr=args$Tr,
  v=solved$v
)

initial <- list(
  lambda=lg$lambda,
  gamma=lg$gamma,
  z=args$z
)

# parameters <- list(
#   sigma_l=0.02,
#   sigma_g=0.02,
#   sigma_z=0.02
# )

# paratemers = c(sigma_l, sigma_g, sigma_z)
parameters <- c(0.02,0.02,0.02)

T <- 100

covariance <- function(parameters,R,N,T,invariant,initial) {
print(parameters)
  sigma_l <- parameters[1]
  sigma_g <- parameters[2]
  sigma_z <- parameters[3]

  g <- matrix(0,nrow=N,ncol=T)
  for (t in 1:T) {
#    print(t)
    # don't use log growth rates.

    lambda_prime <- initial$lambda
    lambda_prime@x <- lambda_prime@x * (1 + rnorm(length(lambda_prime@x),0,sigma_l))

    gamma_prime <- initial$gamma
    gamma_prime@x <- gamma_prime@x * (1 + rnorm(length(gamma_prime@x),0,sigma_g))

#print(gamma_prime[84,83])

    z_prime <- initial$z * (1 + rnorm(N,0,sigma_z))

    current <- list(
      lambda=lambda_prime,
      gamma=gamma_prime,
      z=z_prime
      )

    solved_prime <- solve_v(R,N,args=c(sigma=sigma,invariant,current))

    g_prime = log(solved_prime$v) - log(v)
    g[,t] <- g_prime[,1]
  }

  # say this is the data.
  return(cov(t(g)) %>% as("sparseMatrix") %>% summary() %>% tbl_df())
}

# do T=25 times.
# g <- matrix(0,nrow=N,ncol=T)
# for (t in 1:T) {
#   print(t)
#   lambda_prime <- lambda
#   lambda_prime@x <- lambda_prime@x * (1 + rnorm(length(lambda_prime@x),0,sigma_l))
#
#   gamma_prime <- gamma
#   gamma_prime@x <- gamma_prime@x * (1 + rnorm(length(gamma_prime@x),0,sigma_g))
#
#   z_prime <- z * (1 + rnorm(N,0,sigma_z))
#
#   # do this:
#   lg_prime <- list(gamma=gamma_prime,lambda=lambda_prime)
#   args_prime <- args
#   args_prime$z <- z_prime
#
#   solved_prime <- solve_v(R,N,args=c(sigma=sigma,args_prime,lg_prime,solved))
#
#   g_prime = log(solved_prime$v) - log(v)
#   g[,t] <- g_prime[,1]
# }

# say this is the data.
#X <- cov(t(g)) %>% as("sparseMatrix") %>% summary() %>% tbl_df()
print("calculate `data`")
X <- covariance(parameters,R,N,T,invariant,initial)

# the covariance function can be parallelled.

objective <- function(parameters,R,N,T,invariant,initial,X) {
  Xhat <- covariance(parameters,R,N,T,invariant,initial)
  obj <- (X %>% inner_join(Xhat,by=c("i","j")) %>% mutate(obj=(x.x-x.y)^2) %>% summarize(obj=sum(obj)))[[1]]
  return(obj)
}

# X <- covariance(parameters,R,N,25,invariant,initial)
# # say this is data. certain sigmas, all 0.02 for now.
# # then do SA?
#
# Xhat <- covariance(parameters,R,N,25,invariant,initial)
# # now do this a million times, trying to get sigmas.
# XX <- X %>% left_join(Xhat,by=c("i","j"))
# ok cool cool.
# summary(lm(x.x ~ x.y, data=XX %>% filter(i<=j)))

library(GenSA)
# f <- function(x) {
#   x[1]^2 + x[2]^2 - 4
# }

# now change covar to be a vector of sigmas. l, g, z.
# they have to be first elements.

# ans <- GenSA(c(1,1),f,lower=c(-1,-1),upper=c(1,1))
# str(ans)

# dang.

ans <- GenSA(par=c(0.01,0.01,0.01),fn=objective,lower=c(1e-3,1e-3,1e-3),upper=c(0.2,0.2,0.2),R=R,N=N,T=T,invariant=invariant,initial=initial,X=X,control=list(verbose=TRUE))
str(ans)
