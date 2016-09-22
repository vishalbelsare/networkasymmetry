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
library(tidyr)
library(ggplot2)

# A, G, Ti, Tr, beta, z.
# A: matrix of observed region-plant demand shares, RxN
# beta: vector of plant labour shares, Nx1
# G: matrix of observed plant-plant demand shares (IO table), NxN
# N: number of plants, scalar
# R: number of regions, scalar
# sigma: elasticity of substitution, scalar
# Ti: trade costs matrix for plants, NxN
# Tr: trade costs matrix for region-plants, RxN
# z: vector of plant productivity, Nx1


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

set.seed(10) # set.seed(10), R=2, N=300 screws up.
R <- 73
N <- 10000
sigma <- 2 # here.

print("fake edges")
argsx <- initialize_fake_links(R,N)

print("fake shares")
args <- initialize_fakes(R,N,args=c(sigma=sigma,argsx))

print("solve for lambda, gamma; use A and G as initial guesses?")
lg <- solve_lambda_gamma(R,N,args=c(sigma=sigma,args))

print("solve for s, A, G") # should try to use existing s, A, G, as intial arguments.
solved <- solve_v(R,N,args=c(sigma=sigma,args,lg))

dat <- tibble(v=solved$v[,1],s=args$s)
dat <- dat %>% rownames_to_column() %>% mutate(ss=sum(s),s=s/ss) %>% select(-ss) %>% gather(type,size,v:s)

#plot(log(solved$v),log(args$s))
ggplot(dat %>% spread(type,size),aes(y=s,x=v)) + geom_point() + scale_x_log10() + scale_y_log10()
# doesn't work at all.
ggplot(dat , aes(size,colour=type)) + geom_density() + scale_x_log10()

# returns v, p_r, p_i, A, G.

# shocks:
v <- solved$v
z <- args$z
lambda <- lg$lambda
gamma <- lg$gamma

# sigma_z <- 0.02
# sigma_l <- 0.02
# sigma_g <- 0.02

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

# parameters <- c(0.02,0.02,0.02)
parameters <- c(0.02)

T <- 25

# the covariance function can be parallelled.
covariance <- function(parameters,R,N,T,invariant,initial,randoms) {
  print(str_c("sigma_z: ",parameters))
#  sigma_l <- parameters[1]
  # sigma_g <- parameters[2]
  sigma_z <- parameters[1]

  g <- matrix(0,nrow=N,ncol=T)
  for (t in 1:T) {
#    print(t)

    lambda_prime <- initial$lambda
    # lambda_prime@x <- lambda_prime@x * (1 + qnorm(randoms$lamda_runif,0,sigma_l))

    gamma_prime <- initial$gamma
    # gamma_prime@x <- gamma_prime@x * (1 + qnorm(randoms$gamma_runif,0,sigma_g))

    z_prime <- initial$z * (1 + qnorm(randoms$z_runif[t,] %>% as.vector(),0,sigma_z))

    current <- list(
      lambda=lambda_prime,
      gamma=gamma_prime,
      z=z_prime
    )

    solved_prime <- solve_v(R,N,args=c(invariant,current))

    g_prime = log(solved_prime$v) - log(invariant$v)
    g[,t] <- g_prime[,1]
  }

  # say this is the data.
  return(cov(t(g)) %>% as("sparseMatrix") %>% summary() %>% tbl_df())
}

# say this is the data.

print("calculate ```data```")
# frig, this has to be N^2xT, N^2xT and NxT instead.
randoms_data <- list(
  # lambda_runif=runif(length(lambda@x),0,1),
  # gamma_runif=runif(length(gamma@x),0,1),
  z_runif=matrix(runif(N*T,0,1),nrow=T,ncol=N)
)

X <- covariance(parameters,R,N,T,invariant,initial,randoms_data)

objective <- function(parameters, R, N, T, invariant, initial, X, randoms) {
  # last thing I have to do is change moments.
  # maybe aggregate volatility?

  # invariant$v.
#  invariant$v %>% df_to_s(dims=c(N,N))
  agg_vol <- sum((invariant$v %>% to_sdiag()) %*% (X %>% df_to_s(dims=c(N,N))) %*% (invariant$v %>% to_sdiag()))

  Xhat <- covariance(parameters, R, N, T, invariant, initial, randoms)
  agg_vol_prime <- sum((invariant$v %>% to_sdiag()) %*% (Xhat %>% df_to_s(dims=c(N,N))) %*% (invariant$v %>% to_sdiag()))

  obj <- abs(agg_vol - agg_vol_prime)
  #obj <- (X %>% inner_join(Xhat,by=c("i","j")) %>% mutate(obj=(x.x-x.y)^2) %>% summarize(obj=sum(obj)))[[1]]
  print(str_c("obj: ",obj," vol: ",agg_vol," vol': ",agg_vol_prime))
  return(obj)
}

library(GenSA)
# f <- function(x) {
#   x[1]^2 + x[2]^2 - 4
# }
# ans <- GenSA(c(1,1),f,lower=c(-1,-1),upper=c(1,1))
# str(ans)

# randoms for the simulations.
randoms <- list(
  # lambda_runif=runif(length(lambda@x),0,1),
  # gamma_runif=runif(length(gamma@x),0,1),
  z_runif=matrix(runif(N*T,0,1),nrow=T,ncol=N)
)

ans <- GenSA(par=c(0.01), fn=objective, lower=c(1e-3), upper=c(0.2), R=R, N=N, T=T, invariant=invariant, initial=initial, X=X, randoms=randoms)
# ans <- GenSA(par=c(0.01,0.01,0.01), fn=objective, lower=c(1e-3,1e-3,1e-3), upper=c(0.2,0.2,0.2), R=R, N=N, T=T, invariant=invariant, initial=initial, X=X)
str(ans)
