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

set.seed(10) # set.seed(10), R=2, N=300 screws up.
R <- 5
N <- 500
sigma <- 2
print("fake edges")
argsx <- initialize_fake_links(R,N)
print("fake shares")
args <- initialize_fakes(R,N,args=c(sigma=sigma,argsx))

print("solve for lambda, gamma; use A and G as initial guesses?")
print(system.time(lg <- solve_lambda_gamma(R,N,args=c(sigma=sigma,args))))
xxx

print("solve for s, A, G") # should try to use existing s, A, G, as intial arguments.
solved <- solve_v(R,N,args=c(sigma=sigma,args,lg))
# returns v, p_r, p_i, A, G.

# ok, ok. now?
plot(log(args$s),log(solved$v))

# 0. guess volatility parameters.
# 1. take given parameters, z, gamma, lambda, draw random shocks to them.
# 2. solve model again
# 3. repeat 1-2, say 25 times.
# 4. calculate: g_{it} = log(s_{it}) - log(s_i), etc, t=1,...,T
# 5. calculate cov(g_{it},g_{jt})
# 6. see how close that is to data. then go to step 0.
