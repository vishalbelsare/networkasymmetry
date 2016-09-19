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

system.time(sourceDir(paste0(getwd(),"/R")))

# strat: will have to give up on estimating elasticities of substitution. 
# bc: would have to solve for lambda, gamma, every time. doubles computation, really.
# go straight to covariance estimation. 

# 0. guess volatility parameters.
# 1. take given parameters, z, gamma, lambda, draw random shocks to them.
# 2. solve model again
# 3. repeat 1-2, say 25 times.
# 4. calculate: g_{it} = log(s_{it}) - log(s_i), etc, t=1,...,T
# 5. calculate cov(g_{it},g_{jt})
# 6. see how close that is to data. then go to step 0.





