########################
# main.R
# Clean and load relevant files
####
# For Granularity, Network Asymmetry and Aggregate Volatility
# Jesse Tweedle
# Sep 15, 2016
########################

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

# Load all files in ./R subdirectory
# helpers.R contains some help functions.
# the others solve things.
sourceDir(paste0(getwd(),"/R"))

# how it works. start with R, N, need to generate Er, En, beta, z, s, ir.
R = 2
N = 5
# initialize fake data; generate Er, En, beta, z, s, ir
# args <- initialize_fakes(R,N)
# args <- solve_lambda_gamma(R=R,N=N)
# Otherwise, I have data:
# args <- input()
# should give beta, z, s, Er, En, ir (plant characteristic), R=75, N=30,000-ish
# then call benchmark. then call solve_lambda_gamma.

# could use better organization here. initialize, save the data. pass it in, etc.
# need a flat heirarchy.