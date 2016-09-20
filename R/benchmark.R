########################################################################
# benchmark.R
# Function to benchmark A and G consistent with given parameters
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

benchmark <- function(R,N,args) {

# to debug:
  # rm(list=ls(all=TRUE))
  # gc()
  # library(dplyr)
  # library(tibble)
  # library(Matrix)
  # library(stringr)
  # source("R/initialize_functions.R")
  # source("R/helpers.R")
  # R <- 2
  # N <- 50
  # args <- initialize_fake_links(R,N)

  # conclusion: benchmarking is hard. no guarantee it will work, really.

  if (missing(args)) {
    args <- initialize_fake_links(R,N)
  }

  # Region-plant demand edges
  Er <- args$Er

  # Plant-plant demand edges
  En <- args$En

  # Region income
  I <- args$I

  # Plant's region
  ir <- args$ir

  # Value-added share
  beta <- args$beta

  # Output
  s <- args$s

  tol <- 1e-5 / N # should be a function of N, because of normalization by sum(y1).

  A0 <- Er * 0
  G0 <- En * 0
  A1 <- Er
  G1 <- En

  # initialize y0, y1, objective.
  y0 <- rep_len(0,N)
  y1 <- rep_len(1/N,N)
  obj <- norm(A1-A0,"f") + norm(G1-G0,"f") # sum((y1-y0)^2) %>% sqrt()

  # While y1 is too far from y0:
  while (obj > tol) {
#    print(str_c("obj: ",obj))

    A0 <- A1
    G0 <- G1

    # Save old y1 in y0
    y0 <- y1

    # Intermediate calculations
    u <- I/(Er %*% y0) %>% as.matrix()
    v <- s/(En %*% y0) %>% as.matrix()

    # Update y1
    y1 <- s * (t(Er) %*% u + t(En) %*% v)^(-1)

    xr <- (Er %*% y1)^(-1)
    xn <- (1-beta) * (En %*% y1)^(-1)
    xn[!is.finite(xn)] <- 0 # try to fix xn->infty problems.

    # Outer product gives A and G. a_{ri} = x^R_r * y_i, g_{ij} = x^N_i * y_j
    A1 <- (outer(X=xr,Y=y1) * Er) %>% as("sparseMatrix")
    G1 <- (outer(X=xn,Y=y1) * En) %>% as("sparseMatrix")
    # this one might be too big. should split into xn * En and y1 * En

# 321, 347. xn has only one y1 with En turned on, and that guy's small?
# just give this up for now. work on actual details later.
# reset G1 to work.

    # A1[A1<1e-25] <- 0
    # G1[G1<1e-25] <- 0

    obj <- norm(A1-A0,"f") + norm(G1-G0,"f")
    if (!is.finite(obj)) {
      break
    }
  }

  b <- rowSums(A1)
  c <- colSums(A1)
  bx <- b==0 | is.na(b)
  cx <- c==0 | is.na(c)

  A1[bx, cx] <- matrix(runif(length(bx[bx])*length(cx[cx]),0,1),nrow=length(bx[bx]),ncol=length(cx[cx]))

  A1 <- A1 / rowSums(A1) # re-normalize.

  b <- rowSums(G1)
  c <- colSums(G1)
  bx <- b==0 | is.na(b)
  cx <- c==0 | is.na(c)

  G1[bx, cx] <- matrix(runif(length(bx[bx])*length(cx[cx]),0,1),nrow=length(bx[bx]),ncol=length(cx[cx]))

  G1 <- (1-beta) * G1 / rowSums(G1) # re-normalize.

  return(list(A=A1,G=G1))
}
