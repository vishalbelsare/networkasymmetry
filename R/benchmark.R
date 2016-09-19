########################################################################
# benchmark.R
# Function to benchmark A and G consistent with given parameters
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

benchmark <- function(R,N,args) { #=initialize_fake_links(R,N)) {

# to debug:
  # rm(list=ls(all=TRUE))
  # gc()
  # source("R/initialize_functions.R")
  # source("R/helpers.R")
  # R <- 50
  # N <- 500

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
  y1 <- rep_len(1,N) 
  obj <- norm(A1-A0,"f") + norm(G1-G0,"f") # sum((y1-y0)^2) %>% sqrt()
  
  # While y1 is too far from y0:
  while (obj > tol) {
    if (obj < 1.03e-5) {
#      print(sum(rowSums(G1)+beta))
    }
    A0<-A1
    G0<-G1
    # Save old y1 in y0
    y0 <- y1
  
    # Intermediate calculations
    u <- I/(Er %*% y0) %>% as.matrix()
    v <- s/(En %*% y0) %>% as.matrix()
  
    # Update y1  
    y1 <- s * (t(Er) %*% u + t(En) %*% v)^(-1)
    
    # Normalize y1; otherwise y1 and y0 both approach 0 and satisfies obj < tol automatically.
#    y1 <- y1 / mean(y1)

    xr <- (Er %*% y1)^(-1)
    xn <- (1-beta) * (En %*% y1)^(-1)

    # xr <- xr * sum(y1)
    # xn <- xn * sum(y1)
    
    # Calculate the other bit of the demand matrices. # why are some 1e-239??
    # xr <- (Er %*% y1)^(-1)
    # xn <- (1-beta) * (En %*% y1)^(-1)
    
    # Outer product gives A and G. a_{ri} = x^R_r * y_i, g_{ij} = x^N_i * y_j  
    A1 <- (outer(X=xr,Y=y1) * Er) %>% as("sparseMatrix")
    G1 <- (outer(X=xn,Y=y1) * En) %>% as("sparseMatrix") # this one might be too big. should split into xn * En and y1 * En
    
    A1[A1<1e-25] <- 0
    G1[G1<1e-25] <- 0
    # A1[A1>1] <- 1
    # G1[G1>1] <- 1
    # A1[!is.finite(A1)] <- 1
    # G1[!is.finite(G1)] <- 1
    # Recalculate objective.
    obj <- norm(A1-A0,"f") + norm(G1-G0,"f") # sum((y1-y0)^2) %>% sqrt()

  }
  
  # rowSums(A1)
  # rowSums(G1) + beta
  # xxx
  
  # so benchmarking sets some to 0. why?
  # should try to normalize xr, xn instead.
  # Calculate the other bit of the demand matrices. # why are some 1e-239??
  # xr <- (Er %*% y1)^(-1)
  # xn <- (1-beta) * (En %*% y1)^(-1)

  # Outer product gives A and G. a_{ri} = x^R_r * y_i, g_{ij} = x^N_i * y_j  
  # A <- (outer(X=xr,Y=y1) * Er) %>% as("sparseMatrix")
  # G <- (outer(X=xn,Y=y1) * En) %>% as("sparseMatrix") # this one might be too big. should split into xn * En and y1 * En
  return(list(A=A1,G=G1))
}
