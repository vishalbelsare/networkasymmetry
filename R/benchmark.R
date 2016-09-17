#### Benchmarking.


# Make up income by region. Will have to do this and the next one at the same time,
# so that beta, s, regional income, are all consistent.

  
benchmark <- function(R,N,args=initialize_fake_links(R,N)) {

  # Generate fake edges, regional incomes, plant characteristics
  args <- initialize_fake_links(R,N)
  
  # Region-plant demand edges
  Er <- args$Er
  # Plant-plant demand edges
  En <- args$En
  # Region income
  I <- args$I
  
  # Plant characteristics:
  
  # Plant's region
  ir <- args$ir
  # Value-added share
  beta <- args$beta
  # Output
  s <- args$s

  tol <- 1e-5 / N # should be a function of N, because of normalization by sum(y1).

  y0 <- rep_len(0,N) 
  y1 <- rep_len(1,N) 
  obj <- sum((y1-y0)^2) %>% sqrt()
  
  while (obj > tol) {
    y0 <- y1
  
    u <- I/(Er %*% y0) %>% as.matrix()
    v <- s/(En %*% y0) %>% as.matrix()
    
    y1 <- s * (t(Er) %*% u + t(En) %*% v)^(-1)
    y1 <- y1 / sum(y1)
  
    obj <- sum((y1-y0)^2) %>% sqrt()
  }
  
  xr <- (Er %*% y1)^(-1)
  xn <- (1-beta) * (En %*% y1)^(-1)
  
  A <- (outer(X=xr,Y=y1) * Er) %>% as("sparseMatrix")
  G <- (outer(X=xn,Y=y1) * En) %>% as("sparseMatrix") # this one too big. should split into xn * En and y1 * En
  return(list(A=A,G=G))
}
