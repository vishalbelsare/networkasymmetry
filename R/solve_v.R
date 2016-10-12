########################################################################
# solve_v.R
# Solve model, given parameters.
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

solve_v <- function(R,N,args) {

  beta <- args$beta
  C <- args$C
  ir <- args$ir
  s <- args$s
  Ti <- args$Ti
  Tr <- args$Tr
  v <- args$v
  z <- args$z

  lambda <- args$lambda
  gamma <- args$gamma
  sigma <- args$sigma

  # initial guesses: prices from last time we solved it.
  p_r0 <- p_r1 <- args$p_r
  p_i0 <- p_i1 <- args$p_i

  tol <- 1e-5 / N # might need to make this even smaller might need to calculate optimal tolerance.

  # set objective: Frobenius norms of all matrices
  # (note: the `p`s are diagonal matrices for technical reasons)
  obj = tol + 1 #norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")

  while (obj > tol) {

    # save last iteration of parameters
    p_i0 <- p_i1
    p_r0 <- p_r1

    # calculate new p_r1
    temp <- (Tr %*% p_i0)
    temp@x <- temp@x^(1-sigma)

    p_r1 <- rowSums(lambda * temp)^(1/(1-sigma)) %>% to_sdiag()

    # calculate new eta ( = unit intermediate cost)
    temp <- (Ti %*% p_i0)
    temp@x <- temp@x^(1-sigma)
    m2 <- rowSums(gamma * temp)^((1-beta)/(1-sigma))
#    eta <- (m2 * C) %>% to_sdiag()
    eta <- m2

    # calculate new p_i1
    p_i1 <- (C * eta / z) %>% to_sdiag()
#    print(rowSums(eta))

    obj = norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
    if (!is.finite(obj)) {
      break
    }
  }

  # calculate new observed region-plant demand shares

  temp.1 <- p_r1
  temp.1@x <- temp.1@x^(sigma-1)

  temp.2 <- Tr %*% p_i1
  temp.2@x <- temp.2@x^(1-sigma)

  A <- temp.1 %*% (lambda * temp.2)

  temp.3 <- ((1-beta) %>% to_sdiag()) %*% (eta %>% to_sdiag())^(sigma-1)
  temp.4 <- Ti %*% p_i1
  temp.4@x <- temp.4@x^(1-sigma)
  G <- temp.3 %*% (gamma * temp.4)


# maybe this should have it's own initial level. yes, def.
  if (is.null(v)) {
    v0 <- v1 <- rep_len(1,N)
  } else {
    v0 <- v1 <- v #rep_len(0,N)
  }

  tol <- 1e-5 / N
  obj <- tol + 1

  while (obj > tol) {
    v0 <- v1
    I <- ir %*% (t(t(beta))*v0)
    v1 <- t(A) %*% I + t(G) %*% v0
    v1 <- v1 / sum(v1) # normalize weights.
    obj <- sum(abs(v1-v0))
  }

  return(list(v=v1,p_r=p_r1,p_i=p_i1,A=A,G=G))
}
