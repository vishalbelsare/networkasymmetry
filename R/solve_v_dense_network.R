########################################################################
# solve_v.R
# Solve model, given parameters.
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

solve_v_dense_network <- function(R,N,args) {

  beta <- args$beta
  C <- args$C
  ir <- args$ir
  s <- args$s
  v <- args$v
  z <- args$z

  sigma <- args$sigma

  # initial guesses: prices from last time we solved it.
  if (is.null(args$p_r)) {
    p_r0 <- p_r1 <- rep_len(1,R)
  } else {
    p_r0 <- p_r1 <- rowSums(args$p_r)
  }
  if (is.null(args$p_i)) {
    p_i0 <- p_i1 <- rep_len(1,N)
  } else {
    p_i0 <- p_i1 <- rowSums(args$p_i)
  }

  tol <- 1e-5 #/ N # might need to make this even smaller might need to calculate optimal tolerance.

  # set objective: Frobenius norms of all matrices
  # (note: the `p`s are diagonal matrices for technical reasons)
  obj = tol + 1 #norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")

  while (obj > tol) {

    # save last iteration of parameters
    p_i0 <- p_i1
    p_r0 <- p_r1

    # calculate new p_r1
    #rowSums(lambda * temp)^(1/(1-sigma)) %>% to_sdiag()
    p_r1 <- sum(p_i0^(1-sigma))^(1/(1-sigma)) # could be epsilon instead.

    p_m <- sum(p_i0^(1-sigma))^(1/(1-sigma))

    # calculate new eta ( = unit intermediate cost)
#    p_i1 <- z^(-1) * (1-beta)^(beta-1) * beta^(-beta) * (sum(p_i0^(1-sigma))^((1-beta)/(1-sigma)))
    p_i1 <- z^(-1) * (1-beta)^(beta-1) * beta^(-beta) * p_m^(1-beta)

    # calculate new p_i1
#    p_i1 <- eta / z
#    print(rowSums(eta))

    obj = sum((p_r1-p_r0)^2) + sum((p_i1-p_i0)^2)
    if (!is.finite(obj)) {
      break
    }
  }

  # calculate new observed region-plant demand shares. nah. may need to do some algebra.
  # A is now an R length vector
  A <- p_i1^(1-sigma) / p_r1^(1-sigma) # I think it's the same but without beta
#  A <- temp.1 %*% (lambda * temp.2)

  # G is now an N length vector.
#  G <- temp.3 %*% (gamma * temp.4)
  G <- (1-beta) * p_i1^(1-sigma) / p_m^(1-sigma) #?

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
#    v1 <- t(A) %*% I + t(G) %*% v0
    v1 <- A * sum(I) + G * sum(v0)
    v1 <- v1 / sum(v1) # normalize weights.
    obj <- sum(abs(v1-v0))
  }

  return(list(v=v1,p_r=p_r1,p_i=p_i1,A=A,G=G))
}
