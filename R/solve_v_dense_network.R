########################################################################
# solve_v_dense_network.R
# Solve model, given parameters.
# License: MIT

# ""
# Jesse Tweedle
# Oct 6, 2016
########################################################################

solve_v_dense_network <- function(R,N,args) {

  beta <- args$beta
  z <- args$z
  ir <- args$ir
  # need gamma, lambda, here. single vector, multiplies by p_i0 in p_r, p_m eqs.
  lambda <- args$lambda
  gamma <- args$gamma

  epsilon <- args$epsilon
  eta <- args$eta

  p_i0 <- p_i1 <- rep_len(1,N)

  tol <- 1e-5 / sqrt(N)

  obj = tol + 1

  while (obj > tol) {

    # save last iteration of parameters
    p_i0 <- p_i1

    # calculate new p_r1
    p_r <- sum(lambda*p_i0^(1-epsilon))^(1/(1-epsilon)) # could be epsilon instead.
    p_m <- sum(gamma*p_i0^(1-eta))^(1/(1-eta))

    # calculate new eta ( = unit intermediate cost)
    p_i1 <- (beta^(-beta) * (1-beta)^(beta-1)) * z^(-1) * p_m^(1-beta)

    obj <- sqrt(sum((p_i1-p_i0)^2))
    print(obj)
  }

  obj = tol + 1
  v0 <- v1 <- rep_len(1,N)
  while (obj > tol) {
    v0 <- v1
    v1 <- lambda * p_i1^(1-epsilon) * sum(p_r^(epsilon-1) * (ir %*% (beta*v0))) + gamma * p_i1^(1-eta) * sum(p_m^(eta-1)*(1-beta)*v0)
    v1 <- v1 / sum(v1)
    obj <- (v1-v0)^2 %>% sum() %>% sqrt()
    print(obj)
  }

  return(list(v=v1,p_r=p_r,p_i=p_i1))
}
