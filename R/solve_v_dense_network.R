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

    obj = sum((p_i1-p_i0)^2)
  }

  A <- lambda * p_i1^(1-epsilon) / p_r^(1-epsilon)
  G <- gamma * p_i1^(1-eta) / p_m^(1-eta) # But needs to be multiplied by (1-beta), really.

  obj = tol + 1
  v0 <- v1 <- rep_len(1/N,N)
  while (obj > tol) {
    print(obj)
    v0 <- v1
    v1 <- A * sum(beta*v0) + G * sum((1-beta)*v0)
    obj <- (v1-v0)^2 %>% sum() %>% sqrt()
  }

  return(list(v=v1,p_r=p_r,p_i=p_i1))
}
