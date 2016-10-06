########################################################################
# solve_v_dense_network.R
# Solve model, given parameters.
# License: MIT

# ""
# Jesse Tweedle
# Oct 6, 2016
########################################################################

solve_v_dense_network <- function(R,N,args) {

  C <- args$C
  beta <- args$beta
  z <- args$z

  epsilon <- args$epsilon
  eta <- args$eta

  p_i0 <- p_i1 <- rep_len(1,N)

  tol <- 1e-5 / N

  obj = tol + 1

  while (obj > tol) {

    # save last iteration of parameters
    p_i0 <- p_i1

    # calculate new p_r1
    p_r <- sum((1/N)*p_i0^(1-epsilon))^(1/(1-epsilon)) # could be epsilon instead.
    p_m <- sum((1/N)*p_i0^(1-eta))^(1/(1-eta))

    # calculate new eta ( = unit intermediate cost)
    p_i1 <- C * z^(-1) * p_m^(1-beta)

    obj = sum((p_i1-p_i0)^2)
  }

  A <- p_i1^(1-epsilon) / p_r^(1-epsilon)
  G <- p_i1^(1-eta) / p_m^(1-eta) # But needs to be multiplied by (1-beta), really.

  v0 <- v1 <- rep_len(1/N,N)
  while ((v1-v0)^2 %>% sum() %>% sqrt()) {
    v0 <- v1
    v1 <- A * sum(beta*v0) + G * sum((1-beta)*v0)
  }

  return(list(v=v1,p_r=p_r,p_i=p_i1))
}
