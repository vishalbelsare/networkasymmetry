########################################################################
# solve_lambda_gamma.R
# Function to solve for unobserved \Lambda and \Gamma,given observed A and G.
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

solve_lambda_gamma <- function(R,N,args) {

  A <- args$A
  beta <- args$beta
  C <- args$C
  G <- args$G
  ir <- args$ir
  sigma <- args$sigma
  Ti <- args$Ti
  Tr <- args$Tr
  s <- args$s
  z <- args$z

  tol <- 1e-5

  # region prices
  p_r0 <- p_r1 <- .sparseDiagonal(n=R,x=1)

  # plant prices
  p_i0 <- p_i1 <- .sparseDiagonal(n=N,x=1)

  LAM0 <- LAM1 <- A
  GAM0 <- GAM1 <- G

  # set objective: Frobenius norms of all matrices
  # (note: the `p`s are diagonal matrices for technical reasons)
  obj = tol + 1

  # while the difference between iterations is greater than tolerance
  while (obj > tol) {

    # save last iteration of parameters
    p_i0 <- p_i1
    p_r0 <- p_r1
    GAM0 <- GAM1
    LAM0 <- LAM1

    # calculate new p_r1
    # when doing elementwise multiplication with sparse diagonal matrices, the i-row dimension must be pre-multiplied,
    # and the j-column dimension diagonal must be post-multiplied.
    m1 <- Tr %*% p_i0
    m1@x <- m1@x^(1-sigma)

    p_r1 <- rowSums(LAM0 * m1)^(1/(1-sigma)) %>% to_sdiag()

    # calculate new eta ( = unit intermediate cost)
    m2 <- Ti %*% p_i0
    m2@x <- m2@x^(1-sigma)
    m3 <- rowSums(GAM0 * m2)^((1-beta)/(1-sigma))

    eta <- (m3 * C) %>% to_sdiag()

    # calculate new p_i1
    p_i1 <- eta / z

    # calculate new region-plant demand shares
    temp.1 <- p_r1
    temp.1@x <- temp.1@x^(1-sigma)

    temp.2 <- Tr %*% p_i1
    temp.2@x <- temp.2@x^(sigma-1)

    LAM1 <- temp.1 %*% (A * temp.2)

    # calculate new plant-plant demand shares
    temp.eta <- eta
    temp.eta@x <- temp.eta@x^(1-sigma)
    temp.3 <- ((1/(1-beta)) %>% to_sdiag()) %*% temp.eta

    temp.4 <- Ti %*% p_i1
    temp.4@x <- temp.4@x^(sigma-1)
    GAM1 <- temp.3 %*% (G * temp.4)

    # solve for w, normalize p and p?
    obj = norm(LAM1-LAM0,"f") + norm(GAM1-GAM0,"f") + norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
    if (!is.finite(obj)) {
      break
    }
  }

  # return the region-plant and plant-plant demand shares matrices
  return(list(lambda=LAM1,gamma=GAM1,p_i=p_i1,p_r=p_r1))
}
