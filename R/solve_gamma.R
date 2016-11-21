########################################################################
# solve_lambda_gamma.R
# Function to solve for unobserved \Lambda and \Gamma, given observed A and G.
# License: MIT

# ""
# Jesse Tweedle
# , 2016
########################################################################

solve_gamma <- function(R,N,args) {

  beta <- args$beta
  C <- args$C
  A <- args$A
  G <- args$G
  ir <- args$ir
  eta <- args$eta
  epsilon <- args$epsilon
  Ti <- args$Ti
  Tr <- args$Tr
  s <- args$s
  z <- args$z

  tol <- 1e-5

  # plant prices
  p_i0 <- p_i1 <- .sparseDiagonal(n=N,x=1)

  GAM0 <- GAM1 <- G 
  
  obj = tol + 1
  obj_0 <- obj + 1
  counter <- 0
  
  # while the difference between iterations is greater than tolerance
  while (obj > tol) {

    if (obj > obj_0 | (log(obj_0) - log(obj)) < 0.005) {
      counter <- counter+1
      if (counter>3) {
        break
      }
    } else {
      counter <- 0
    }
    obj_0 <- obj
    
    # save last iteration of parameters
    p_i0 <- p_i1
    GAM0 <- GAM1

    # calculate new p_mi ( = unit intermediate cost)
    m2 <- Ti %*% p_i0
    m2@x <- m2@x^(1-eta)
  
    mxx <- rowSums(GAM0 * m2)
    p_mi <- mxx^(1/(1-eta)) #^((1-beta)/(1-sigma))

    # calculate new p_i1
    p_i1 <- (C * p_mi^(1-beta) / z) %>% to_sdiag()

    # really need p_mi^(eta-1), but that gives mxx back, so use that directly
    temp.3 <- (mxx / (1-beta)) %>% to_sdiag() 
    
    temp.4 <- Ti %*% p_i1
    temp.4@x <- temp.4@x^(eta-1)
    
    # Calculate new estimate of gamma.
    GAM1 <- temp.3 %*% (G * temp.4)

    # solve for w, normalize p and p?
    obj <- sqrt(sum((diag(p_i1) - diag(p_i0))^2))
    print(obj)
  }

  p_r0 <- p_r1 <- .sparseDiagonal(n=R,x=1)
  LAM0 <- LAM1 <- A
  obj = tol + 1
  obj_0 <- obj + 1
  counter <- 0
  
  # while the difference between iterations is greater than tolerance
  while (obj > tol) {
    
    if (obj > obj_0 | (log(obj_0) - log(obj)) < 0.005) {
      counter <- counter+1
      if (counter>3) {
        break
      }
    } else {
      counter <- 0
    }
    obj_0 <- obj
    
    # save last iteration of parameters
    p_r0 <- p_r1
    LAM0 <- LAM1
    
    m1 <- Tr %*% p_i0
    m1@x <- m1@x^(1-epsilon)
    
    p_r1 <- rowSums(LAM0 * m1)^(1/(1-epsilon)) %>% to_sdiag()
    
    temp.1 <- p_r1
    temp.1@x <- temp.1@x^(1-epsilon)
    
    temp.2 <- Tr %*% p_i1
    temp.2@x <- temp.2@x^(epsilon-1)
    
    LAM1 <- temp.1 %*% (A * temp.2)
    obj <- sqrt(sum((diag(p_r1) - diag(p_r0))^2))
    print(obj)
    
  }
  
  # return the region-plant and plant-plant demand shares matrices
  return(list(lambda=LAM1,gamma=GAM1,p_r=p_r1,p_i=p_i1))
}
