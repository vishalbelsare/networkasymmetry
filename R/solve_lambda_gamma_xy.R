########################################################################
# solve_lambda_gamma_xy.R
# Function to solve for unobserved \Lambda and \Gamma,given observed A and G.
# License: MIT

# ""
# Jesse Tweedle
# Oct 10, 2016
########################################################################

solve_lambda_gamma_xy <- function(R,N,args) {

  o <- function(x1,x2) {
    (x1-x2)^2 %>% sum() %>% sqrt()
  }
  
  # debug:
  # args <- initialize_fake_links(R,N)
  # xy <- benchmark_xy(R,N,args=args)
  # args$xn = xy$xn
  # args$xr = xy$xr
  # args$y = xy$y
  # args$sigma <- 2
  # Generate trade cost matrix, plant-plant.

  beta <- args$beta
  C <- args$C
  ir <- args$ir
  sigma <- args$sigma
  Ti <- args$Ti
  Tr <- args$Tr
  s <- args$s
  z <- args$z

  # C <- beta^(-beta) * (1-beta)^(beta-1)
  # z = rlnorm(N,0,1)
  
  xn <- args$xn %>% as.vector()
  xr <- args$xr %>% as.vector()
  y <- args$y %>% as.vector()
  
  eta <- sigma
  epsilon <- sigma
  
  tol <- 1e-5 / sqrt(N)

  # region prices
  p_r0 <- p_r1 <- rep_len(1,R) #.sparseDiagonal(n=R,x=1)

  # plant prices
  p_i0 <- p_i1 <- rep_len(1,N) #.sparseDiagonal(n=N,x=1)

  # initial conditions should be xr, xn, y1, adjusted for productivity.
  u0 <- u1 <- rep_len(1,N)
  v0 <- v1 <- rep_len(1,N)
  r0 <- r1 <- rep_len(1,R)
  i0 <- i1 <- rep_len(1,N)
  
  # so calculate pi using p_i0.
  p_mi <- rep_len(1,N)
  
  Er <- Tr
  Er@x <- rep_len(1,length(Er@x))
  En <- Ti
  En@x <- rep_len(1,length(En@x))

  temp.ti <- Ti
  temp.ti@x <- temp.ti@x^(eta-1)
  temp.tr <- Tr
  temp.tr@x <- temp.tr@x^(epsilon-1)
  
  # set objective: Frobenius norms of all matrices
  # (note: the `p`s are diagonal matrices for technical reasons)
  obj = tol + 1

  # while the difference between iterations is greater than tolerance
  while (obj > tol) {
    # r0 <- r1
    # i0 <- i1
    u0 <- u1
    v0 <- v1
    
#    p_r1 <- (r0^(1/(1-epsilon)) * (temp.tr %*% (i0 * p_i0^(1-epsilon)))^(1/(1-epsilon))) %>% as.vector()
    p_mi <- (u0^((1-beta)/(1-eta)) * (temp.ti %*% ((v0 * p_i0^(1-eta)) %>% matrix(nrow=N,ncol=1)))^((1-beta)/(1-eta))) %>% as.vector()
    p_i1 <- (C * p_mi / z) %>% as.vector()
    
    # temp.r <- (temp.tr %*% ((y * p_i1^(epsilon-1)) %>% matrix(nrow=N,ncol=1))) %>% as.vector()
    # temp.ri <- (((xr * p_r1^(1-epsilon)) %>% matrix(nrow=1,ncol=R)) %*% temp.tr) %>% as.vector()
    # 
    # r1 <- (temp.r * p_r1^(1-epsilon) * xr / sum(i0)) %>% as.vector()
    # i1 <- (temp.ri * p_i1^(1-epsilon) * y / sum(r0)) %>% as.vector()
    
    temp.j <- (temp.ti %*% ((y * p_i1^(eta-1)) %>% matrix(nrow=N,ncol=1))) %>% as.vector()
    temp.i <- (t(temp.ti) %*% (((1-beta)^(-1) * xn * p_mi^(1-eta)) %>% matrix(nrow=N,ncol=1))) %>% as.vector()
    
    # print(dim(En))
    # print((v0) %>% as.matrix() %>% dim())
    u1 <- (1-beta)^(-1) * xn * p_mi^(1-eta) * temp.j / (En %*% (v0 %>% matrix(nrow=N,ncol=1)))
    v1 <- p_i1^(eta-1) * y * temp.i / (u0 %>% matrix(nrow=1,ncol=N) %*% En)
#    v1 <- v1/sum(v1)
#    i1 <- i1/sum(i1)
    # might need to normalize something.
    
    # solve for w, normalize p and p?
    obj = o(v1,v0) + o(u1,u0)
    print(obj)
  }

  print(u1 %>% as.vector())
  print(v1 %>% as.vector())
  
  # return the region-plant and plant-plant demand shares matrices
  return(list(lambda=LAM1,gamma=GAM1,p_i=p_i1,p_r=p_r1))
}
