########################################################################
# solve_v.R
# Solve model, given parameters.
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

solve_v <- function(R,N,args) {

  # to debug:
  # rm(list=ls(all=TRUE))
  # gc()
  # source("R/initialize_functions.R")
  # source("R/helpers.R")
  # source("R/benchmark.R")
  # source("R/solve_lambda_gamma.R")
  # set.seed(10) # set.seed(10), R=2, N=300 screws up.
  # R <- 5
  # N <- 500
  # argsx <- initialize_fake_links(R,N)
  # print("get fakes")
  # args <- initialize_fakes(R,N,args=argsx)
  #
  # Ax <- args$A
  # Gx <- args$G
  #
  # beta <- args$beta
  # C <- args$C
  # ir <- args$ir
  # s <- args$s
  # Ti <- args$Ti
  # Tr <- args$Tr
  # z <- args$z
  #
  # true_sigma <- 2
  # print("lambda, gamma")
  # lg <- solve_lambda_gamma(R=R,N=N,sigma=true_sigma,args=args)
  #
  # lambda <- lg$lambda
  # gamma <- lg$gamma
  #
  # sigma <- true_sigma

  # now what. solve some shib. want to solve for s, given lambda, gamma, beta, z, Ti, Tr, C, ir.
  # first, price equations?
  # then solve for A, G. and then s is more obvious.

  beta <- args$beta
  C <- args$C
  ir <- args$ir
  s <- args$s
  Ti <- args$Ti
  Tr <- args$Tr
  z <- args$z

  lambda <- args$lambda
  gamma <- args$gamma
  sigma <- args$sigma

# initial guesses: prices from last time we solved it.
  p_r1 <- args$p_r
  p_i1 <- args$p_i

  tol <- 1e-5 / N^2 # might need to make this even smaller might need to calculate optimal tolerance.

  # region prices
  p_r0 <- .sparseDiagonal(n=R,x=0)
#  p_r1 <- .sparseDiagonal(n=R,x=1)

  # plant prices
  p_i0 <- .sparseDiagonal(n=N,x=0)
#  p_i1 <- .sparseDiagonal(n=N,x=1)

  lambda.df <- lambda %>% s_to_df()
  gamma.df <- gamma %>% s_to_df()

  # set objective: Frobenius norms of all matrices
  # (note: the `p`s are diagonal matrices for technical reasons)
  obj = norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")

  # while the difference between iterations is greater than tolerance
  while (obj > tol) {
    # 83, 84.
#    print(p_i1[83:84,83:84])

    # save newest iteration of parameters
    p_i0 <- p_i1
    p_r0 <- p_r1

    # calculate new p_r1
    m1 <- (Tr %*% p_i0) %>%
      s_to_df() %>%
      mutate(x = x^(1-sigma)) %>%
      df_to_s(dims=c(R,N))

    p_r1 <- rowSums(lambda * m1) %>%
      as("matrix") %>%
      s_to_df() %>%
      mutate(x=x^(1/(1-sigma))) %>%
      df_to_s(dims=c(R,R)) %>%
      to_sdiag()

    # calculate new eta ( = unit intermediate cost)
    m2 <- (Ti %*% p_i0) %>%
      s_to_df() %>%
      mutate(x=x^(1-sigma)) %>%
      df_to_s(dims=c(N,N))

    # i=84, j=83. gamma
#    print(rowSums(gamma*m2)[83:84]) # one of these is negative? one of them really is negative. can't be m2?
#    x <- which(gamma<0)
#    print(gamma %>% summary() %>% tbl_df() %>% filter(x<0))

    m2.0 <- rowSums(gamma * m2) %>%
      as("matrix") %>%
      s_to_df() %>%
      mutate(x=x^(1/(1-sigma)))

    m2.1 <- (m2.0 %>%
               df_to_s(dims=c(N,1))
             )^(1-beta)

    eta <- m2.1 %>%
      to_sdiag() *
      (C %>% to_sdiag())

      # 83, 84.
#      print(eta[83:84,83:84])

    # calculate new p_i1
    p_i1 <- eta / z

    # solve for w, normalize p and p?
    # 83, 84.
#    print(p_i1[83:84,83:84])
    obj = norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
    if (!is.finite(obj)) {
      break
    }
  }

  # calculate new observed region-plant demand shares
  A <- (p_r1 %>%
          s_to_df() %>%
          mutate(x=x^(sigma-1)) %>%
          df_to_s(dims=c(R,R))) %*%
    lambda *
    (Tr %*%
       p_i1 %>%
       s_to_df() %>%
       mutate(x=x^(1-sigma)) %>%
       df_to_s(dims=c(R,N))
     )

  # calculate new plant-plant demand shares
  m3 <- ((1-beta) %>% to_sdiag()) %*%
    (eta %>%
       s_to_df() %>%
       mutate(x=x^(sigma-1)) %>%
       df_to_s(dims=c(N,N)))

  r3 <- gamma *
    (Ti %*%
       p_i1 %>%
       s_to_df() %>%
       mutate(x=x^(1-sigma)) %>%
       df_to_s(dims=c(N,N)))

  G <- m3 %*% r3

  # 83, 84.
#  print(colSums(G))

  v0 <- rep_len(0,N)
  v1 <- rep_len(1,N)
  tol <- 1e-5 / N^2

  while (sum(abs(v1-v0)) > tol) {
#    print(sum(abs(v1-v0)))
    v0 <- v1
    I <- ir %*% (t(t(beta))*v0)
    v1 <- t(A) %*% I + t(G) %*% v0
    v1 <- v1 / sum(v1) # normalize weights.
  }

  return(list(v=v1,p_r=p_r1,p_i=p_i1,A=A,G=G))
}
