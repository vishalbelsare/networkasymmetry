########################################################################
# solve_lambda_gamma.R
# Function to solve for unobserved \Lambda and \Gamma,given observed A and G.
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

solve_lambda_gamma <- function(R,N,args) {

  # inputs:
  # A, G, Ti, Tr, beta, z.
  # A: matrix of observed region-plant demand shares, RxN
  # beta: vector of plant labour shares, Nx1
  # G: matrix of observed plant-plant demand shares (IO table), NxN
  # N: number of plants, scalar
  # R: number of regions, scalar
  # sigma: elasticity of substitution, scalar
  # Ti: trade costs matrix for plants, NxN
  # Tr: trade costs matrix for region-plants, RxN
  # z: vector of plant productivity, Nx1

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

#### initialize ####

  # region prices
  p_r0 <- .sparseDiagonal(n=R,x=0)
  p_r1 <- .sparseDiagonal(n=R,x=1)

  # plant prices
  p_i0 <- .sparseDiagonal(n=N,x=0)
  p_i1 <- .sparseDiagonal(n=N,x=1)

  A.df <- A %>% s_to_df()
  G.df <- G %>% s_to_df()

  # region-plant demand shares
  LAM0 <- sparseMatrix(i=A.df$i, j=A.df$j, x=0, dims=c(R,N))
  LAM1 <- A
  # LAM1 <- sparseMatrix(i=A.df$i, j=A.df$j, x=1, dims=c(R,N))

  # plant-plant demand shares
  GAM0 <- sparseMatrix(i=G.df$i, j=G.df$j, x=0, dims=c(N,N))
  GAM1 <- G
  # GAM1 <- sparseMatrix(i=G.df$i, j=G.df$j, x=1, dims=c(N,N))

  # set objective: Frobenius norms of all matrices
  # (note: the `p`s are diagonal matrices for technical reasons)
  obj = norm(LAM1-LAM0,"f") + norm(GAM1-GAM0,"f") + norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")

### Solve for \Lambda and \Gamma:

  # while the difference between iterations is greater than tolerance
  while (obj > tol) {
print(obj)
    # save newest iteration of parameters
    p_i0 <- p_i1
    p_r0 <- p_r1
    GAM0 <- GAM1
    LAM0 <- LAM1

    # calculate new p_r1
    # when doing elementwise multiplication with sparse diagonal matrices, the i-row dimension must be pre-multiplied,
    # and the j-column dimension diagonal must be post-multiplied.
    m1 <- Tr %*%
      p_i0 %>%
      s_to_df() %>%
      mutate(x = x^(1-sigma)) %>%
      df_to_s(dims=c(R,N))

    p_r1 <- rowSums(LAM0 * m1) %>%
      as("matrix") %>%
      s_to_df() %>%
      mutate(x=x^(1/(1-sigma))) %>%
      df_to_s(dims=c(R,R)) %>%
      to_sdiag()

    # calculate new eta ( = unit intermediate cost)
    m2 <- Ti %*%
      p_i0 %>%
      s_to_df() %>%
      mutate(x=x^(1-sigma)) %>%
      df_to_s(dims=c(N,N))

    m2.0 <- rowSums(GAM0 * m2) %>%
      as("matrix") %>%
      s_to_df() %>%
      mutate(x=x^(1/(1-sigma)))

    m2.1 <- (m2.0 %>%
               df_to_s(dims=c(N,1))
             )^(1-beta)

    eta <- m2.1 %>%
      to_sdiag() *
      (C %>% to_sdiag())

    # calculate new p_i1
    p_i1 <- eta / z

    # calculate new region-plant demand shares
    LAM1 <- (p_r1 %>%
               s_to_df() %>%
               mutate(x=x^(1-sigma)) %>%
               df_to_s(dims=c(R,R))) %*%
      A *
      (Tr %*%
         p_i1 %>%
         s_to_df() %>%
         mutate(x=x^(sigma-1)) %>%
         df_to_s(dims=c(R,N))
       )

    # calculate new plant-plant demand shares
    m3 <-
      ((1/(1-beta)) %>%
         to_sdiag()) %*%
      (eta %>%
         s_to_df() %>%
         mutate(x=x^(1-sigma)) %>%
         df_to_s(dims=c(N,N))
       )

    r3 <- G *
      (Ti %*%
         p_i1 %>%
         s_to_df() %>%
         mutate(x=x^(sigma-1)) %>%
         df_to_s(dims=c(N,N))
       )

    GAM1 <- m3 %*% r3

    # solve for w, normalize p and p?
    obj = norm(LAM1-LAM0,"f") + norm(GAM1-GAM0,"f") + norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
    if (!is.finite(obj)) {
      break
    }
  }

  # return the region-plant and plant-plant demand shares matrices
  return(list(lambda=LAM1,gamma=GAM1,p_i=p_i1,p_r=p_r1))
}
