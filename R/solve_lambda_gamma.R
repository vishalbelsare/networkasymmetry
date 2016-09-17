solve_lambda_gamma <- function(R,N,args=initialize_fakes(R,N)) { #A,beta,G,N,R,Ti,Tr,z) {

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
  Ti <- args$Ti
  Tr <- args$Tr
  z <- args$z
  sigma <- 2
  
  tol <- 1e-5

#### initialize ####
  
  # region prices
  p_r0 <- .sparseDiagonal(n=R,x=0.1)
  p_r1 <- .sparseDiagonal(n=R,x=1.0)
  
  # plant prices
  p_i0 <- .sparseDiagonal(n=N,x=0.1)
  p_i1 <- .sparseDiagonal(n=N,x=1.0)

  A.df <- A %>% s_to_df()
  G.df <- G %>% s_to_df()
  
  # region-plant demand shares
  LAM0 <- sparseMatrix(i=A.df$i, j=A.df$j, x=0.1, dims=c(max(A.df$i),max(A.df$j)))
  LAM1 <- sparseMatrix(i=A.df$i, j=A.df$j, x=1.0, dims=c(max(A.df$i),max(A.df$j)))
  
  # plant-plant demand shares
  GAM0 <- sparseMatrix(i=G.df$i, j=G.df$j, x=0.1, dims=c(max(G.df$i),max(G.df$j)))
  GAM1 <- sparseMatrix(i=G.df$i, j=G.df$j, x=1.0, dims=c(max(G.df$i),max(G.df$j)))

  # set objective: Frobenius norms of all matrices 
  # (note: the `p`s are diagonal matrices for technical reasons)
  obj = norm(LAM1-LAM0,"f") + norm(GAM1-GAM0,"f") + norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
  
### Solve for \Lambda and \Gamma ###   
  # while the difference between iterations is greater than tolerance
  while (obj > tol) {
#    print(str_c("obj: ",obj))
    
    # save newest iteration of parameters
    p_i0 <- p_i1
    p_r0 <- p_r1
    GAM0 <- GAM1
    LAM0 <- LAM1

    # calculate new p_r1
    m1 <- Tr %*% p_i0 %>% s_to_df() %>% mutate(x = x^(1-sigma)) %>% df_to_s()
    p_r1 <- rowSums(LAM0 * m1) %>% as("matrix") %>% s_to_df() %>% mutate(x=x^(1/(1-sigma))) %>% df_to_s() %>% to_sdiag()

    # calculate new eta ( = unit intermediate cost)
    m2 <- Ti %*% p_i0 %>% s_to_df() %>% mutate(x=x^(1-sigma)) %>% df_to_s()
    m2.1 <- rowSums(GAM0 * m2) %>% as("matrix") %>% s_to_df() %>% mutate(x=x^((1-beta)/(1-sigma))) %>% df_to_s()
    eta <- m2.1 %>% to_sdiag() * (C %>% to_sdiag())

    # calculate new p_i1    
    p_i1 <- eta / z
    
    # calculate new region-plant demand shares
    LAM1 <- (p_r1 %>% s_to_df() %>% mutate(x=x^(1-sigma)) %>% df_to_s()) %*% A * (Tr %*% p_i1 %>% s_to_df() %>% mutate(x=x^(sigma-1)) %>% df_to_s())
    
    # calculate new plant-plant demand shares
    m3 <- (eta %>% s_to_df() %>% mutate(x=x^(1-sigma)) %>% df_to_s())/(1-beta)
    r3 <- G * (Ti %*% p_i1 %>% s_to_df() %>% mutate(x=x^(sigma-1)) %>% df_to_s())
    GAM1 <- m3 %*% r3
    
    # normalize LAM1, GAM1
    LAM1 <- LAM1 / rowSums(LAM1)
    GAM1 <- GAM1 / rowSums(GAM1)

    # solve for w, normalize p and p?
    obj = norm(LAM1-LAM0,"f") + norm(GAM1-GAM0,"f") + norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
  }

### Return \Lambda and \Gamma
  # return the region-plant and plant-plant demand shares matrices
  return(list(lambda=LAM1,gamma=GAM1))
}