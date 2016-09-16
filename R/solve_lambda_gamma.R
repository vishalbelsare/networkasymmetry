library(dplyr)
library(tibble)
library(Matrix)
library(stringr)

# function to convert sparseMatrix to tbl_df.
s_to_df <- function(m) {
  m %>% as("sparseMatrix") %>% summary() %>% tbl_df()
}
## Region price equations
df_to_s <- function(tdf) {
  sparseMatrix(i=tdf$i,j=tdf$j,x=tdf$x,dims=c(max(tdf$i),max(tdf$j))) # dims here.
}
to_sdiag <- function(x) { # an Mx1 vector to sparse diagonal
  if (class(x)=="dgCMatrix") {
    return(.sparseDiagonal(n=length(x[,1]),x=x[,1]))
  } else if (is.vector(x)) {
    return(.sparseDiagonal(n=length(x),x=x))
  }
}

initialize_fakes <- function(R,N) {
  sigma <- 2
  # epsilon <- 2
  
  # keep track of plant-region id. could be sparse matrix. could be tbl_df. or both.
  ir <- tibble(j=1:N,x=1,i=sample(1:R,N,replace=TRUE)) %>% df_to_s()
  
  # should make sure each region has at least one demand parameter.
  A <- rsparsematrix(R,N,density=0.75, rand.x=runif) 
  G <- rsparsematrix(N,N,density=0.75, rand.x=runif)
  
  A.df <- A %>% s_to_df()
  G.df <- G %>% s_to_df()
  
  Ti <- sparseMatrix(i=G.df$i, j=G.df$j, x=1, dims=c(max(G.df$i),max(G.df$j))) * (1+matrix(rlnorm(N^2,0,0.1),nrow=N,ncol=N)) %>% as("sparseMatrix")
  Tr <- sparseMatrix(i=A.df$i, j=A.df$j, x=1, dims=c(max(A.df$i),max(A.df$j))) * (1+matrix(rlnorm(R*N,0,0.1),nrow=R,ncol=N)) %>% as("sparseMatrix")
  
  # use random sparsematrix here instead.
  beta <- runif(N,0.4,0.6)
  C <- beta^beta * (1-beta)^(beta-1)
  z <- rlnorm(N,0,1)
  
  return(list(A=A,beta=beta,C=C,G=G,ir=ir,sigma=sigma,Ti=Ti,Tr=Tr,z=z))
}

solve_lambda_gamma <- function(R=2,N=10,x=initialize_fakes(R,N)) { #A,beta,G,N,R,Ti,Tr,z) {
#  if (missing(x)) {
#    x <- initialize_fakes(R,N)
#  }
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
  
  A <- x$A
  beta <- x$beta
  C <- x$C
  G <- x$G
  ir <- x$ir
  sigma <- x$sigma
  Ti <- x$Ti
  Tr <- x$Tr
  z <- x$z

  tol=1e-5

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

  # set objective: Frobenius norms of all parameters.
  obj = norm(LAM1-LAM0,"f") + norm(GAM1-GAM0,"f") + norm(p_r1-p_r0,"f") + norm(p_i1-p_i0,"f")
  
  
### Solve for \Lambda and \Gamma ###   
  # while the difference between iterations is greater than tolerance
  while (obj > tol) {
    print(str_c("obj: ",obj))
    
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

### Return \Lambda and \Gamma ###  
  # return the region-plant and plant-plant demand shares matrices
  return(list(lambda=LAM1,gamma=GAM1))
}