initialize_fake_links <- function(R,N) {
  # Plant value-added share
  beta <- runif(N,0.3,0.7)
  # Plant output
  s <- rlnorm(N,0,1)
  # Plant's region
  ir <- tibble(j=1:N,i=sample(1:R,N,replace=TRUE),x=1) %>% df_to_s()
  # Regional income
  I <- ir %*% (beta * s)
  # Non-zero edges of region-plant demand matrix
  Er <- rsparsematrix(R,N,density=0.75,rand.x=function(n) 1)
  # Non-zero edges of plant-plant demand matrix
  En <- rsparsematrix(N,N,density=0.75,rand.x=function(n) 1)
  
  # Return generated parameters.
  return(list(beta=beta,Er=Er,En=En,I=I,ir=ir,s=s))
}

initialize_fakes <- function(R,N) {

  args <- initialize_fake_links(R,N)
  
  # Region-plant demand edges
  Er <- args$Er
  # Plant-plant demand edges
  En <- args$En
  # Region income
  I <- args$I
  
  # Plant characteristics:

  # Plant's region
  ir <- args$ir
  # Value-added share
  beta <- args$beta
  # Output
  s <- args$s
  # Constant (varies by plant)
  C <- beta^beta * (1-beta)^(beta-1)
  # Productivity
  z <- rlnorm(N,0,1)
  
  # Solve for consistent region-plant and plant-plant demand share matrices.
  demand <- benchmark(R,N,args=args)
  
  # Region-plant demand shares
  A <- demand$A
  # Plant-plant demand shares.
  G <- demand$G

  # Change these to data frames.
  A.df <- A %>% s_to_df()
  G.df <- G %>% s_to_df()

  # Generate trade cost matrix, plant-plant.  
  Ti <- sparseMatrix(i=G.df$i, j=G.df$j, x=1, dims=c(max(G.df$i),max(G.df$j))) * (1+matrix(rlnorm(N^2,0,0.1),nrow=N,ncol=N)) %>% as("sparseMatrix")
  # Generate trade cost matrix, region-plant.  
  Tr <- sparseMatrix(i=A.df$i, j=A.df$j, x=1, dims=c(max(A.df$i),max(A.df$j))) * (1+matrix(rlnorm(R*N,0,0.1),nrow=R,ncol=N)) %>% as("sparseMatrix")

  # Return fake data; don't need non-zero edge matrices anymore.
  return(list(A=A,beta=beta,C=C,G=G,I=I,ir=ir,Ti=Ti,Tr=Tr,z=z))
}

