########################################################################
# initialize_functions.R
# Functions to generate fake parameters and data
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

initialize_fake_links <- function(R,N) {
  
  # Plant value-added share
  beta <- runif(N,0.3,0.7)
  
  # Plant output
  s <- rlnorm(N,0,1) * N 
  
  # Plant's region
  ir <- tibble(j=1:N,i=sample(1:R,N,replace=TRUE),x=1) %>% df_to_s()
  
  # Regional income
  I <- ir %*% (beta * s)
  
  # Non-zero edges of region-plant demand matrix
  Er <- rsparsematrix(R,N,density=0.01,rand.x=function(n) 1)

  while (min(colSums(Er))==0 | min(rowSums(Er))==0) {
    if (min(colSums(Er))==0) {
      q <- colSums(Er) == 0
      Er[,q] <- tibble(j=1:length(q[q]),
                       i=sample(1:R,length(q[q]),replace=TRUE),
                       x=1) %>% 
        df_to_s(dims=c(R,length(q[q])))
    }
    if (min(rowSums(Er))==0) {
      q <- rowSums(Er) == 0
      Er[q,] <- tibble(i=1:length(q[q]),
                       j=sample(1:N,length(q[q]),replace=TRUE),
                       x=1) %>% 
        df_to_s(dims=c(length(q[q]),N))
    }
    Er <- Er %>% summary() %>% tbl_df() %>% filter(x>0) %>% df_to_s(dims=c(R,N))
  }

  # Non-zero edges of plant-plant demand matrix
  En <- rsparsematrix(N,N,density=0.001,rand.x=function(n) 1)
  while ( min(colSums(En))==0 | min(rowSums(En))==0 | sum(diag(En))>0 ) {

    diag(En) <- 0
    # for columns:
    if (min(colSums(En))==0) {
      q <- colSums(En) == 0
      
      En[,q] <- tibble(j=1:length(q[q]),
                       i=sample((1:N)[!(1:N %in% j)],length(q[q]),replace=TRUE),
                       x=1) %>% 
        df_to_s(dims=c(N,length(q[q])))
    }
    if (min(rowSums(En))==0) {
      q <- rowSums(En) == 0
      En[q,] <- tibble(i=1:length(q[q]),
                       j=sample((1:N)[!(1:N %in% i)],length(q[q]),replace=TRUE),
                       x=1) %>% 
        df_to_s(dims=c(length(q[q]),N))
    }
    En <- En %>% summary() %>% tbl_df() %>% filter(x>0) %>% df_to_s(dims=c(N,N))
  }
  
  
  # Return generated parameters.
  return(list(beta=beta,Er=Er,En=En,I=I,ir=ir,s=s))
}

initialize_fakes <- function(R,N,args) {

  if (missing(args)) {
    args <- initialize_fake_links(R,N)
  } 
  
  # to debug:
  # rm(list=ls(all=TRUE))
  # gc()
  # source("R/initialize_functions.R")
  # source("R/benchmark.R")
  # source("R/helpers.R")
  # R <- 50
  # N <- 500
  # args <- initialize_fake_links(R,N)
  
  # Region-plant demand edges
  Er <- args$Er
  
  # Plant-plant demand edges
  En <- args$En
  
  # Region income
  I <- args$I
  
  # Plant's region
  ir <- args$ir
  
  # Value-added share
  beta <- args$beta
  
  # Output
  s <- args$s
  
  # Constant (varies by plant)
  C <- beta^(-beta) * (1-beta)^(beta-1)
  
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
  Ti <- sparseMatrix(i=G.df$i, j=G.df$j, x=1, dims=c(N,N)) * (1+matrix(rlnorm(N^2,-1,1),nrow=N,ncol=N)) %>% as("sparseMatrix")

  # Generate trade cost matrix, region-plant.
  Tr <- sparseMatrix(i=A.df$i, j=A.df$j, x=1, dims=c(R,N)) * (1+matrix(rlnorm(R*N,-1,1),nrow=R,ncol=N)) %>% as("sparseMatrix")

  # Return fake data; don't need non-zero edge matrices anymore.
  return(list(A=A,beta=beta,C=C,G=G,I=I,ir=ir,s=s,Ti=Ti,Tr=Tr,z=z))
}

