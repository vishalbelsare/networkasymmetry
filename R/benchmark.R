########################################################################
# benchmark.R
# Function to benchmark A and G consistent with given parameters
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

benchmark <- function(R,N,args) {

  if (missing(args)) {
    args <- initialize_fake_links(R,N)
  }

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

  tol <- 1e-5 #/ sqrt(N)  # should be a function of N, because of normalization by sum(y1).

  A0 <- A1 <- Er
  G0 <- G1 <- En

  # initialize y0, y1, objective.
  y0 <- y1 <- rep_len(1,N)
  obj <- tol + 1

  while (obj > tol) {
    print(obj)
    # save last iteration's solutions
    A0 <- A1
    G0 <- G1
    y0 <- y1

    # xr <- -log(Er %*% exp(y0))
    # xn <- log(1-beta) - log(En %*% exp(y0))
    #
    # u <- (log(I) + xr) %>% as.matrix()
    # v <- (log(s) + xn) %>% as.matrix()
    #
    # y1 <- log(s) - log(t(Er) %*% exp(u) + t(En) %*% exp(v))
    # A1 <- ((exp(xr) %>% to_sdiag()) %*% Er) * (Er %*% (exp(y1) %>% to_sdiag()))
    # G1 <- ((exp(xn) %>% to_sdiag()) %*% En) * (En %*% (exp(y1) %>% to_sdiag()))

    # Intermediate calculations
    xr <- (Er %*% y0)^(-1)
    xn <- (1-beta) * (En %*% y0)^(-1)
#    xn[!is.finite(xn)] <- 0 # try to fix xn->infty problems.

    u <- I/(Er %*% y0) %>% as.matrix()
    v <- s/(En %*% y0) %>% as.matrix()

    # Update y1
    y1 <- s * (t(Er) %*% u + t(En) %*% v)^(-1)

    # Outer product gives A and G. a_{ri} = x^R_r * y_i, g_{ij} = x^N_i * y_j
    A1 <- ((xr %>% to_sdiag()) %*% Er) * (Er %*% (y1 %>% to_sdiag()))
    G1 <- ((xn %>% to_sdiag()) %*% En) * (En %*% (y1 %>% to_sdiag()))

    A1@x[!is.finite(A1@x)] <- 0
    G1@x[!is.finite(G1@x)] <- 0

    obj <- norm(A1-A0,"f") + norm(G1-G0,"f")
    print(str_c(norm(A1-A0,"f"),", G: ",norm(G1-G0,"f")))

    # if (!is.finite(obj)) {
    #   break
    # }
  }

  # Some ugly stuff to correct the elements that go wacko.
  # Find the ones that go to 0 or Inf or NaN, and pick a new one at random.
  # b <- rowSums(A1)
  # c <- colSums(A1)
  # bx <- b==0 | is.na(b)
  # cx <- c==0 | is.na(c)
  #
  # A1[bx, cx] <- matrix(runif(length(bx[bx])*length(cx[cx]),0,1),nrow=length(bx[bx]),ncol=length(cx[cx]))

  A1 <- A1 / rowSums(A1) # re-normalize.

  # b <- rowSums(G1)
  # c <- colSums(G1)
  # bx <- b==0 | is.na(b)
  # cx <- c==0 | is.na(c)
  #
  # G1[bx, cx] <- matrix(runif(length(bx[bx])*length(cx[cx]),0,1),nrow=length(bx[bx]),ncol=length(cx[cx]))

  G1 <- (1-beta) * G1 / rowSums(G1) # re-normalize.

  return(list(A=A1,G=G1))
}
