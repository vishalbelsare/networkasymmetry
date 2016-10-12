########################################################################
# benchmark_xy.R
# Function to benchmark A and G consistent with given parameters
# License: MIT

# ""
# Jesse Tweedle
# Sep 15, 2016
########################################################################

benchmark_xy <- function(R,N,args) {

  if (missing(args)) {
    args <- initialize_fake_links(R,N)
  }

  Er <- args$Er # Region-plant demand edges
  En <- args$En # Plant-plant demand edges
  I <- args$I # Region income

  ir <- args$ir # Plant's region
  beta <- args$beta # Value-added share
  s <- args$s # Output

  tol <- 1e-10 / sqrt(N)  # should be a function of N, because of normalization by sum(y1).

  xr0 <- xr1 <- rep_len(1,R) # initialize y0, y1, objective.
  xn0 <- xn1 <- rep_len(1,N) # initialize y0, y1, objective.
  y0 <- y1 <- rep_len(1,N) # initialize y0, y1, objective.
  obj <- tol + 1

  while (obj > tol) {
    # save last iteration's solutions
    y0 <- y1
    xr0 <- xr1
    xn0 <- xn1

    # Intermediate calculations
    xr1 <- (Er %*% y0)^(-1)
    xn1 <- (1-beta) * (En %*% y0)^(-1)

    u <- I/(Er %*% y0) %>% as.matrix()
    v <- s/(En %*% y0) %>% as.matrix()

    # Update y1
    y1 <- s * (t(Er) %*% u + t(En) %*% v)^(-1)
    y1 <- y1/sum(y1)
    
    obj <- sqrt(sum((y1-y0)^2)) # + sum((xn1-xn0)^2) + sum((xr1-xr0)^2))
    print(str_c("y: ",obj))
  }

  # now we getting somewhere. this is quick, when matrix is dense.
  # sx <- (outer(xr1,y1) * Er) %>% t() %*% I + (En * outer(xn1,y1)) %>% t() %*% s
  # (s[(R+1):N]/sx[(R+1):N])  %>% as.vector() %>% summary()
  # plot(s[(R+1):N] %>% log(), sx[(R+1):N] %>% log())
  # tibble(s=s,sx=sx[,1]) %>% mutate(d=s/sx) %>% View()
  return(list(xr=xr1,xn=xn1,y=y1))
}
