########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# ""
# Jesse Tweedle
# Sep 15, 2016
########################################################################

library(dplyr)
library(tibble)
library(Matrix)
library(stringr)
library(tidyr)
library(ggplot2)

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

rm(list=ls(all=TRUE))
# Sometimes R or RStudio (?) doesn't automatically clean removed objects. I think.
# That's a problem with a NxN matrix when N really big. So clean just in case:
gc()

# Tip from http://jeromyanglim.tumblr.com/post/33418725712/how-to-source-all-r-files-in-a-directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir(paste0(getwd(),"/R"))


set.seed(10) # set.seed(10), R=2, N=300 screws up.
fakes <- TRUE
if (fakes) {
  R <- 2
  N <- 10
  sigma <- 2 # here.
  eta <- epsilon <- sigma

  print("fake edges")
  argsx <- initialize_fake_links(R,N)
  
  print("fake shares")
  args <- initialize_fakes(R,N,args=c(sigma=sigma,argsx))

} else {
  # Load data.
  # Get E, etc., and benchmark to get A, G, using firm-network-lasso code.
}

print("solve for lambda, gamma; use A and G as initial guesses?")
lg <- solve_lambda_gamma(R,N,args=c(sigma=sigma,args))

print("solve for s, A, G") # should try to use existing s, A, G, as intial arguments.
solved <- solve_v(R,N,args=c(sigma=sigma,args,lg))
plot(args$s %>% log(),solved$v[,1] %>% log())
print("done")

xxx

# counterfactuals
#^((beta-1)/(1-epsilon))
xlp <- rowSums(lg$lambda)^(-1) %>% to_sdiag()
lp_aug <- xlp %*% lg$lambda
xgp <- rowSums(lg$gamma)^(-1) %>% to_sdiag()
gp_aug <- xgp %*% lg$gamma
#^((beta-1)/(1-eta))
z_aug <- args$z * rowSums(lg$gamma)^((1-args$beta)/(eta-1))

# htf can I keep love-of-variety constant? Adding a million products to production function
# want the augmented bit to stay constant. So go back to production function, not the other thing.
# Normalize gamma, or gamma 1/N.

# 1. set everything but z to 1 (or 1/N)
ext_args <- list(
  beta=rep_len(0.5,N), #args$beta,
  eta=eta,
  epsilon=epsilon,
  lambda=(1/N),
  gamma=(1/N),
  z=z_aug #* (1/N)^((1-args$beta)/eta)
)

z_ext <- solve_v_dense_network(R,N,args=ext_args)

# 2. set everything but gamma, lambda to 1
Tip <- args$Ti
Trp <- args$Tr
Tip@x <- rep_len(1,length(args$Ti@x))
Trp@x <- rep_len(1,length(args$Tr@x))

demand_args <- list(
  beta=args$beta,
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  s=args$s,
  sigma=sigma,
  Ti=Tip,
  Tr=Trp,
  v=solved$v,
  lambda=lg$lambda,
  gamma=lg$gamma,
  z=rep_len(mean(args$z),N) # ok, so this does work. only question is whether it matters in the data. extensive margin shouldn't matter here, yea? why does productivity matter so much here? dunno. hope it not like this in the data.
)

d_c <- solve_v(R,N,args=demand_args)

beta_args <- list(
  beta=rep_len(0.5,N),
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  s=args$s,
  sigma=sigma,
  Ti=args$Ti,
  Tr=args$Tr,
  v=solved$v,
  lambda=lg$lambda,
  gamma=lg$gamma,
  z=args$z
)

beta_c <- solve_v(R,N,args=beta_args)

augment_args <- list(
  beta=args$beta,
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  s=args$s,
  sigma=sigma,
  Ti=args$Ti,
  Tr=args$Tr,
  v=solved$v,
  lambda=lp_aug,
  gamma=gp_aug,
  z=z_aug
)

aug_c <- solve_v(R,N,args=augment_args)
#solved_z_aug_prime <- solve_v(R,N,args=c(invariant,list(lambda=lp_aug,gamma=gp_aug,z=rep_len(1,N))))

# Higher order interconnections---set demand network to constant by j.
# also extensive margin things? To get there, I'll have to re-write solve_v_dense.

#lp <- lg$lambda
#lp@x <- rep_len(1,length(lp@x))
#lp@x <- rep_len(row)
#lp <- (N / R * rowSums(lp)^(-1) %>% to_sdiag()) %*% lp %*% (colSums(lg$lambda) %>% to_sdiag())
lamd <- colSums(lp_aug)
gamd <- colSums(gp_aug)
lamd <- lamd / sum(lamd)
gamd <- gamd / sum(gamd)
#rowSums(lg$lambda)

#lp <- lp %*% (colSums(lg$lambda) %>% to_sdiag())

# gp <- lg$gamma
# gp@x <- rep_len(1,length(gp@x))
# #gp <- gp %*% (colSums(lg$gamma) %>% to_sdiag())
# gp <- (rowSums(gp)^(-1) %>% to_sdiag()) %*% gp %*% (colSums(lg$gamma) %>% to_sdiag())
#(rowSums(gp) - rowSums(lg$gamma)) %>% summary()

#gamma_ext <- rowSums(lg$gamma)

higher_order_args <- list(
  beta=args$beta,
  eta=eta,
  epsilon=epsilon,
  lambda=lamd,
  gamma=gamd,
  z=args$z
)

ho_c <- solve_v_dense_network(R,N,args=higher_order_args)

dat <- tibble(
  z=args$z,
  # beta=args$beta,
  s=args$s,
  Data=solved$v[,1],
  Productivity=z_ext$v,
  Augmented=aug_c$v[,1],
  # Beta=beta_c$v[,1],
  "Higher Order"=ho_c$v,
  Demand=d_c$v[,1]
)

#dat <- dat %>% rownames_to_column() %>% mutate(vv=sum(vp),vp=vp/vv) %>% select(-vv) %>% gather(type,size,v:vp)
dat <- dat %>% rownames_to_column() %>% gather(type,size,Data:Demand)

p <- ggplot(dat) + # %>% filter(type=="v" | type=="vext")) +
  stat_density(position="dodge",trim=TRUE,geom="line",aes(x=size,colour=type)) +
  scale_x_log10() + labs(x="Size (normalized, log scale)",y="Density")

px <- p + theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.justification=c(0,1),
        legend.position=c(0,1)
        )
px

# see
ggplot(dat) + geom_point(aes(x=s,y=size,colour=type),alpha=0.3) + scale_x_log10() + scale_y_log10()

print(
  dat %>% group_by(type) %>% mutate(size=size^2) %>% summarize(size=sum(size)) %>% mutate(size=sqrt(size))
)
