########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
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
R <- 50
N <- 1000
sigma <- 2 # here.

print("fake edges")
argsx <- initialize_fake_links(R,N)

print("fake shares")
args <- initialize_fakes(R,N,args=c(sigma=sigma,argsx))

print("solve for lambda, gamma; use A and G as initial guesses?")
lg <- solve_lambda_gamma(R,N,args=c(sigma=sigma,args))

print("solve for s, A, G") # should try to use existing s, A, G, as intial arguments.
solved <- solve_v(R,N,args=c(sigma=sigma,args,lg))

print("done")

# counterfactuals

# 1. set all z to 1

# 2. set all gamma, lambda to 1

# 3. set all tau to 1

# 4. all intensive margin counterfactuals

# 5. extensive margin lambda and tau.

# returns v, p_r, p_i, A, G.

v <- solved$v
z <- args$z
lambda <- lg$lambda
gamma <- lg$gamma

# beta, C, ir, lambda, gamma, Ti, Tr, v, z
invariant <- list(
  beta=args$beta,
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  sigma=sigma,
  Ti=args$Ti,
  Tr=args$Tr,
  v=solved$v
)

initial <- list(
  lambda=lg$lambda,
  gamma=lg$gamma,
  z=args$z
)


current <- list(
  lambda=lambda,
  gamma=gamma,
  z=rep_len(1,N)
)

ext_args <- list(
  beta=args$beta,
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  sigma=sigma,
  v=solved$v,
  z=args$z
)

# ok, decent.
solved_z_ext <- solve_v_dense_network(R,N,args=c(sigma,ext_args))

solved_z_prime <- solve_v(R,N,args=c(invariant,current))

lp = lambda
lp@x = rep_len(1,length(lp@x))
gp = gamma
gp@x = rep_len(1,length(gp@x))
solved_d_prime <- solve_v(R,N,args=c(invariant,list(lambda=lp,gamma=gp,z=args$z)))

xlp = rowSums(lambda)^(-1) %>% to_sdiag()
lp_aug = xlp %*% lambda
xgp = rowSums(gamma)^(-1) %>% to_sdiag()
gp_aug = xgp %*% gamma
solved_d_aug_prime <- solve_v(R,N,args=c(invariant,list(lambda=lp_aug,gamma=gp_aug,z=args$z)))

solved_z_aug_prime <- solve_v(R,N,args=c(invariant,list(lambda=lp_aug,gamma=gp_aug,z=rep_len(1,N))))


# and another one that sets z and z_aug to 0.

Tip=args$Ti
Trp=args$Tr
Tip@x=rep_len(1,length(args$Ti@x))
Trp@x=rep_len(1,length(args$Tr@x))
solved_t_prime <- solve_v(R,N,args=c(list(
  beta=args$beta,
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  sigma=sigma,
  Ti=Tip,
  Tr=Trp,
  v=solved$v
),initial))

solved_all_prime <- solve_v(R,N,args=c(list(
  beta=args$beta,
  C=args$C,
  ir=args$ir,
  p_i=solved$p_i,
  p_r=solved$p_r,
  sigma=sigma,
  Ti=Tip,
  Tr=Trp,
  v=solved$v
),list(
  lambda=lp,
  gamma=gp,
  z=rep_len(1,N)
)))

dat <- tibble(
  v=solved$v[,1],
  vz=solved_z_prime$v[,1],
  vd=solved_d_prime$v[,1],
  vt=solved_t_prime$v[,1],
  vdaug=solved_d_aug_prime$v[,1],
  vzaug=solved_z_aug_prime$v[,1],
  vext=solved_z_ext$v,
  vall=solved_all_prime$v[,1]
)

#dat <- dat %>% rownames_to_column() %>% mutate(vv=sum(vp),vp=vp/vv) %>% select(-vv) %>% gather(type,size,v:vp)
dat <- dat %>% rownames_to_column() %>% gather(type,size,v:vall)

#ggplot(dat , aes(size,colour=type)) + stat_density(geom="line") + scale_x_log10()

p <- ggplot(dat %>% filter(type!="vall")) +
  stat_density(position="dodge",geom="line",aes(x=size,colour=type)) +
  scale_x_log10() + #scale_colour_brewer(palette="Dark2") +
  # scale_colour_brewer(palette="Dark2", name="Counterfactual",
  #                       breaks=c("v","vd","vt","vz","vall","vext"),
  #                       labels=c("Data", "Demand","Geography","Productivity","All","Extensive")) +
  labs(x="Size (normalized, log scale)",y="Density")

px <- p + #theme_bw() +
  theme(panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.justification=c(0,1),
        legend.position=c(0,1)
        )

px

print(
dat %>% group_by(type) %>% mutate(size=size^2) %>% summarize(size=sum(size)) %>% mutate(size=sqrt(size))
)
#ggsave("plot-distributions.png", px, width=7, height=5, device="png")


# so, get correlations in data between productivity, gamma, lambda, etc. correlation between productivity of customers and demand parameters?
# idea:
# so, compare lambda, gamma, z. in data. gamma -> outdegree, etc.
#
# lambda <- lg$lambda
# gamma <- lg$gamma
lambda <- lp_aug
gamma <- gp_aug
# outdegree <-
outdegree <- colSums(lambda) + colSums(gamma)
z_aug <- rowSums(gamma) #apply(gamma,1,FUN=mean) #m/lapply or whatever.
# then compare to z, outd, v, z_aug, and?
# correlations, dist, scatters.

cor(z,outdegree)
cor(z,z_aug)

# ok.



# definitely negative in fake data, since for a given s, drawing high z means inferring low values for other parameters.
