#library(ggplot2)
rm(list=ls(all=TRUE))
#library("reshape2")

setwd("/Users/jessetweedle/Documents/github/networkasymmetry")

#set.seed(11)
# suppose that's gamma.
# need to solve for p.
# need data on distance for each nonzero thing in gamma.
# sigma.
# w. 
# vector of betas.
# alphas.
# distance to markets.
# need some practice with data frame. and probably won't use Matrix for now, too hard to get on server.
# should try to do it in logs too.

rename = function(frame, old, new) {
  names(frame)[names(frame)==old]=new  
  head(frame)
  return(frame)
}

## a function that calculates something and normalizes it.
## also keep a separate copy of prices, ya?

N=5
nonzero=N^2 #ceiling(N^1.5)
sigma=2
epsilon=2
theta=4

w0=0
w1=1

i=seq(N)
j=seq(N)
p0=numeric(N)
p1=numeric(N)+1
beta=numeric(N)+1

z=exp(rnorm(N)) # k. useful. productivities.

data = data.frame(
  i = rep(1:N,N), #sample(1:N, nonzero, replace = TRUE),
  j = sample(1:N, nonzero, replace = TRUE), # kk.
  gamma = exp(rnorm(nonzero)), #sample(1:N, nonzero, replace = TRUE),
  tau = exp(rnorm(nonzero)), #sample(1:N, nonzero, replace = TRUE),
#  beta = sample(1:N, nonzero, replace = TRUE),
  g0 = 0,
  g1 = 1)

# say i=1 is labour?
# should rewrite all my equations and stuff in CES, and also add distance.

data = aggregate(data, by=list(data$i,data$j), FUN=mean, na.rm=TRUE)
data$Group.1=NULL
data$Group.2=NULL
#kk. total number of plants not fully represented.
# remove i=j.
data=data[data$i!=data$j,]
 
p0=data.frame(i,p0)
p1=data.frame(i,p1)
#beta=data.frame(i,beta)

pj=data.frame(j,p1$p1)
pj=rename(pj,"p1.p1","pj")
data=merge(data, pj, by="j")

head(data)

# and labour is in there somewhere.

# learn to normalize. sum within i, normalize gamma?
#tapply(data$gamma1, data$i, FUN=sum)

normalize = function(frame, var, id) {

  temp=aggregate(reformulate(c(id), response=var), frame, FUN=sum)
  data=merge(data, temp, by=id)

  x=paste0(var,".x")
  y=paste0(var,".y")

  data[,x] = data[,x] / data[,y]
  data[,y]=NULL
  data=rename(data,x,var)
  return(data)
}

beta=0.5
#pconstant=

data=normalize(data,"gamma","i")
# but can just do one and then normalize.
data$p1=data$gamma*data$pj^(1-sigma)*(1+data$tau)^(theta*(1-sigma))
p1$p1=(beta^(-beta)*(1-beta)^(beta-1)*w1)*z^(-1)
p1=merge(p1,aggregate(reformulate(c("i"), response="p1"), data, FUN=sum),by="i", all.x=TRUE)
data$p1=NULL
p1$p1=p1$p1.x * p1$p1.y^((1-beta)/(1-sigma))
p1=p1[,c("i","p1")]

# then merge back on as pj
data$pj=NULL
pj=rename(p1,"p1","pj")
pj=rename(pj,"i","j")
data=merge(data, pj, by="j")
data$g1=data$gamma*data$pj^(1-sigma)
data$pj=NULL
data=normalize(data,"g1","i")
data=data[order(data$i,data$j),]

# beta. labour. but first get loop? nah.
