#library(ggplot2)
library("reshape2")

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
normalize = function(frame, var, id) {
  temp=aggregate(var ~ id, frame, FUN=sum)
  #print(id)
  head(temp)
  #  temp=rename(temp, string(var),"gammasum")
  data=merge(data, temp, by=id)
  data$var.x = data$var.x / data$var.y
  data$var.y=NULL
  return(data)
}

N=5
nonzero=ceiling(N^1.5)
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
  i = sample(1:N, nonzero, replace = TRUE),
  j = sample(1:N, nonzero, replace = TRUE), # kk.
  gamma = exp(rnorm(nonzero)), #sample(1:N, nonzero, replace = TRUE),
  tau = exp(rnorm(nonzero)), #sample(1:N, nonzero, replace = TRUE),
#  beta = sample(1:N, nonzero, replace = TRUE),
  g0 = 0,
  g1 = 1)

p0=data.frame(i,p0)
p1=data.frame(i,p1)
beta=data.frame(i,beta)
#data=merge(data, p0, by="i")
#data=merge(data, p1, by="i")

pj=data.frame(j,p1$p1)
pj=rename(pj,"p1.p1","pj")
data=merge(data, pj, by="j")
data=merge(data, pj, by="j")

# calculate p1.
# need the temp, depends on pj and gamma.
#data$p1 = 

head(data)

# and labour is in there somewhere.

# learn to normalize. sum within i, normalize gamma?
#tapply(data$gamma1, data$i, FUN=sum)


x=normalize(data,"gamma","i")
head(x)

#xxxxx=84903291830290)(*(@)#!(
#                       
#                       
# # hmm, kind of annoying.
# temp=aggregate(gamma ~ i, data, sum)
# temp=rename(temp,"gamma","gammasum")
# data=merge(data, temp, by="i")
# data$gamma = data$gamma / data$gammasum
# data$gammasum=NULL
# # drop columns I don't need, or merge properly. need to replace columns in data.
# # doesn't look like there is one. must drop before merge.
# 
# # same idea, merge on price on j.
# #temp=aggregate(p1 ~ i, data, mean) # but would rather just get "first" with really big ones, ya?
# # just use p1 itself.
# temp=p1
# temp=rename(temp,"i","j")
# temp=rename(temp,"p1","pj")
# # drop old pj? or actuall don't have to redo this until loop.
# #data=data[ , !names(data) %in% c("pj")]
# data$pj=NULL
# data=merge(data, temp, by="j")
# head(data)
# 
# # now get gamma * pj^(1-sigma) * (1 + trade cost)^(elasticity)
# data$temp=data$gamma*data$pj^(1-sigma)*(1+data$tau)^(-theta)
# x=aggregate(temp ~ i, data, sum)
# data$temp=NULL
# data=merge(data,x,by="i")
# # normalize that. multiply by beta.
# 
# 
# head(data)
# 
# 
