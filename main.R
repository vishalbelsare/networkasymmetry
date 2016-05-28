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

N=5
nonzero=ceiling(N^1.5)

w0=0
w1=1

i=seq(N)
j=seq(N)
p0=numeric(N)
p1=numeric(N)+1

z=exp(rnorm(N)) # k. useful. productivities.

data = data.frame(
  i = sample(1:N, nonzero, replace = TRUE),
  j = sample(1:N, nonzero, replace = TRUE), # kk.
  gamma = sample(1:N, nonzero, replace = TRUE),
  distance = sample(1:N, nonzero, replace = TRUE),
  g0 = 0,
  g1 = 1)

#data$gamma*p0 # it's really within an i.
data[data$i == 1, ] # hmm. might just say p0 p1 pj? or get used to merging?

p0=data.frame(i,p0)
p1=data.frame(i,p1)

data=merge(data, p0, by="i")
data=merge(data, p1, by="i")
pj=data.frame(j,p1$p1)
names(pj)[names(pj)=="p1.p1"] <- "pj"
#names(pj)[names(pj)=="p1"] <- "pj"
data=merge(data, pj, by="j")