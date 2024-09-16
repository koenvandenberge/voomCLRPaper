library(DirichletReg) 
library(tidyverse)

### for testing:
# n0=5
# n1=n0
# K=10
# frac.daPop = 0.25
# abs.treat.effect.range = c(0.5, 2)
# intercept = 1
# repeated=FALSE
# phi_0 = 0.99
# phi_1 = phi_0
# parent.size=1e5
# scale.dirichlet = 1
# B0t = 0
# B1t = B0t
# count.out.format = "long"
# seed=round(runif(1, 1, 1e6))

DM_simCP <- function(sample_size=c(10, 10), 
                           K=10,  
                           frac.daPop = 0.1,
                           #abs.treat.effect.range = c(0, 10), 
                           # abs.treat.effect.mean = 1,
                           # abs.treat.effect.sd = 2,
                           mean.intercept = 1, sd.intercept = 0.1, sd.slope = 2, 
                     
                           parent.size=rep(1e5, 2), 
                           scale.dirichlet = 1,   
                           makeCompositional=TRUE,
                           seed=round(runif(1, 1, 1e6))){
   
  
  set.seed(seed)
  
  # sample size
  n0 <- sample_size[1]
  n1 <- sample_size[2]
   
  
  # DA populations
  n.daPop = min(max(1, ceiling(frac.daPop*K)), K)
  which.pop.da = sample(K, n.daPop)
  
  
  # parent size per sample
  if(length(parent.size)==2){
    N0 <- rep(parent.size[1], n0)
    N1 <- rep(parent.size[1], n1)
  }else{
    N0 <- parent.size[1:n0]
    N1 <- parent.size[(n0+1):(n0+n1)]
  }
  
  # effect size for DA populations
  # delt = runif(n.daPop, abs.treat.effect.range[1], 
  #              abs.treat.effect.range[2])*sample(c(-1, 1), n.daPop, replace = TRUE)
  delt = rnorm(n.daPop, 0, sd.slope)  
  B1 = numeric(K)
  B1[which.pop.da] = delt 
  
  m0 = rnorm(n=K, mean=mean.intercept, sd=sd.intercept)
  m1 = m0 + B1 
  
   
  if(makeCompositional){
    ## dirichlet parameters 
    dr.alpha0   = exp(m0)*scale.dirichlet
    dr.alpha1   = exp(m1)*scale.dirichlet
    
    ## multinomial probabilities
    pi0   = rdirichlet(n0, dr.alpha0)
    pi1   = rdirichlet(n1, dr.alpha1)
    
    # simulate counts
    count0    = t(sapply(1:n0, function(i){
      rmultinom(1, size = N0[i], prob = pi0[i,])
    }))
    count1    = t(sapply(1:n1, function(i){
      rmultinom(1, size = N1[i], prob = pi1[i,])
    }))  
  }else{
    # beta parameters
    mu0 <- exp(m0)/(1+exp(m0))
    mu1 <- exp(m1)/(1+exp(m1))
    
    # simulate binomial probabilities
    vv=n0+n1
    pi0   = sapply(1:K, function(k){
      rbeta(n0, shape1 = vv*mu0[k], shape2 = vv*(1-mu0[k]))
    }) 
    pi1   = sapply(1:K, function(k){
      rbeta(n1, shape1 = vv*mu1[k], shape2 = vv*(1-mu1[k]))
    })
    
    # simulate counts
    count0    = t(sapply(1:n0, function(i){
      sapply(1:K, function(k){
        rbinom(1, size = N0[i], prob =  pi0[i,k]) 
      })
    }))
    count1    = t(sapply(1:n1, function(i){
      sapply(1:K, function(k){
        rbinom(1, size =  N1[i], prob = pi1[i,k]) 
      })
    }))  
  }
  
  
  # combine and return
  count.all = bind_rows(as.data.frame(count0), as.data.frame(count1))
  sample_id = 1:(n0+n1)
  rownames(count.all) <- sample_id
  group = rep(0:1, c(n0, n1))
  colnames(count.all) =  paste0("ppl_", 1:K)
  dif.abundant.comps=paste0("ppl_", which.pop.da)
  
  list(X=count.all, group=group, DApops = matrix(dif.abundant.comps, ncol=1),
       true.ES=B1, dr.alpha0=dr.alpha0, dr.alpha1=dr.alpha1, N0=N0, N1=N1)
}

# simdat1 = simComposData2(count.out.format="longer", n0 = 100, scale.dirichlet=1,
#                          phi_0 = 0.99, B0t = 0, B1t = 2, abs.treat.effect.range = c(1.75, 2))
# simdatDat1 = simdat1$count %>%
#   group_by(sample, x, t) %>%
#   mutate(frac = count/N, CLR = log((count)/(prod(count)^(1/10))))
# 
# simdat1$dif.abundant.comps
# 
# ggplot(simdatDat1, aes(x=as.factor(t), y=frac, colour=as.factor(x)))+
#   geom_boxplot()+
#   facet_wrap(~population)+
#   theme_bw()


# simdatDat1_t0 = simdatDat1 %>% subset(t==0) %>%
#   dplyr::select(c(sample, population, x, count, CLR)) %>%
#   dplyr::rename(sample0=sample, population0=population,
#                 x0=x, count0=count, CLR0=CLR, t0=t)
# 
# simdatDat1_t1 = simdatDat1 %>% subset(t==1) %>%
#   dplyr::select(c(sample, population, x, count, CLR)) %>%
#   dplyr::rename(sample1=sample, population1=population,
#                 x1=x, count1=count, CLR1=CLR, t1=t)
# 
# simdatDat1_v = cbind(simdatDat1_t0, simdatDat1_t1)
# 
# ggplot(simdatDat1_v, aes(x=CLR0, y=CLR1))+
#   geom_point()+
#   facet_grid(x0~population0)+
#   theme_bw()
