#' @param X a matrix of raw counts from samples in the reference group. 
#' Samples in rows and cell populations in columns
#' @param nDA a positive integer for the number of populations with differential abundant (DA)
#' between the reference group and simulated group
#' @param sample_size an integer vector of dimension 2 with the first and second value are the required 
#' sample size for simulated group 1 and 2, respectively.
#' @param rltv.effect a character c("random", "low", "median", "high") for the degree of relative 
#' difference of DA populations
#' @param frac.DA.samples a value between 0 and 1 for the fraction of samples where DA occurs. Default 1. 
#' @param seed seed
#' 
#' @return
#' A list of three items
#' \itemize{
#' \item{"X"}{A combined matrix of the reference and simulated counts}
#' \item{"group"}{A vector of indicator variable for the groups: 0=reference group amd 1=simulated group}
#' \item{"DApops"}{A vector of population names for DA populations}
#' \item{"DAsamples"}{A vector of samples where DA populations occured}
#' }
#' 
simCP <- function(X, frac.daPop=0.1, sample_size = NULL,  makeCompositional=TRUE,
                  rltv.effect = "random", weight=TRUE, mix_distance_group=3,  seed=3545){  
  
  set.seed(seed)
  
  # split samples into two groups
  if(is.null(sample_size)){
    n0 <- max(1, nrow(X)/2)
    n1 <- nrow(X)-n0
  }else{
    n0 <- sample_size[1]
    n1 <- sample_size[2]
    if(n0+n1 > nrow(X)) stop("the total sample size must be at mostn the rows of X")
  }
  
  n0.samples <- sample(seq_len(nrow(X)), n0)
  n1.samples <- sample(seq_len(nrow(X))[-n0.samples], n1)
  X0 <- X[n0.samples, ]
  X1 <- X[n1.samples, ]
  
  
  # number of DA populations
  K <- ncol(X)
  nDA = min(max(1, ceiling(frac.daPop*K)), floor(K/2))
  
  # rank populations according to their mean rank
  normX  <- t(apply(X, 1, function(x) as.numeric(clr(x))))
  normX0 <- t(apply(X0, 1, function(x) as.numeric(clr(x))))
  normX1 <- t(apply(X1, 1, function(x) as.numeric(clr(x))))
  colnames(normX) <- colnames(normX0) <-  colnames(normX1) <- colnames(X)
  
  dist_pp <- combn(colnames(normX), 2, function(x){
    y1=normX[, x[1]]
    y2=normX[, x[2]]
    d = sqrt(sum((y1-y2)^2)) 
    data.frame(ppl1=x[1], ppl2=x[2], comb_ppl=paste0(x, collapse = "_&_"), d=d)
  }, simplify = FALSE) %>% bind_rows() %>%
    arrange(desc(d)) 
   
  # if(rltv.effect != "random"){
  #   q.prob  <- seq(0, 1, length.out = mix_distance_group+1)
  #   q.break <- quantile(dist_pp$d, probs =  q.prob, names = FALSE) 
  #   mix_dist_group = cut(dist_pp$d, breaks = q.break, include.lowest=TRUE) 
  #   dist_pp$mix_dist_group =mix_dist_group
  # }else{
  #   dist_pp$mix_dist_group = NA   
  # }
  
  
  # ggplot(dist_pp, aes(x=ppl1, y=ppl2, fill=d))+
  #   geom_tile()+ 
  #   theme_bw()
  
  # ggplot(dist_pp, aes(x=ppl1, y=ppl2, fill=d))+
  #   geom_tile()+
  #   facet_grid(~mix_dist_group)+
  #   theme_bw()
  # 
  # ggplot(dist_pp, aes(x=d))+
  #   geom_density()+
  #   facet_grid(~mix_dist_group)+
  #   theme_bw()

  # plot(density(dist_pp$d))
  
  #popl.rank <- rank(apply(normX1, 2, function(x) mean(x)))
  
  # # select samples where DA takes place
  # h <- min(round(n1*frac.DA.samples), n1)
  # DA.samples <- sample(seq_len(nrow(X1)), h)
  
  # select DA populations
  DApops <- matrix(NA, nrow=nDA, ncol=2)
  
  if(rltv.effect=="random"){
    slcted_cmbn_indx <- sample(1:nrow(dist_pp), size = nDA)
    DApops[,1] <-dist_pp[slcted_cmbn_indx, "ppl1"]
    DApops[,2] <-dist_pp[slcted_cmbn_indx, "ppl2"]
  }else if(rltv.effect=="lowDiffWeighted.random"){
    d = dist_pp$d
    # de = exp(-sqrt(d))
    # 
    # min.de = min(de)*0.99
    # max.de = max(de)*1.01
    dist_pp$whgt = sum(d)/d
    
    for(i in 1:nDA){
      uupdate.ppls <- dist_pp[!(dist_pp$ppl1 %in% DApops[,1]),]
      slcted_cmbn_indx <- sample(1:nrow(uupdate.ppls), size = 1, prob = uupdate.ppls$whgt)
      DApops[i,1] <- uupdate.ppls[slcted_cmbn_indx, "ppl1"]
      DApops[i,2] <- uupdate.ppls[slcted_cmbn_indx, "ppl2"]
    }
    # slcted_cmbn_indx <- sample(1:nrow(dist_pp), size = nDA, prob = whgt)
    # DApops[,1] <-dist_pp[slcted_cmbn_indx, "ppl1"]
    # DApops[,2] <-dist_pp[slcted_cmbn_indx, "ppl2"]
    
    #check for d
  }else if(rltv.effect=="highDiffWeighted.random"){
    d = dist_pp$d
    # de = exp(sqrt(d))
    # 
    # min.de = min(de)*0.99
    # max.de = max(de)*1.01
    # whgt = (de-min.de)/(max.de-min.de)
    dist_pp$whgt = d/sum(d)
    
    for(i in 1:nDA){
      uupdate.ppls <- dist_pp[!(dist_pp$ppl1 %in% DApops[,1]),]
      slcted_cmbn_indx <- sample(1:nrow(uupdate.ppls), size = 1, prob = uupdate.ppls$whgt)
      DApops[i,1] <- uupdate.ppls[slcted_cmbn_indx, "ppl1"]
      DApops[i,2] <- uupdate.ppls[slcted_cmbn_indx, "ppl2"]
    } 
    # slcted_cmbn_indx <- sample(1:nrow(dist_pp), size = nDA, prob = whgt)
    # DApops[,1] <-dist_pp[slcted_cmbn_indx, "ppl1"]
    # DApops[,2] <-dist_pp[slcted_cmbn_indx, "ppl2"]
  }
  
  
  
  # create DA populations
  pps <- colnames(X1)
  X1tmp <- X1
  
  for(i in seq_len(nrow(DApops))){  
    ui <- as.numeric(unlist(as.vector(X1[, DApops[i, 1]])))
    vi <- as.numeric(unlist(as.vector(X1[, DApops[i, 2]])))
    
    X1tmp[, DApops[i, 1]] <- vi 
    if(makeCompositional){
      d <- vi - ui
      X1_woi <- X1tmp[, which(pps != DApops[i, 1])]
      if(weight){
        wws <- t(t(X1_woi) %*% diag(1/rowSums(X1_woi)))  
        dmat <- d*wws 
        #mwws <- colMeans(wws) 
        #dmat <- d%*%t(mwws) 
      }else{
        dmat <- matrix(d/(K-1), length(d), length(pps)-1, byrow = FALSE)
      } 
      X1tmp[, which(pps != DApops[i, 1])] <- round(X1tmp[, which(pps != DApops[i, 1])] - dmat) 
      X1tmp[X1tmp<0] <- 0
    } 
  }  
  
  # combine and return
  Xsim <- bind_rows(X0, X1tmp)
  group = rep(0:1, c(n0, n1))
  list(X=Xsim, group=group, DApops = DApops)
}





















