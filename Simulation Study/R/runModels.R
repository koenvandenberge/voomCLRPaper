# helper function to run models

library(DirichletReg)
library(compositions)
library(VGAM)
library(lmtest)
library(lme4)
library(lmerTest)
library(MASS)
library(glmmTMB)
library(limma)
library(voomCLR) #
#system("cd \"/Users/aassefa/Library/CloudStorage/OneDrive-JNJ/JNJ (ATA)/Projects/FlowCytometryand CyTOF data analysis methods/voomCLR/\"; R CMD INSTALL \"/Users/aassefa/Library/CloudStorage/OneDrive-JNJ/JNJ (ATA)/Projects/FlowCytometryand CyTOF data analysis methods/voomCLR\"")

library(LinDA)
library(betareg)
library(nlme)
library(modeest)
library(edgeR)
library(DESeq2)
library(mixtools)
library(DCATS)
library(sccomp)  # added after revision
# install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz", repos = NULL)
# BiocManager::install("ggtree")
# BiocManager::install("phyloseq")
# remotes::install_github("ZRChao/adaANCOM")
# remotes::install_github("liudoubletian/phyloMDA")
# devtools::install_github('zhouhj1994/LinDA')
#library(phyloseq)
# library(phyloMDA)
# library(ape)

calcMode <- function(beta, n){
  suppressMessages(mode <- modeest::mlv(sqrt(n) * beta, 
                                        method = "meanshift", kernel = "gaussian")/sqrt(n))
  mode
}


runMoldel <- function(dat, model = "NBGLM"){
  
  dat.wid <- cbind(id = rownames(dat$Y), total = rowSums(dat$Y),  group=dat$group, dat$Y)
  
  dat.wid.clr <- cbind(id = rownames(dat$Y), group=dat$group, dat$Y) %>%
    pivot_longer(cols=colnames(dat$Y)) %>%
    group_by(id, group) %>%
    mutate(clr = as.numeric(compositions::clr(value+1))) %>% ungroup() %>%
    pivot_wider(id_cols = c(id, group), 
                names_from=name, values_from=clr)
  
  dat.wid.frac <- cbind(id = rownames(dat$Y), group=dat$group, dat$Y) %>%
    pivot_longer(cols=colnames(dat$Y)) %>%
    group_by(id, group) %>%
    mutate(frac = value/sum(value)) %>% ungroup() %>%
    pivot_wider(id_cols = c(id, group), 
                names_from=name, values_from=frac)
  
  dat.wid.logit <- cbind(id = rownames(dat$Y), group=dat$group, dat$Y) %>%
    pivot_longer(cols=colnames(dat$Y)) %>%
    group_by(id, group) %>%
    mutate(frac = value/sum(value),
           lgtFrac= log((frac+0.0001)/(1-(frac+0.0001)))) %>% ungroup() %>%
    pivot_wider(id_cols = c(id, group), 
                names_from=name, values_from=lgtFrac)
  
  ppls.vec <- colnames(dat$Y)
  
  
  if(model == "NBGLM"){ 
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(count = dat.wid %>% pull(all_of(ppl)),
                           N = dat.wid$total,
                           id=as.character(dat.wid$id),
                           x=as.factor(dat.wid$group))
      
      suppressMessages({ 
        res <- try({
          glm.nb(count~x+offset(log(N)), data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))[2, "Estimate"], 
                             test.stat = coef(summary(res))[2, "z value"], 
                             pval = coef(summary(res))[2, "Pr(>|z|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>% data.frame()  
  }
  else if(model == "NBGLM_gmOffset"){ 
    dat.wid.gm <- apply(dat.wid [, -c(1:3)], 1, function(y){
      exp(mean(log(y+1)))
    })
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(count = dat.wid %>% pull(all_of(ppl)),
                           gm = dat.wid.gm,
                           id=as.character(dat.wid$id),
                           x=as.factor(dat.wid$group))
      
      suppressMessages({ 
        res <- try({
          glm.nb(count~x+offset(log(gm)), data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))[2, "Estimate"], 
                             test.stat = coef(summary(res))[2, "z value"], 
                             pval = coef(summary(res))[2, "Pr(>|z|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>% data.frame()  
  }
  else if(model == "NBGLM_TMB"){ 
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(count = dat.wid %>% pull(all_of(ppl)),
                           total = dat.wid$total,
                           id=as.character(dat.wid$id),
                           x=as.factor(dat.wid$group))
      
      suppressMessages({ 
        res <- try({
          glmmTMB(count~x+offset(log(total)), family = nbinom2, data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))$cond[2, "Estimate"], 
                             test.stat = coef(summary(res))$cond[2, "z value"], 
                             pval = coef(summary(res))$cond[2, "Pr(>|z|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "logit_prop"){ 
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(lgtFrac = dat.wid.logit %>% pull(all_of(ppl)),  
                           id=as.character(dat.wid.logit$id),
                           x=as.factor(dat.wid.logit$group))
      
      suppressMessages({ 
        res <- try({
          glm(lgtFrac~x, family = "gaussian", data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))[2, "Estimate"], 
                             test.stat = coef(summary(res))[2, "t value"], 
                             pval = coef(summary(res))[2, "Pr(>|t|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "wlogit_prop"){ 
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(lgtFrac = dat.wid.logit %>% pull(all_of(ppl)),  
                           id=as.character(dat.wid.logit$id),
                           x=as.factor(dat.wid.logit$group))
      
      suppressMessages({ 
        mdl1 <- lm(lgtFrac~x,  data = dat.ppl)
        wts <- 1/(lm(abs(mdl1$residuals) ~ mdl1$fitted.values)$fitted.values^2)
        res <- try({ 
          lm(lgtFrac~x, weights = wts, data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))[2, "Estimate"], 
                             test.stat = coef(summary(res))[2, "t value"], 
                             pval = coef(summary(res))[2, "Pr(>|t|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "quasibinomial"){ 
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(frac = dat.wid.frac %>% pull(all_of(ppl)),  
                           id=as.character(dat.wid.frac$id),
                           x=as.factor(dat.wid.frac$group))
      
      suppressMessages({ 
        res <- try({
          glm(frac~x, family = "quasibinomial",  data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))[2, "Estimate"], 
                             test.stat = coef(summary(res))[2, "t value"], 
                             pval = coef(summary(res))[2, "Pr(>|t|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "betareg_1"){ 
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(frac = dat.wid.frac %>% pull(all_of(ppl)),   
                           id=as.character(dat.wid.frac$id),
                           x=as.factor(dat.wid.frac$group)) %>%
        mutate(frac2 = case_when(frac==0~0.0001, frac==1~1-0.0001, TRUE~frac))
      
      suppressMessages({ 
        res <- try({
          betareg(frac2~x, data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))$mean[2, "Estimate"], 
                             test.stat = coef(summary(res))$mean[2, "z value"], 
                             pval = coef(summary(res))$mean[2, "Pr(>|z|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "betareg_2"){ 
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(frac = dat.wid.frac %>% pull(all_of(ppl)),   
                           id=as.character(dat.wid.frac$id),
                           x=as.factor(dat.wid.frac$group)) %>%
        mutate(frac2 = case_when(frac==0~0.0001, frac==1~1-0.0001, TRUE~frac))
      
      suppressMessages({ 
        res <- try({
          betareg(frac2~x|x, data = dat.ppl)
        }, silent = TRUE)
      })
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))$mean[2, "Estimate"], 
                             test.stat = coef(summary(res))$mean[2, "z value"], 
                             pval = coef(summary(res))$mean[2, "Pr(>|z|)"])  
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
      out.df
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "LMclr"){
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(y = dat.wid.clr %>% pull(all_of(ppl)),
                           id=as.character(dat.wid.clr$id),
                           x=as.factor(dat.wid.clr$group))
       
      res = try(lm(y~x, data = dat.ppl),
                     silent = TRUE)
      
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))[2, "Estimate"], 
                             test.stat = coef(summary(res))[2, "t value"], 
                             pval = coef(summary(res))[2, "Pr(>|t|)"]) 
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "wLMclr"){
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(y = dat.wid.clr %>% pull(all_of(ppl)),
                           id=as.character(dat.wid.clr$id),
                           x=as.factor(dat.wid.clr$group))
      
      mdl1 <- lm(y~x,  data = dat.ppl)
      wts <- 1/(lm(abs(mdl1$residuals) ~ mdl1$fitted.values)$fitted.values^2)
      res = try(lm(y~x,weights =wts,  data = dat.ppl),
                silent = TRUE)
      
      if (class(res)[1] != "try-error"){ 
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res))[2, "Estimate"], 
                             test.stat = coef(summary(res))[2, "t value"], 
                             pval = coef(summary(res))[2, "Pr(>|t|)"]) 
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame()  
  }
  else if(model == "wbcLMclr"){
    df.B <- sapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(y = dat.wid.clr %>% pull(all_of(ppl)),
                           id=as.character(dat.wid.clr$id),
                           x=as.factor(dat.wid.clr$group))
      
      # calculate weights
      mdl1 <- lm(y~x,  data = dat.ppl)
      wts <- 1/(lm(abs(mdl1$residuals) ~ mdl1$fitted.values)$fitted.values^2)
      res = lm(y~x, weights =wts,  data = dat.ppl) 
      res$coefficients[2]
    })
    
    res.df <- lapply(ppls.vec, function(ppl){
      dat.ppl = data.frame(y = dat.wid.clr %>% pull(all_of(ppl)),
                           id=as.character(dat.wid.clr$id),
                           x=as.factor(dat.wid.clr$group))
      
      mdl1 <- lm(y~x,  data = dat.ppl)
      wts <- 1/(lm(abs(mdl1$residuals) ~ mdl1$fitted.values)$fitted.values^2)
      res = try(lm(y~x,weights =wts,  data = dat.ppl),
                silent = TRUE)
      
      if (class(res)[1] != "try-error"){ 
        # bias correction
        res2 <- res
        res2$coefficients[2] <- res2$coefficients[2]-as.numeric(calcMode(beta = df.B, n = length(df.B)))
        out.df  = data.frame(population = ppl,
                             coef.est = coef(summary(res2))[2, "Estimate"], 
                             test.stat = coef(summary(res2))[2, "t value"], 
                             pval = coef(summary(res2))[2, "Pr(>|t|)"]) 
      } else {
        out.df  = data.frame(population = ppl, 
                             coef.est = NA, 
                             test.stat = NA, 
                             pval = NA)
      }
    }) %>% bind_rows() %>%
      mutate(padj= p.adjust(pval, method = "BH"), method=model) %>%data.frame() 
    
  }
  else if(model == "voomCLR"){
    #source("R/voomCLR.R")
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL) 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_bootNonParametric"){
    require(voomCLR)
    #source("R/voomCLR.R")
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL) 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf, bootstrap="nonparametric")
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_NPboot_NBweight"){
    require(voomCLR)
    #source("R/voomCLR.R")
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL,
                         varCalc = "analytical", varDistribution = "NB") 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf, bootstrap="nonparametric")
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_NPboot_PoisWeight"){
    require(voomCLR)
    #source("R/voomCLR.R")
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL,
                         varCalc = "analytical", varDistribution = "poisson") 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf, bootstrap="nonparametric")
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_bootParametric"){
    require(voomCLR)
    #source("R/voomCLR.R")
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL) 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf, bootstrap="parametric",  voomWeights=v$weights)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_Pboot_PoisWeight"){
    require(voomCLR)
    #source("R/voomCLR.R")
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL,
                         varCalc = "analytical", varDistribution = "poisson") 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf, bootstrap="parametric",  voomWeights=v$weights)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_Pboot_NBweight"){
    require(voomCLR)
    #source("R/voomCLR.R")
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL,
                         varCalc = "analytical", varDistribution = "NB") 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf, bootstrap="parametric",  voomWeights=v$weights)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_woB"){ 
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL) 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTable(fit, coef = 2, number = Inf)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_woH"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL) 
    v$weights <- NULL #matrix(1, nrow=nrow(v$weights), ncol=ncol(v$weights))
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_woE"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL)  
    fit <- lmFit(v, design.mat)
    #fite <- eBayes(fit) 
    fit <- applyBiasCorrection(fit)
    fit$t <- fit$coef[,2]/fit$stdev.unscaled[,2]/fit$sigma
    fit$p.value <- 2*pt(-abs(fit$t), df=fit$df.residual)
    fit$padj <- p.adjust(fit$p.value, "BH")
  
    #voomCLR.res <- topTableBC(fit, coef = 2, number = Inf)
    voomCLR.res <- data.frame(logFC = fit$coefficients[,2], t=fit$t, 
                              P.Value = fit$p.value, adj.P.Val=fit$padj)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_PoissonAnalyticalVarCalc"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  
                         lib.size = NULL, varCalc = "analytical", varDistribution = "poisson") 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_NBAnalyticalVarCalc"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  
                         lib.size = NULL, varCalc = "analytical", varDistribution = "NB") 
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTableBC(fit, coef = 2, number = Inf)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_woB_woH"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = 1)
    v$weights <- matrix(1, nrow=nrow(v$weights), ncol=ncol(v$weights))
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTable(fit, coef = 2, number = Inf)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_woB_woE"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = 1)
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit$t <- fit$coef[,2]/fit$stdev.unscaled[,2]/fit$sigma
    fit$p.value <- 2*pt(-abs(fit$t), df=fit$df.residual)
    fit$padj <- p.adjust(fit$p.value, "BH")
    
    voomCLR.res <- data.frame(logFC = fit$coefficients[,2], t=fit$t, 
                              P.Value = fit$p.value, adj.P.Val=fit$padj)  
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_woE_woH"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL)
    v$weights <- matrix(1, nrow=nrow(v$weights), ncol=ncol(v$weights))
    fit <- lmFit(v, design.mat)
    fit <- applyBiasCorrection(fit)
    fit$t <- fit$coef[,2]/fit$stdev.unscaled[,2]/fit$sigma
    fit$p.value <- 2*pt(-abs(fit$t), df=fit$df.residual)
    fit$padj <- p.adjust(fit$p.value, "BH")
    
    voomCLR.res <- data.frame(logFC = fit$coefficients[,2], t=fit$t, 
                              P.Value = fit$p.value, adj.P.Val=fit$padj)  
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_woB_woH_woE"){
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = 1)
    v$weights <- matrix(1, nrow=nrow(v$weights), ncol=ncol(v$weights))
    fit <- lmFit(v, design.mat)
    #fit <- applyBiasCorrection(fit)
    fit$t <- fit$coef[,2]/fit$stdev.unscaled[,2]/fit$sigma
    fit$p.value <- 2*pt(-abs(fit$t), df=fit$df.residual)
    fit$padj <- p.adjust(fit$p.value, "BH")
    
    voomCLR.res <- data.frame(logFC = fit$coefficients[,2], t=fit$t, 
                              P.Value = fit$p.value, adj.P.Val=fit$padj)  
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "voomCLR_randomWeight"){ 
    design.mat <- model.matrix(~group, data = dat.wid) 
    v = voomCLR::voomCLR(counts = t(dat.wid[,-c(1:3)]), design = design.mat,  lib.size = NULL)
    rw <- sample(as.vector(v$weights), length(as.vector(v$weights)))  #shuffle weights
    v$weights <- matrix(rw, nrow=nrow(v$weights), ncol=ncol(v$weights))
    fit <- lmFit(v, design.mat)
    fit <- applyBiasCorrection(fit)
    fit <- eBayes(fit) 
    voomCLR.res <- topTable(fit, coef = 2, number = Inf)
    
    res.df <- lapply(ppls.vec, function(ppl){ 
      data.frame(population = ppl, 
                 coef.est = voomCLR.res[ppl, "logFC"],
                 test.stat = voomCLR.res[ppl, "t"],  
                 pval = voomCLR.res[ppl, "P.Value"],
                 padj = voomCLR.res[ppl, "adj.P.Val"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "linDA"){ 
    linda.res = linda(otu.tab = as.matrix(t(dat.wid[,-c(1:3)])), 
                      meta = dat.wid[, c(1:3)],
                      formula = '~ group', type = "count", 
                      adaptive=FALSE, imputation = FALSE)
    res.df <- lapply(ppls.vec, function(ppl){
      data.frame(population = ppl, 
                 coef.est = linda.res$output$group[ppl, "log2FoldChange"],
                 test.stat = linda.res$output$group[ppl, "stat"],
                 pval = linda.res$output$group[ppl, "pvalue"],
                 padj = linda.res$output$group[ppl, "padj"]) 
    }) %>% bind_rows() %>% mutate(method=model) %>% data.frame()  
  }
  else if(model == "edgeR"){ 
    require(edgeR)
    counts <- t(dat.wid[, -c(1:3)])
    group  <- factor(dat.wid$group)
    y      <- DGEList(counts=counts, group=group)
    y      <- calcNormFactors(y)
    design <- model.matrix(~group)
    y      <- estimateDisp(y,design)
    fit    <- glmQLFit(y,design) # quasi-likelihood F-tests
    qlf    <- glmQLFTest(fit,coef=2)
    topTab <- topTags(qlf, sort.by = "none", n = Inf)
    res.df <- data.frame(population = rownames(topTab$table), 
             coef.est = topTab$table$logFC,
             test.stat = topTab$table$`F`,
             pval = topTab$table$PValue,
             padj = topTab$table$FDR, 
             method=model) 
    
  }
  else if(model == "edgeR_woLSnorm"){ 
    counts <- t(dat.wid[, -c(1:3)])
    group  <- factor(dat.wid$group)
    y      <- DGEList(counts=counts,group=group)
    #y      <- calcNormFactors(y)
    design <- model.matrix(~group)
    y      <- estimateDisp(y,design)
    fit    <- glmQLFit(y,design) # quasi-likelihood F-tests
    qlf    <- glmQLFTest(fit,coef=2)
    topTab <- topTags(qlf, n = Inf)
    res.df <- data.frame(population = rownames(topTab$table), 
                         coef.est = topTab$table$logFC,
                         test.stat = topTab$table$`F`,
                         pval = topTab$table$PValue,
                         padj = topTab$table$FDR, 
                         method=model) 
    
  }
  else if(model == "DESeq2"){ 
    counts  <- t(dat.wid[, -c(1:3)])
    group   <- factor(dat.wid$group) 
    coldata <- data.frame(group=group, row.names = colnames(counts))
    dds     <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design = ~ group)
    dds <- DESeq(dds, quiet=TRUE, fitType='local')
    res <- results(dds)
    res.df <- data.frame(population = rownames(res), 
                         coef.est = res$log2FoldChange,
                         test.stat = res$stat,
                         pval = res$pvalue,
                         padj = res$padj, 
                         method=model) 
    
  }
  else if(model == "DESeq2_woLSnorm"){ 
    counts  <- t(dat.wid[, -c(1:3)])
    group   <- factor(dat.wid$group) 
    coldata <- data.frame(group=group, row.names = colnames(counts))
    dds     <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = coldata,
                                      design = ~ group) 
    normalizationFactors(dds) <- matrix(1, nrow(counts), ncol(counts))
    dds <- DESeq(dds, quiet=TRUE, fitType='local') 
    res <- results(dds)
    res.df <- data.frame(population = rownames(res), 
                         coef.est = res$log2FoldChange,
                         test.stat = res$stat,
                         pval = res$pvalue,
                         padj = res$padj, 
                         method=model) 
    
  }
  else if(model == "limmaVoom"){
    require(limma)
    counts <- t(dat.wid[, -c(1:3)])
    group  <- factor(dat.wid$group)
    y      <- DGEList(counts=counts, group=group)
    y      <- calcNormFactors(y)
    design <- model.matrix(~group)
    y   <- voom(y, design, plot = FALSE)
    fit <- lmFit(y, design)
    fit <- eBayes(fit)
    topTab <- topTable(fit, coef = 2, sort.by = "none", n = Inf)
    res.df <- data.frame(population = rownames(topTab), 
                         coef.est = topTab$logFC,
                         test.stat = topTab$`t`,
                         pval = topTab$P.Value,
                         padj = topTab$adj.P.Val, 
                         method=model) 
  }
  else if(model == "limmaVoom_woLSnorm"){
    require(limma)
    counts <- t(dat.wid[, -c(1:3)])
    group  <- factor(dat.wid$group)
    y      <- DGEList(counts=counts, group=group)
    #y      <- calcNormFactors(y)
    design <- model.matrix(~group)
    y   <- voom(y, design, plot = FALSE)
    fit <- lmFit(y, design)
    fit <- eBayes(fit)
    topTab <- topTable(fit, coef = 2, sort.by = "none", n = Inf)
    res.df <- data.frame(population = rownames(topTab), 
                         coef.est = topTab$logFC,
                         test.stat = topTab$`t`,
                         pval = topTab$P.Value,
                         padj = topTab$adj.P.Val, 
                         method=model) 
  }
  else if(model == "propeller_logit"){ 
    source("R/fitPropeller.R")
    counts  <- t(dat.wid[, -c(1:3)])
    group   <- factor(dat.wid$group) 
    dat.prop <- dat.wid %>%
      pivot_longer(cols = -c(1:3), names_to = "population", values_to = "count") %>%
      dplyr::rename("sample"="id", "x"="group") %>%
      mutate(t=1)
    dat.prop.list <- list(count=dat.prop)
    
    res=fitPropeller(simDat = dat.prop.list, transform = "logit")
    
    res.df <- data.frame(population = rownames(res), 
                         coef.est = res$logFC,
                         test.stat = res$t,
                         pval = res$P.Value,
                         padj = res$adj.P.Val, 
                         method=model) 
    
  }
  else if(model == "propeller_asin"){ 
    source("R/fitPropeller.R")
    counts  <- t(dat.wid[, -c(1:3)])
    group   <- factor(dat.wid$group) 
    dat.prop <- dat.wid %>%
      pivot_longer(cols = -c(1:3), names_to = "population", values_to = "count") %>%
      dplyr::rename("sample"="id", "x"="group") %>%
      mutate(t=1)
    dat.prop.list <- list(count=dat.prop)
    
    res=fitPropeller(simDat = dat.prop.list, transform = "asin")
    
    res.df <- data.frame(population = rownames(res), 
                         coef.est = res$logFC,
                         test.stat = res$t,
                         pval = res$P.Value,
                         padj = res$adj.P.Val, 
                         method=model) 
    
  }
  else if(model == "dcats"){ 
    require(DCATS)
    counts  <- dat.wid[, -c(1:3)]
    group   <- data.frame(x=factor(dat.wid$group))
     
    res=dcats_GLM(count_mat = counts, design_mat = group)
    
    res.df <- data.frame(population = rownames(res$ceoffs), 
                         coef.est = res$ceoffs[,"x"],
                         test.stat = res$LR_vals[,"x"],
                         pval = res$LRT_pvals[,"x"],
                         padj = res$fdr[,"x"], 
                         method=model) 
    
  }else if(model=="SCCOMP"){
    require(sccomp)
    sccomp.ounts_obj = dat.wid %>%
      pivot_longer(cols = c(-id, -total, -group), names_to = "ppls", values_to = "count") %>%
      mutate(count=as.integer(count))
    
    sccomp_result = 
      sccomp.ounts_obj |>
      sccomp_estimate( 
        formula_composition = ~ group, 
        .sample = id,
        .cell_group = ppls,
        .abundance = count, 
        bimodal_mean_variability_association  = TRUE,
        cores = 1, verbose = FALSE
      ) |> 
      sccomp_test() |>
      subset(parameter != "(Intercept)")
    
    res.df <- data.frame(population = sccomp_result$ppls, 
                         coef.est = sccomp_result$c_effect,
                         test.stat = NA,
                         pval = NA,
                         padj = sccomp_result$c_FDR, 
                         method="sccomp") 
  }
  else{
    error("Worng method")
  }
  
  res.df
}
