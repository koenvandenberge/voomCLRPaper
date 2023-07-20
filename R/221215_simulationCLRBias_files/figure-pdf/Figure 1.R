# load Lups data
lupsPoplCount <- readRDS("Data/lupsPoplCountV2_LongIndividualBatchID.rds") %>%
  mutate(group=case_when(group=="Healthy"~"Healthy", TRUE~"SLE"))
# head(lupsPoplCount)
 
lupsPoplCountCh4 <- lupsPoplCount %>%
  subset(Processing_Cohort == "4.0")

lupsPoplCountCh4_2 <- lupsPoplCountCh4 %>%  
  pivot_wider(id_cols = c(patient, group), 
              names_from = celltype, values_from = nCells) %>%
  replace(is.na(.), 0) %>% 
  pivot_longer(cols = unique(lupsPoplCountCh4$celltype), names_to = "celltype", values_to = "nCells")


lupsPoplCountCh4_2_wide <- lupsPoplCountCh4_2 %>%
  dplyr::select(patient, group, celltype, nCells) %>%
  pivot_wider(id_cols = c(patient, group), names_from = celltype, values_from = nCells) %>%
  data.frame()
rownames(lupsPoplCountCh4_2_wide) <- lupsPoplCountCh4_2_wide$patient

Xhealthy <- lupsPoplCountCh4_2_wide[lupsPoplCountCh4_2_wide$group == "Healthy", -c(1:2)]


# run simulation 
source("R/simCP_v3.R")
source("R/runModels.R")
source("R/DM_simCP.R")
 
models <- c("voomCLR_woH", "voomCLR_woB_woH")

seed = 4513#64868
var_scale = 1
nn=20

simData <- DM_simCP(sample_size = c(nn, nn), 
                    K = ncol(Xhealthy), frac.daPop = 0.01, abs.treat.effect.range = c(2, 5), 
                    parent.size = rowSums(Xhealthy)[sample(nrow(Xhealthy), nn*2)], 
                    sd.intercept = 1,
                    scale.dirichlet = var_scale,  seed = seed)

# apply models 
resAll <- lapply(models, function(mdl){
  simInput <- list(Y=simData$X, group=simData$group) 
  res.sim <- runMoldel(dat = simInput, model = mdl) 
  res.sim$trueDA <- res.sim$population %in% simData$DApops[,1]
  res.sim
}) %>% bind_rows() %>%
  mutate( seed=seed) %>%
  mutate(biasCorrected=case_when(method=="voomCLR_woH"~TRUE, TRUE~FALSE)) 

df.trueES = data.frame(population = colnames(simData$X), true.ES=simData$true.ES)
resAll <- merge(resAll, df.trueES, by="population")

####. Figure 1 a
ppa=resAll %>%
  dplyr::select(population, coef.est, true.ES, biasCorrected) %>% 
  ggplot(aes(x=as.numeric(gsub("ppl_", "", population)), y=coef.est-true.ES, 
             color=biasCorrected))+
  scale_x_continuous(breaks = 1:11)+
  geom_point(size=2)+ 
  theme_classic() +
  geom_hline(data = resAll %>%
               group_by(biasCorrected) %>%
               summarise(m=mean(coef.est-true.ES)),
             aes(yintercept=m, color=biasCorrected),
             lty=2)+
  geom_hline(yintercept = 0)+
  labs(x="cell population", y="Difference between estimated vs true LFC",
       color="Bias corrected", title = "a")+
  theme(legend.position = "top")
  
####. Figure 1 b & c 
seed = 83270 #round(runif(1, 1, 1e5))
seed
var_scale = 1.5
nn=20

simData <- DM_simCP(sample_size = c(nn, nn), 
                    K = ncol(Xhealthy), frac.daPop = 0, abs.treat.effect.range = c(0.1, 0.5), 
                    parent.size = rep(10000, 2*nn), mean.intercept = 20, sd.intercept = 1,
                    scale.dirichlet = var_scale,  seed = seed)


MNvar <- function(K, frac, N){
  ((1-1/K)^2)*((1-frac)/(N*frac))
}
simdata.df <- simData$X %>%  
  rownames_to_column(var = "id") %>%
  #dplyr::select(-simData$DApops[,1]) %>%
  pivot_longer(cols = -id, names_to = "population", values_to = "y") %>%
  group_by(id) %>%
  mutate(CLR=as.numeric(clr(y)), N=sum(y)) %>%
  ungroup() %>% 
  mutate(VarY = mean(y), 
         varCLR2 = MNvar(K=11, N = N, frac=y/N)) %>%
  group_by(population) %>% 
  mutate(meanCount=mean(y), varCount=var(y), 
            meanCLR=mean(CLR), varCLR=var(CLR),
            n=n(),
            varCount2 = mean(VarY),
            varCLR2=1/meanCount) %>%
  ungroup()

ppb = ggplot(simdata.df, aes(x=meanCount, y=varCount))+
  geom_point(size=2)+
  geom_line(aes(y=meanCount),lwd=1, lty=2,  col="blue")+
  #stat_smooth(method = "gam", se=FALSE)+
  theme_classic() +
  labs(x="Mean of counts", y="Variance of counts", title = "b")

ppc = ggplot(simdata.df, aes(x=meanCLR, y=varCLR))+
  geom_point(size=2)+
  geom_line(aes(y=varCLR2), lwd=1, lty=2, col="blue")+
  #stat_smooth(method = "gam", se=FALSE)+
  theme_classic() +
  labs(x="Mean of CLR", y="Variance of CLR", title = "c")
  
#gridExtra::grid.arrange(grobs=list(ppb, ppc), ncol=2)



####. Figure 1 d & e 

lupus.df <- lupsPoplCountCh4_2 %>%  
  mutate(nCells=nCells+1) %>%
  subset(group=="Healthy") %>%
  group_by(patient) %>%
  mutate(CLR=as.numeric(clr(nCells)), N=sum(nCells)) %>%
  ungroup() %>%
  # mutate(VarY = nCells*(1-nCells/N), 
  #        varCLR2 = MNvar(K=11, N = N, frac=nCells/N)) %>%
  group_by(celltype) %>% 
  summarise(meanCount=mean(nCells), varCount=var(nCells), 
            meanCLR=mean(CLR), varCLR=var(CLR)) %>%
  ungroup()

x=lupus.df$meanCLR
names(x) = lupus.df$celltype

library(glmmTMB)
estVec <- sapply(names(x), function(ppl){ 
  df = lupsPoplCountCh4_2[lupsPoplCountCh4_2$celltype==ppl,]
  mnb <- suppressWarnings(glmmTMB::glmmTMB(df$nCells ~ 1, 
                                           family=nbinom2(link = "log")))
  c(disp=sigma(mnb), meanCount=exp(mnb$fit$par[1]))
})

K=11
y <- ((1-1/K)^2)*(1/estVec[2,] + 1/(estVec[1,])^2)
# library(mgcv)
# mdlCLR <- gam(y~-1+s(x))
# xgrid=seq(-3.5, 3.5, 0.1)
# predVarCLR <- predict(mdlCLR, newdata=data.frame(x=xgrid))


ppd = ggplot(lupus.df, aes(x=meanCount, y=varCount))+ 
  geom_point(size=2)+
  #geom_line(aes(y=varCount2), col="blue")+
  stat_smooth(method = "loess", span = 2, lwd=1, lty=2, se=FALSE)+
  theme_classic() +
  labs(x="Mean of counts", y="Variance of counts", title = "d")

ppe = ggplot(lupus.df, aes(x=meanCLR, y=varCLR))+
  geom_point(size=2)+
  stat_smooth(data=data.frame(x, y),
              aes(x, y), method = "loess", formula = y~-1+I(x),  method.args = list(degree=1, span=2),
              se = FALSE, lwd=1, lty=2, col="blue")+
  # stat_smooth(data=data.frame(meanCLR=x, varCLR=predVarCLR),
  #             method = "loess", span = 2, se=FALSE)+
  theme_classic() +
  labs(x="Mean of CLR", y="Variance of CLR", title = "e")

#gridExtra::grid.arrange(grobs=list(ppd, ppe), ncol=2)

png("results/Figure1.png", width = 26, height = 14, units = "cm", res = 150)
gridExtra::grid.arrange(grobs=list(ppa, ppb, ppc, ppd, ppe), 
                        layout_matrix=matrix(c(1, 2, 3, 1,4, 5), ncol = 3, nrow = 2, 
                                             byrow = TRUE), 
                        widths = c(0.4, 0.3, 0.3))

dev.off()

