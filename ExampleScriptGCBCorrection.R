
load("DataGRPCorrection.RData")
library(nimble)
library(nimbleEcology)




MA1_P1<-nimbleCode({
  for (e in 1:2){ ###detection 'fixed' effects
    mu_a[e]~dnorm(0, sd=1) 
    sig_a[e]~T(dnorm(0, sd = 1), 0, )
    for (i in 1:nspec){
      a[e,i]~dnorm(mu_a[e], sd=sig_a[e])
    }
  }
  
  sig_det~T(dnorm(0, sd = 1), 0, ) 
  
  for (i in 1:30){
    mu_b[i]~dnorm(0, sd=0.75)
    sig_b[i]~T(dnorm(0, sd = 0.75), 0, )
    for (s in 1:nspec){
      b[i, s]~dnorm(mu_b[i], sd=sig_b[i])
    }
  }
  
  for (e in 1:2){
    theta[1, e]~dgamma(5, 1)
    theta[2, e]~dgamma(10, 1)
    theta[3, e]~dgamma(10, 1)
    
    for (k in 1:3){
      tau[k, e]<-prod(theta[1:k, e])
      for (i in 1:nspec){
        eta[i, k, e]~dt(0, tau=tau[k, e], 3)
      }  
    }  
  }  
  
  for (s in 1:nsites){
    for (i in 1:nspec){
      psi[i, s] <- iprobit(inprod(b[1:8, i], X1[s,1:8])+inprod(eta[i, 1:3, 1], lambda[s, 1:3, 1]))
      phi[i, s] <- iprobit(inprod(b[9:19, i], X2[s,1:11])+inprod(eta[i, 1:3, 2], lambda[s, 1:3, 2]))
      gamma[i, s] <- iprobit(inprod(b[20:30, i], X2[s,1:11])+inprod(eta[i, 1:3, 2], lambda[s, 1:3, 2]))
      y[s,1:2,1:5,i]~dDynOcc_ssm(init=psi[i, s], probPersist=phi[i, s],
                                 probColonize = gamma[i, s], p=p[s,1:2, 1:5,i],
                                 start=Starts2[s,1:2], end=Ends[s,1:2]) ###note, new ends.some=0
    }
  }
  
  for (s in 1:nsites){
    for (e in 1:2){
      lambda[s, 1, e]~dnorm(0, 1)
      lambda[s, 2, e]~dnorm(0, 1)
      lambda[s, 3, e]~dnorm(0, 1)
      for (j in 1:5){
        eps[s, e, j]~dnorm(0, sd=sig_det)
        for (i in 1:nspec){
          logit(p[s, e, j, i])<-eps[s, e, j]+a[e, i] 
        }
      }
    }
  }
  
})


Constants<-list(Starts=matrix(ifelse(N_Surveys>0, 1, 0),  320, 2), Ends=ifelse(N_Surveys>5, 5, N_Surveys),
                nspec=185, nsites=320,
                X1=cbind(rep(1, 320),
                         as.numeric(scale(Site_Data$MeanHistoricPPT)),
                         as.numeric(scale(Site_Data$MeanHistoricTMin)),
                         as.numeric(scale(Site_Data$HistoricTMinFPC2)),
                         as.numeric(scale(Site_Data$Prop_Ag_Hist)),
                         as.numeric(scale(Site_Data$PropUrbanHist)),
                         as.numeric(scale(Site_Data$MeanHistoricTMin))^2,
                         as.numeric(scale(Site_Data$HistoricTMinFPC2))^2),
                X2=cbind(rep(1, 320),
                         as.numeric(scale(Site_Data$MeanModernPPT)),
                         as.numeric(scale(Site_Data$MeanModernTMin)),
                         as.numeric(scale(Site_Data$ModernTMinFPC2)),
                         as.numeric(scale(Site_Data$PropAgMod)),
                         as.numeric(scale(Site_Data$PropUrbanMod)),
                         as.numeric(scale(Site_Data$MeanModernTMin))^2,
                         as.numeric(scale(Site_Data$ModernTMinFPC2))^2,
                         as.numeric(scale(Site_Data$MeanModernPPT-Site_Data2$MeanHistoricPPT)),
                         as.numeric(scale(Site_Data$MeanModernTMin-Site_Data2$MeanHistoricTMin)),
                         as.numeric(scale(Site_Data$ModernTMinFPC2-Site_Data2$HistoricTMinFPC2))))

Data<-list(y=y[, , 1:5, ])


###This shouldn't be neccessary, but will avoid getting a warning.
yIn<-array(NA, dim=dim(y))
yIn[is.na(y)]<-0


###Here, just setting very precise starting values that guarantee the model will start sampling.
Inits <- list(y=yIn[,,1:5,], mu_b=rep(0, 30), sig_b=rep(.25, 30), 
              a=matrix(rnorm(370, 0, .01), 2, 185), mu_a=rep(0, 2), sig_a=rep(.25, 2),
              b=matrix(rnorm(5550, 0, .01), 30, 185), eps=array(rnorm(3200, 0, .02), dim=c(320, 2, 5)),
              eta=array(0, dim=c(185, 3, 2)), lambda=array(0, dim=c(320, 3, 2)),
              theta=matrix(rep(c(5, 10, 10), 2), 3, 2))



GRP_Mod <- nimbleModel(code = MA1_P1, name = 'A1_P1', constants = Constants,
                         data=Data, inits=Inits, calculate=FALSE)

GRP_ModConf <- configureMCMC(GRP_Mod,
                               monitors = c("mu_a", "mu_b", "sig_a", "sig_b", "a", "b", "sig_det", "eps", "lambda", "eta"), UseConjugacy=FALSE) ###z, cp, etc.

GRP_ModConf$removeSamplers(c("sig_b[5]", "sig_b[13]", "sig_b[24]"))
GRP_ModConf$addSampler(target ='sig_b[5]', type = 'RW', control=list(log=TRUE))
GRP_ModConf$addSampler(target ='sig_b[13]', type='RW', control=list(log=TRUE))
GRP_ModConf$addSampler(target = 'sig_b[24]', type= 'RW', control=list(log=TRUE))

Rmcmc<-buildMCMC(GRP_ModConf)
compMCMC <- compileNimble(Rmcmc, GRP_Mod)

###Note, in practice the model was run on an HPC with a limited run window, and so we ran short chains
###on independent nodes, took the last values and used them as starting values for the next set of chains,
###and so forth. I roughly estimate this would work in a single go (although it would take a while and 
###eat up a bunch of RAM).
samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=300000, nburnin=200000, thin=50, 
               nchains=3)







###Post processing
###Assumes the inputs/constants/samples are still in the workspace.
library(loo)

N_Surveys[268,]<-c(4, 4)

CSamps<-rbind(samps[[1]], samps[[2]], samps[[3]])


grab <- function(x, y) {x[,grep(y,colnames(x))]}
theta1<-as.matrix(grab(CSamps, "^theta\\["))

eps1<-as.matrix(grab(CSamps, "^eps\\["))
eps<-array(NA, dim=c(320, 2, 5, 6000))
lambda1<-as.matrix(grab(CSamps, "^lambda\\["))
lambda<-array(NA, dim=c(320, 3, 2, 6000))
eta1<-as.matrix(grab(CSamps, "^eta\\["))
eta<-array(NA, dim=c(185,3, 2, 6000))
for (i in 1:6000){
  eps[,,,i]<-eps1[i, ]
  lambda[,,,i]<-lambda1[i, ]
  eta[,,,i]<-eta1[i, ]
}

#cor1<-array(NA, dim=c(185, 185, 9000))
#cor2<-array(NA, dim=c(185, 185, 9000))

#for (i in 1:9000){
#  cor1[,,i]<-eta[,,1,i]%*%t(eta[,,1,i])
#  cor2[,,i]<-eta[,,2,i]%*%t(eta[,,2,i])
#}




a<-array(t(as.matrix(grab(CSamps, "^a\\["))), dim=c(2,185,6000))
b<-array(t(as.matrix(grab(CSamps, "^b\\["))), dim=c(30,185,6000))




###Derive the LOO-PSIS. Presumably the Constants are still in the workspace.
loglike<-matrix(NA, 185*320, 6000)


for (l in 1:185){
  p<-array(NA, dim=c(320, 2, 5, 6000))
  
  for (i in 1:320){
    for (m in 1:6000){
      psi<-pnorm(Constants$X1[i,] %*% b[1:8,l,m]+eta[l, ,1,m] %*% lambda[i, ,1,m])
      phi<-pnorm(Constants$X2[i, ] %*% b[9:19,l,m]+eta[l, ,2,m] %*% lambda[i, ,2,m])
      gamma<-pnorm(Constants$X2[i, ] %*% b[20:30,l,m]+eta[l, ,2,m] %*% lambda[i, ,2,m])
      
      for (j in 1:2){
        for (k in 1:5){
          p[i, j, k, m]<-plogis(a[j, l, m]+eps[i, j, k, m])
        }
      }
      loglike[(i-1)*185+l, m]<-dDynOcc_ssm(y[i, 1:2, 1:5, l], init=psi, probPersist=phi,
                                           probColonize = gamma, p=p[i ,1:2, 1:5,m],
                                           start=Constants$Starts[i,1:2], end=Constants$Ends[i,1:2], log=1)
      
      
      ###run this if you want a PPC
      #ysim<-rDynOcc_ssm(1, .init=psi, probPersist=phi,
      #                  probColonize = gamma, p=p[i ,1:2, 1:5,m],
      #                 start=Constants$Starts[i,1:2], end=Constants$Ends[i,1:2])
      # 
      # remember to create a loglike sim array around L179 if you want this
      #loglikesim[(i-1)*185+l, m]<-dDynOcc_ssm(ysim, init=psi, probPersist=phi,
      # probColonize = gamma, p=p[i ,1:2, 1:5,m],
      # start=Constants$Starts[i,1:2], end=Constants$Ends[i,1:2], log=1)
    }
  }
}
loomod<-loo::loo(t(loglike))

dev<-loglike*-2
devsim<-loglikesim*-2
BPV<-mean(colSums(dev)>colSums(devsim))
