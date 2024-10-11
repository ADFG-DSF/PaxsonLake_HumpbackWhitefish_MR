library(tidyverse)
library(dsftools)
library(jagsUI)
library(jagshelper)



# reading data
spawn_sample <- read_csv("FDS_2024/flat_data/spawn_sample.csv")%>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(Date = as.Date(Date, format="%m/%d/%Y"))

# defining breaks for age and length
agebreaks <- c(5, 10, 15 ,20, 30, 40)
lengthbreaks <- c(250, 350, 400, 500)
spawn_sample$agecut <- cut(spawn_sample$Age, breaks=agebreaks, right=FALSE)
spawn_sample$lengthcut <- cut(spawn_sample$Length, breaks=lengthbreaks, right=FALSE)

# tabulating length and age (binned)
lengthage_tab <- with(spawn_sample, 
     table(agecut,
           lengthcut))
lengthage_tab
mosaicplot(lengthage_tab, 
           xlab="Age", ylab="Length (mm FL)", main="",
           col=rev(grey.colors(ncol(lengthage_tab))))

# tabulating mean length & proportions per age bin
mnlength_tab <- with(spawn_sample, 
                      ASL_table(age = agecut,
                                length = Length))
mnlength_tab
with(spawn_sample,
     boxplot(Length ~ agecut, 
             at=(agebreaks[-1]+agebreaks[-length(agebreaks)])/2/5,
             xlab="Age", ylab="Length (mm FL)", main=""))

with(spawn_sample,
     plot(Length ~ Age,
          xlab="Age", ylab="Length (mm FL)", main=""))




## trying multi-model   - LOGNORMAL ERROR
cat('model {
  for(i in 1:n) {
    for(j in 1:4) {
      L[i,j] ~ dlnorm(logmu[i,j], tau[j])
      ypp[i,j] ~ dlnorm(logmu[i,j], tau[j])
    }
    logmu[i,1] <- log(L_inf[1]*(1-exp(-k[1]*(a[i]-t0[1]))))      # VBM - VB
    logmu[i,2] <- log(L_inf[2]*(1-exp(-k[2]*(a[i]-t0[2])))^P[2])    # GGM - generalized VB
    logmu[i,3] <- log(L_inf[3]/(1+exp(-k[3]*(a[i]-t0[3]))))      # LGM - logistic growth model
    logmu[i,4] <- log(L_inf[4]*exp((-1/k[4])*exp(-k[4]*(a[i]-t0[4])))) # PGM - Gompertz
  }

  for(j in 1:4) {
    tau[j] <- pow(sig[j], -2)
    sig[j] ~ dunif(0, 10)#00)
    L_inf[j] ~ dnorm(400, .0001)
    k[j] ~ dnorm(0,.1)T(0,)
    t0[j] ~ dnorm(0,.01)T(,6)#<-0#
    P[j] ~ dlnorm(0,1)
  }
  
  for(ifit in 1:nfit) {
    logmufit[ifit,1] <- log(L_inf[1]*(1-exp(-k[1]*(agefit[ifit]-t0[1]))))      # VBM - VB
    logmufit[ifit,2] <- log(L_inf[2]*(1-exp(-k[2]*(agefit[ifit]-t0[2])))^P[2])    # GGM - generalized VB
    logmufit[ifit,3] <- log(L_inf[3]/(1+exp(-k[3]*(agefit[ifit]-t0[3]))))      # LGM - logistic growth model
    logmufit[ifit,4] <- log(L_inf[4]*exp((-1/k[4])*exp(-k[4]*(agefit[ifit]-t0[4])))) # PGM - Gompertz
    for(j in 1:4) {
      ppfit[ifit,j] ~ dlnorm(logmufit[ifit,j], tau[j])
      mufit[ifit,j] <- exp(logmufit[ifit,j])
    }
  }
  
}', file="lvb_jags2_multi_t0_free_lnorm")

cat('model {
  for(i in 1:n) {
    for(j in 1:4) {
      L[i,j] ~ dlnorm(logmu[i,j], tau[j])
      ypp[i,j] ~ dlnorm(logmu[i,j], tau[j])
    }
    logmu[i,1] <- log(L_inf[1]*(1-exp(-k[1]*(a[i]-t0[1]))))      # VBM - VB
    logmu[i,2] <- log(L_inf[2]*(1-exp(-k[2]*(a[i]-t0[2])))^P[2])    # GGM - generalized VB
    logmu[i,3] <- log(L_inf[3]/(1+exp(-k[3]*(a[i]-t0[3]))))      # LGM - logistic growth model
    logmu[i,4] <- log(L_inf[4]*exp((-1/k[4])*exp(-k[4]*(a[i]-t0[4])))) # PGM - Gompertz
  }

  for(j in 1:4) {
    tau[j] <- pow(sig[j], -2)
    sig[j] ~ dunif(0, 10)#00)
    L_inf[j] ~ dnorm(400, .0001)
    k[j] ~ dnorm(0,.1)T(0,)
    t0[j] <-0 #~ dnorm(0,.01)T(,6)#
    P[j] ~ dlnorm(0,1)
  }
  
  for(ifit in 1:nfit) {
    logmufit[ifit,1] <- log(L_inf[1]*(1-exp(-k[1]*(agefit[ifit]-t0[1]))))      # VBM - VB
    logmufit[ifit,2] <- log(L_inf[2]*(1-exp(-k[2]*(agefit[ifit]-t0[2])))^P[2])    # GGM - generalized VB
    logmufit[ifit,3] <- log(L_inf[3]/(1+exp(-k[3]*(agefit[ifit]-t0[3]))))      # LGM - logistic growth model
    logmufit[ifit,4] <- log(L_inf[4]*exp((-1/k[4])*exp(-k[4]*(agefit[ifit]-t0[4])))) # PGM - Gompertz
    for(j in 1:4) {
      ppfit[ifit,j] ~ dlnorm(logmufit[ifit,j], tau[j])
      mufit[ifit,j] <- exp(logmufit[ifit,j])
    }
  }

}', file="lvb_jags2_multi_t0_0_lnorm")




modelnames <- c("VB", "Generalized VB", "Logistic", "Gompertz")


# bundle data to pass into JAGS
nfit <- 100
lvbdata <- subset(spawn_sample, !is.na(Age) & !is.na(Length))

lengthfit <- seq(min(lvbdata$Length), max(lvbdata$Length), length.out=nfit)
lvb_data2_multi <- list(L=matrix(lvbdata$Length,ncol=4,nrow=nrow(lvbdata)),
                        a=lvbdata$`Age`,
                        n=nrow(lvbdata),
                        agefit=seq(5,40,length.out=nfit),
                        nfit=nfit,
                        ppp=rep(1,4))


## check if posts/ files exist, otherwise run

run_anyway <- FALSE  # set this to TRUE to run the models anyway


allfiles <- list.files(recursive=TRUE)

if(run_anyway | !("FDS_2024/posts/lvb_jags_out2_multi_t0_free_lnorm.Rdata" %in% allfiles)) {
  
  # JAGS controls - t0 free
  niter <- 1000*1000 # 49 min at 1000k
  ncores <- min(parallel::detectCores()-1, 10)
  
  {
    tstart <- Sys.time()
    print(tstart)
    lvb_jags_out2_multi_t0_free_lnorm <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free_lnorm", data=lvb_data2_multi,
                                                      parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P","mufit","ppfit","logmu","whichmodel","pp"),
                                                      n.chains=ncores, parallel=T, n.iter=niter,
                                                      n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
  nbyname(lvb_jags_out2_multi_t0_free_lnorm)
  plotRhats(lvb_jags_out2_multi_t0_free_lnorm)
  par(mfrow=c(2,2))
  traceworstRhat(lvb_jags_out2_multi_t0_free_lnorm)
  
  # save(lvb_jags_out2_multi_t0_free_lnorm, file="FDS_2024/posts/lvb_jags_out2_multi_t0_free_lnorm.Rdata")
} else {
  load(file="FDS_2024/posts/lvb_jags_out2_multi_t0_free_lnorm.Rdata")
}




if(run_anyway | !("FDS_2024/posts/lvb_jags_out2_multi_t0_0_lnorm.Rdata" %in% allfiles)) {
  
  # JAGS controls - t0 set to zero
  niter <- 500*1000 # 18 min at 500k
  ncores <- min(parallel::detectCores()-1, 10)
  
  {
    tstart <- Sys.time()
    print(tstart)
    lvb_jags_out2_multi_t0_0_lnorm <- jagsUI::jags(model.file="lvb_jags2_multi_t0_0_lnorm", data=lvb_data2_multi,
                                                   parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P","mufit","ppfit","logmu","whichmodel","pp"),
                                                   n.chains=ncores, parallel=T, n.iter=niter,
                                                   n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
  nbyname(lvb_jags_out2_multi_t0_0_lnorm)
  plotRhats(lvb_jags_out2_multi_t0_0_lnorm)
  par(mfrow=c(2,2))
  traceworstRhat(lvb_jags_out2_multi_t0_0_lnorm)
  
  # save(lvb_jags_out2_multi_t0_0_lnorm, file="FDS_2024/posts/lvb_jags_out2_multi_t0_0_lnorm.Rdata")
} else {
  load(file="FDS_2024/posts/lvb_jags_out2_multi_t0_0_lnorm.Rdata")
}



par(mfrow=c(2,2))
for(j in 1:4) {   # t0 free
  plot(lvbdata$Age, lvbdata$Length, main=paste(modelnames[j], "- t0 free"))
  envelope(lvb_jags_out2_multi_t0_free_lnorm$sims.list$mufit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}
for(j in 1:4) {   # t0 set to zero
  plot(lvbdata$Age, lvbdata$Length, main=paste(modelnames[j], "- t0 set to zero"))
  envelope(lvb_jags_out2_multi_t0_0_lnorm$sims.list$mufit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}

for(j in 1:4) {   # t0 free
  plot(lvbdata$Age, lvbdata$Length, main=paste(modelnames[j], "- t0 free"))
  envelope(lvb_jags_out2_multi_t0_free_lnorm$sims.list$ppfit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}
for(j in 1:4) {   # t0 set to zero
  plot(lvbdata$Age, lvbdata$Length, main=paste(modelnames[j], "- t0 set to zero"))
  envelope(lvb_jags_out2_multi_t0_0_lnorm$sims.list$ppfit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}

for(j in 1:4) {   # t0 free
  qq_postpred(ypp = lvb_jags_out2_multi_t0_free_lnorm$sims.list$ypp[,,j], y=lvbdata$Length, 
              main=paste(modelnames[j], "- t0 free"))
}
for(j in 1:4) {   # t0 set to zero
  qq_postpred(ypp = lvb_jags_out2_multi_t0_0_lnorm$sims.list$ypp[,,j], y=lvbdata$Length, 
              main=paste(modelnames[j], "- t0 set to zero"))
}

par(mfrow=c(2,2))   
# t0 free
caterpillar(lvb_jags_out2_multi_t0_free_lnorm, p="sig", xax=modelnames, xlab="t0 free")
caterpillar(lvb_jags_out2_multi_t0_free_lnorm, p="k", xax=modelnames, xlab="t0 free")
caterpillar(lvb_jags_out2_multi_t0_free_lnorm, p="t0", xax=modelnames, xlab="t0 free")
caterpillar(lvb_jags_out2_multi_t0_free_lnorm, p="L_inf", xax=modelnames, xlab="t0 free")

# t0 set to zero
caterpillar(lvb_jags_out2_multi_t0_0_lnorm, p="sig", xax=modelnames, xlab="t0 set to zero")
caterpillar(lvb_jags_out2_multi_t0_0_lnorm, p="k", xax=modelnames, xlab="t0 set to zero")
caterpillar(lvb_jags_out2_multi_t0_0_lnorm, p="t0", xax=modelnames, xlab="t0 set to zero")
caterpillar(lvb_jags_out2_multi_t0_0_lnorm, p="L_inf", xax=modelnames, xlab="t0 set to zero")

comparecat(list(lvb_jags_out2_multi_t0_free_lnorm, lvb_jags_out2_multi_t0_0_lnorm), p="sig")
comparecat(list(lvb_jags_out2_multi_t0_free_lnorm, lvb_jags_out2_multi_t0_0_lnorm), p="k")
comparecat(list(lvb_jags_out2_multi_t0_free_lnorm, lvb_jags_out2_multi_t0_0_lnorm), p="t0")
comparecat(list(lvb_jags_out2_multi_t0_free_lnorm, lvb_jags_out2_multi_t0_0_lnorm), p="L_inf")


for(j in 1:4) {
  par(mfcol = c(3,3))  
  plot_postpred(ypp = lvb_jags_out2_multi_t0_free_lnorm$sims.list$ypp[,,j], 
                y=lvbdata$Length, x=lvbdata$Age, col = j+1)
}
for(j in 1:4) {
  par(mfcol = c(3,3))  
  plot_postpred(ypp = lvb_jags_out2_multi_t0_0_lnorm$sims.list$ypp[,,j], 
                y=lvbdata$Length, x=lvbdata$Age, col = j+1)
}



doDIC <- TRUE

if(doDIC) {

##### comparing DICs across all candidate models
## trying multi-model
cat('model {
  for(i in 1:n) {
    L[i] ~ dlnorm(logmu[i], tau)
    ypp[i] ~ dlnorm(logmu[i], tau)
    logmu[i] <- log(whichmodel[1]*L_inf*(1-exp(-k*(a[i]-t0))) +      # VBM - VB
               whichmodel[2]*L_inf*(1-exp(-k*(a[i]-t0)))^P +    # GGM - generalized VB
               whichmodel[3]*L_inf/(1+exp(-k*(a[i]-t0))) +      # LGM - logistic growth model
               whichmodel[4]*L_inf*exp((-1/k)*exp(-k*(a[i]-t0)))) # PGM - Gompertz
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)#00)
  L_inf ~ dnorm(400, .0001)
  k ~ dnorm(0,.1)T(0,)
  t0 ~ dnorm(0,.01)T(,6)#<-0#
  P ~ dlnorm(0,1)
  
}', file="lvb_jags2_multi_t0_free_4DIC")

cat('model {
  for(i in 1:n) {
    L[i] ~ dlnorm(logmu[i], tau)
    ypp[i] ~ dlnorm(logmu[i], tau)
    logmu[i] <- log(whichmodel[1]*L_inf*(1-exp(-k*(a[i]-t0))) +      # VBM - VB
               whichmodel[2]*L_inf*(1-exp(-k*(a[i]-t0)))^P +    # GGM - generalized VB
               whichmodel[3]*L_inf/(1+exp(-k*(a[i]-t0))) +      # LGM - logistic growth model
               whichmodel[4]*L_inf*exp((-1/k)*exp(-k*(a[i]-t0)))) # PGM - Gompertz
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)#00)
  L_inf ~ dnorm(400, .0001)
  k ~ dnorm(0,.1)T(0,)
  t0 <-0 #~ dnorm(0,.01)T(,6)#
  P ~ dlnorm(0,1)
  
}', file="lvb_jags2_multi_t0_0_4DIC")

lvb_data2_multi_4DIC <- lvb_data2_multi
lvb_data2_multi_4DIC$L <- lvb_data2_multi_4DIC$L[,1]
modellist <- list()
for(j in 1:4) {
  niter <- 5000*1000             ##############    about 2 hours total at 1000k
  ncores <- min(parallel::detectCores()-1, 10)
  
  lvb_data2_multi_4DIC$whichmodel <- rep(0,4)
  lvb_data2_multi_4DIC$whichmodel[j] <- 1
  {
    tstart <- Sys.time()
    print(tstart)
    modellist[[j]] <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free_4DIC", data=lvb_data2_multi_4DIC,
                                     parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P"),
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
    
    tstart <- Sys.time()
    print(tstart)
    modellist[[j+4]] <- jagsUI::jags(model.file="lvb_jags2_multi_t0_0_4DIC", data=lvb_data2_multi_4DIC,
                                     parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P"),
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
}

par(mfrow=c(3,3))
for(j in 1:8) plotRhats(modellist[[j]])

DICs <- sapply(modellist, \(x) x$DIC)
plot(DICs-min(DICs))
# 3 wins 1:4 (weird)
# 6 wins 5:8 (no surprise)

DICs-min(DICs)
# # at 1000k, normal error
# [1]  16.49792  12.40554   0.00000  49.52038
# [5] 103.98221  23.69256  91.20001 127.61615

# # at 1000k, lognormal error
# [1]  16.72430  10.79417   0.00000  45.04177  
# [5]  98.03871  20.49830  85.91393 121.29690
}




do_kfold <- TRUE

if(do_kfold) {


## honestly the best way to do this will be loocv or kfcv (285 data points so probably kfcv)
# fold <- sample(1:5, length=nrow(lvbdata))
# for model in 1:8 {
#  for j in 1:5 {
#   make new data object where L[fold==j] <- NA
#   preds[fold==j] <- modelname$q50$L[fold==j]    # would need to add L to params (could even make that the only thing)
#  }
# rmse[model] <- rmse(preds, lvbdata$L)
# }

k <- 5
niter <- 1000*1000   # 42 minutes at 100k, 1.7 hours at 250k
fold <- sample(1:k, nrow(lvbdata), replace=TRUE)
rmses <- rep(NA, 8)
rmse <- function(x1, x2) sqrt(mean((x1-x2)^2, na.rm=TRUE))
tstart <- Sys.time()
print(tstart)
for(i_model in 1:4) {
  preds <- rep(NA, nrow(lvbdata))
  for(i_fold in 1:k) {
    lvb_data2_multi_4kfold <- lvb_data2_multi_4DIC
    lvb_data2_multi_4kfold$L[fold==i_fold] <- NA
    lvb_data2_multi_4kfold$whichmodel <- rep(0,4)
    lvb_data2_multi_4kfold$whichmodel[i_model] <- 1
    
    themodel <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free_4DIC", data=lvb_data2_multi_4kfold,
                             parameters.to.save="L", #c("mu","sig","L_inf","k","t0","ypp","P"),
                             n.chains=ncores, parallel=T, n.iter=niter,
                             n.burnin=niter/2, n.thin=niter/2000)
    preds[fold==i_fold] <- themodel$q50$L[fold==i_fold]
  }
  rmses[i_model] <- rmse(preds, lvbdata$Length)
  print(i_model)
}
for(i_model in 1:4) {
  preds <- rep(NA, nrow(lvbdata))
  for(i_fold in 1:k) {
    lvb_data2_multi_4kfold <- lvb_data2_multi_4DIC
    lvb_data2_multi_4kfold$L[fold==i_fold] <- NA
    lvb_data2_multi_4kfold$whichmodel <- rep(0,4)
    lvb_data2_multi_4kfold$whichmodel[i_model] <- 1
    
    themodel <- jagsUI::jags(model.file="lvb_jags2_multi_t0_0_4DIC", data=lvb_data2_multi_4kfold,
                             parameters.to.save="L", #c("mu","sig","L_inf","k","t0","ypp","P"),
                             n.chains=ncores, parallel=T, n.iter=niter,
                             n.burnin=niter/2, n.thin=niter/2000)
    preds[fold==i_fold] <- themodel$q50$L[fold==i_fold]
  }
  rmses[i_model+4] <- rmse(preds, lvbdata$Length)
  print(i_model+4)
}
print(Sys.time() - tstart)

plot(rmses)
# again, 3 wins 1:4 (weird!!)
# again, 6 wins 5:8 (no surprise)

rmses
# # lognormal error at 250k
# [1] 17.82475 17.66177 17.34923 18.51919 
# [5] 20.61240 17.93371 20.19443 21.43202

rmses - min(rmses)
# # normal error at 100k
# [1] 0.4101423 0.3192708 0.0000000 1.2566708
# [5] 3.5312422 0.6970438 3.0755605 4.4452475

# # lognormal error at 250k
# [1] 0.4755244 0.3125416 0.0000000 1.1699629 
# [5] 3.2631753 0.5844788 2.8452035 4.0827890

rmses/min(rmses)
# # lognormal error at 250k
# [1] 1.027409 1.018015 1.000000 1.067436
# [5] 1.188088 1.033689 1.163996 1.235330

}


####### looking at the length & age distrib of the spawning sample by date
par(mfrow=c(2,3))
with(spawn_sample,
     boxplot(Length ~ Date, main="All"))
with(subset(spawn_sample, Sex=="F"),
     boxplot(Length ~ Date, main="Females"))
with(subset(spawn_sample, Sex=="M"),
     boxplot(Length ~ Date, main="Males"))
with(spawn_sample,
     boxplot(Age ~ Date, main="All"))
with(subset(spawn_sample, Sex=="F"),
     boxplot(Age ~ Date, main="Females"))
with(subset(spawn_sample, Sex=="M"),
     boxplot(Age ~ Date, main="Males"))

par(mfrow=c(2,3))
with(spawn_sample, {
  idate <- 1
  for(datei in sort(unique(Date))) {
    plot(ecdf(Length[Date==datei]),
         add = (idate != 1),
         col = idate+1,
         xlab = "Length (mm FL)",
         main = "All")
    idate <- idate+1
  }
})
with(subset(spawn_sample, Sex=="F"), {
  idate <- 1
  for(datei in sort(unique(Date))) {
    plot(ecdf(Length[Date==datei]),
         add = (idate != 1),
         col = idate+1,
         xlab = "Length (mm FL)",
         main = "Females")
    idate <- idate+1
  }
})
with(subset(spawn_sample, Sex=="M"), {
  idate <- 1
  for(datei in sort(unique(Date))) {
    plot(ecdf(Length[Date==datei]),
         add = (idate != 1),
         col = idate+1,
         xlab = "Length (mm FL)",
         main = "Males")
    idate <- idate+1
  }
})

with(spawn_sample, {
  idate <- 1
  for(datei in sort(unique(Date))) {
    plot(ecdf(Age[Date==datei]),
         add = (idate != 1),
         col = idate+1,
         xlab = "Age",
         main = "All")
    idate <- idate+1
  }
})
with(subset(spawn_sample, Sex=="F"), {
  idate <- 1
  for(datei in sort(unique(Date))) {
    plot(ecdf(Age[Date==datei]),
         add = (idate != 1),
         col = idate+1,
         xlab = "Age",
         main = "Females")
    idate <- idate+1
  }
})
with(subset(spawn_sample, Sex=="M"), {
  idate <- 1
  for(datei in sort(unique(Date))) {
    plot(ecdf(Age[Date==datei]),
         add = (idate != 1),
         col = idate+1,
         xlab = "Age",
         main = "Males")
    idate <- idate+1
  }
})


length_kspmat <- age_kspmat <- matrix(nrow=3, ncol=3)

length_kspmat[1,1] <- with(spawn_sample, 
                           ks.test(Length[Date=="2023-12-19"], Length[Date=="2023-12-28"])$p.value)
length_kspmat[1,2] <- with(spawn_sample, 
                           ks.test(Length[Date=="2023-12-28"], Length[Date=="2024-01-03"])$p.value)
length_kspmat[1,3] <- with(spawn_sample, 
                           ks.test(Length[Date=="2023-12-19"], Length[Date=="2024-01-03"])$p.value)
length_kspmat[2,1] <- with(subset(spawn_sample, Sex=="M"), 
                           ks.test(Length[Date=="2023-12-19"], Length[Date=="2023-12-28"])$p.value)
length_kspmat[2,2] <- with(subset(spawn_sample, Sex=="M"), 
                           ks.test(Length[Date=="2023-12-28"], Length[Date=="2024-01-03"])$p.value)
length_kspmat[2,3] <- with(subset(spawn_sample, Sex=="M"), 
                           ks.test(Length[Date=="2023-12-19"], Length[Date=="2024-01-03"])$p.value)
length_kspmat[3,1] <- with(subset(spawn_sample, Sex=="F"), 
                           ks.test(Length[Date=="2023-12-19"], Length[Date=="2023-12-28"])$p.value)
length_kspmat[3,2] <- with(subset(spawn_sample, Sex=="F"), 
                           ks.test(Length[Date=="2023-12-28"], Length[Date=="2024-01-03"])$p.value)
length_kspmat[3,3] <- with(subset(spawn_sample, Sex=="F"), 
                           ks.test(Length[Date=="2023-12-19"], Length[Date=="2024-01-03"])$p.value)

age_kspmat[1,1] <- with(spawn_sample, 
                           ks.test(Age[Date=="2023-12-19"], Age[Date=="2023-12-28"])$p.value)
age_kspmat[1,2] <- with(spawn_sample, 
                           ks.test(Age[Date=="2023-12-28"], Age[Date=="2024-01-03"])$p.value)
age_kspmat[1,3] <- with(spawn_sample, 
                           ks.test(Age[Date=="2023-12-19"], Age[Date=="2024-01-03"])$p.value)
age_kspmat[2,1] <- with(subset(spawn_sample, Sex=="M"), 
                           ks.test(Age[Date=="2023-12-19"], Age[Date=="2023-12-28"])$p.value)
age_kspmat[2,2] <- with(subset(spawn_sample, Sex=="M"), 
                           ks.test(Age[Date=="2023-12-28"], Age[Date=="2024-01-03"])$p.value)
age_kspmat[2,3] <- with(subset(spawn_sample, Sex=="M"), 
                           ks.test(Age[Date=="2023-12-19"], Age[Date=="2024-01-03"])$p.value)
age_kspmat[3,1] <- with(subset(spawn_sample, Sex=="F"), 
                           ks.test(Age[Date=="2023-12-19"], Age[Date=="2023-12-28"])$p.value)
age_kspmat[3,2] <- with(subset(spawn_sample, Sex=="F"), 
                           ks.test(Age[Date=="2023-12-28"], Age[Date=="2024-01-03"])$p.value)
age_kspmat[3,3] <- with(subset(spawn_sample, Sex=="F"), 
                           ks.test(Age[Date=="2023-12-19"], Age[Date=="2024-01-03"])$p.value)

# still need to figure out how to apply sidak correction
# sidak_alpha <- 1-((1-.05)^(1/3))
sidak_alpha <- 1-((1-.1)^(1/3))
length_kspmat < sidak_alpha
age_kspmat < sidak_alpha



## to do still:
# VB growth curve

# VB growth params

# table: mean length for age bins

# KS or AD test?  - lengths and ages
# Variable 1	              Variable 2	              Variable 3	          Uncorrected KS P-value
# All fish from 12/19/23	  All fish from 12/28/23	  All fish from 1/3/24	0.41
# Males from 12/19/23	      Males from 12/28/23	      Males from 1/3/24	    0.20
# Females from 12/19/23	    Females from 12/28/23	    Females from 1/3/24	  0.72
# All males from all dates	All females from all dates		                  X.XX

# cohort thing??