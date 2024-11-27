library(tidyverse)
library(dsftools)
library(jagsUI)
library(jagshelper)


write_output <- FALSE
# write_output <- TRUE


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
write_output
if(write_output) write.csv(mnlength_tab, file="FDS_2024/R_output/Tab10a.csv")

# matching the format of the last report
# Age	n	Mean	SD	SE
mnlength_tab_2 <- data.frame(Age = paste0("'",agebreaks[-length(agebreaks)],# 
                                       "-", agebreaks[-1]-1),
                             n = mnlength_tab$n,
                             Mean = mnlength_tab$mn_length,
                             SD = unname(tapply(spawn_sample$Length, spawn_sample$agecut, sd, na.rm=TRUE)),
                             SE = mnlength_tab$se_length)
mnlength_tab_2
write_output
if(write_output) write.csv(mnlength_tab_2, file="FDS_2024/R_output/Tab10b.csv", quote=TRUE)

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
    # L_inf[j] ~ dnorm(400, .0001)
    L_inf[j] ~ dnorm(381.3, 1/(20.5^2))
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
    # L_inf[j] ~ dnorm(400, .0001)
    L_inf[j] ~ dnorm(381.3, 1/(20.5^2))
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
# run_anyway <- TRUE  # set this to TRUE to run the models anyway


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
  niter <- 1000*1000 # 18 min at 500k
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



##### Writing the output & saving to external files

write_output
# pnames <- c("L_inf","k","P","t0","sig")
if(write_output) {
  for(j in 1:4) {   # t0 free
    # make length age plots for all 4 models
    png(filename = paste0("FDS_2024/R_output/Fig3_", modelnames[j], ".png"), 
        res=300, width=8, height=8, units="in")
    
    par(family="serif")
    plot(lvbdata$Age, lvbdata$Length, xlab = "Age", ylab="Length (mm FL)")
    envelope(lvb_jags_out2_multi_t0_free_lnorm$sims.list$mufit[,,j], x=lvb_data2_multi$agefit, add=TRUE)
    # abline(h=lvb_jags_out2_multi_t0_free_lnorm$q50$L_inf[j], lty=2)
    # abline(h=lvb_jags_out2_multi_t0_free_lnorm$q2.5$L_inf[j], lty=3)
    # abline(h=lvb_jags_out2_multi_t0_free_lnorm$q97.5$L_inf[j], lty=3)
    
    dev.off()
    
  
    # make parameter tables for all 4 models
    outdf <- data.frame(Parameter = c("L_inf","k","P","t0","sig"),
               Estimate = c(lvb_jags_out2_multi_t0_free_lnorm$q50$L_inf[j],
                            lvb_jags_out2_multi_t0_free_lnorm$q50$k[j],
                            lvb_jags_out2_multi_t0_free_lnorm$q50$P[j],
                            lvb_jags_out2_multi_t0_free_lnorm$q50$t0[j],
                            lvb_jags_out2_multi_t0_free_lnorm$q50$sig[j]),
               SE = c(lvb_jags_out2_multi_t0_free_lnorm$sd$L_inf[j],
                            lvb_jags_out2_multi_t0_free_lnorm$sd$k[j],
                            lvb_jags_out2_multi_t0_free_lnorm$sd$P[j],
                            lvb_jags_out2_multi_t0_free_lnorm$sd$t0[j],
                            lvb_jags_out2_multi_t0_free_lnorm$sd$sig[j]),
               CI = c(
                 paste0("(", round(lvb_jags_out2_multi_t0_free_lnorm$q2.5$L_inf[j], 2), ", ",
                        round(lvb_jags_out2_multi_t0_free_lnorm$q97.5$L_inf[j], 2), ")"),
                 paste0("(", round(lvb_jags_out2_multi_t0_free_lnorm$q2.5$k[j], 2), ", ",
                        round(lvb_jags_out2_multi_t0_free_lnorm$q97.5$k[j], 2), ")"),
                 paste0("(", round(lvb_jags_out2_multi_t0_free_lnorm$q2.5$P[j], 2), ", ",
                        round(lvb_jags_out2_multi_t0_free_lnorm$q97.5$P[j], 2), ")"),
                 paste0("(", round(lvb_jags_out2_multi_t0_free_lnorm$q2.5$t0[j], 2), ", ",
                        round(lvb_jags_out2_multi_t0_free_lnorm$q97.5$t0[j], 2), ")"),
                 paste0("(", round(lvb_jags_out2_multi_t0_free_lnorm$q2.5$sig[j], 2), ", ",
                        round(lvb_jags_out2_multi_t0_free_lnorm$q97.5$sig[j], 2), ")")
               ))
    names(outdf)[4] <- "95% Credible Interval"
    if(j != 2) {
      outdf <- outdf[-3,]
    }
    write.csv(outdf, file = paste0("FDS_2024/R_output/Tab11_", modelnames[j], ".csv"))
  }
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



doDIC <- FALSE


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
  # L_inf ~ dnorm(400, .0001)
    L_inf ~ dnorm(381.3, 1/(20.5^2))
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
 # L_inf ~ dnorm(400, .0001)
    L_inf ~ dnorm(381.3, 1/(20.5^2))
  k ~ dnorm(0,.1)T(0,)
  t0 <-0 #~ dnorm(0,.01)T(,6)#
  P ~ dlnorm(0,1)
  
}', file="lvb_jags2_multi_t0_0_4DIC")

lvb_data2_multi_4DIC <- lvb_data2_multi
lvb_data2_multi_4DIC$L <- lvb_data2_multi_4DIC$L[,1]
ncores <- min(parallel::detectCores()-1, 10)



if(doDIC) {
  
modellist <- list()
for(j in 1:4) {
  niter <- 5000*1000             ##############    about 2 hours total at 1000k
  
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

# # at 5000k, lognormal error
# [1]  16.11816  10.85820   0.00000  45.45500  
# [5]  97.88702  20.36136  85.83256 121.43978
}




do_kfold <- FALSE

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
rmses <- maes <- rep(NA, 8)
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
  maes[i_model] <- jagshelper:::mae(preds, lvbdata$Length)
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
  maes[i_model+4] <- jagshelper:::mae(preds, lvbdata$Length)
  print(i_model+4)
}
print(Sys.time() - tstart)

plot(rmses)
# again, 3 wins 1:4 (weird!!)
# again, 6 wins 5:8 (no surprise)

rmses
# # lognormal error at 1000k
# [1] 17.66761 17.57152 17.25255 18.49503 
# [5] 20.68268 17.96346 20.23590 21.58844

rmses - min(rmses)
# # normal error at 100k
# [1] 0.4101423 0.3192708 0.0000000 1.2566708
# [5] 3.5312422 0.6970438 3.0755605 4.4452475

# # lognormal error at 1000k
# [1] 0.4150587 0.3189695 0.0000000 1.2424809 
# [5] 3.4301270 0.7109062 2.9833467 4.3358947

rmses/min(rmses)
# # lognormal error at 1000k
# [1] 1.024058 1.018488 1.000000 1.072017
# [5] 1.198819 1.041206 1.172922 1.251319

maes
# # lognormal error at 1000k
# [1] 12.76432 12.60870 12.30474 13.36021
# [5] 14.63153 12.95392 14.56734 14.83313
}


## making some visuals for model comparison
rmses <- c(17.66761, 17.57152, 17.25255, 18.49503, 20.68268, 17.96346, 20.23590, 21.58844)
maes <- c(12.76432, 12.60870, 12.30474, 13.36021, 14.63153, 12.95392, 14.56734, 14.83313)

par(mfrow=c(1,2))
par(mar=c(7,4,4,2))
barplot(rmses, 
        names.arg=rep(modelnames,2), las=2, 
        col=c(rep(6,4),rep(4,4)),
        ylab="Root Mean Square Prediction Error (mm)")
abline(h=min(rmses), lty=2)

# similar plot of mae
barplot(maes, 
        names.arg=rep(modelnames,2), las=2, 
        col=c(rep(6,4),rep(4,4)),
        ylab="Mean Absolute Prediction Error (mm)")
abline(h=min(maes), lty=2)


# caterpillar(cbind(lvb_jags_out2_multi_t0_free_lnorm$sims.list$sig,
#                   lvb_jags_out2_multi_t0_0_lnorm$sims.list$sig),
#             # ylim=c(0,0.065), 
#             col=c(rep(3,4),rep(4,4)),
#             xax=rep(modelnames,2), las=2, ylab="sigma")

caterpillar(cbind(lvb_jags_out2_multi_t0_free_lnorm$sims.list$t0,
                  lvb_jags_out2_multi_t0_0_lnorm$sims.list$t0),
            # ylim=c(0,0.065), 
            col=c(rep(6,4),rep(4,4)),
            xax=rep(modelnames,2), las=2, ylab="t0 (years)")

caterpillar(cbind(lvb_jags_out2_multi_t0_free_lnorm$sims.list$L_inf,
                  lvb_jags_out2_multi_t0_0_lnorm$sims.list$L_inf),
            # ylim=c(0,0.065), 
            col=c(rep(6,4),rep(4,4)),
            xax=rep(modelnames,2), las=2, ylab="L_inf (mm)")

# length - age with all 8 models, overlayed with Linf envelope
par(mfrow=c(2,2))
par(mar=c(5,4,4,2))
for(j in 1:4) {   # t0 free
  plot(lvbdata$Age, lvbdata$Length, main=paste(modelnames[j], "- t0 free"), 
       ylab="Fork Length (mm)", xlab="Age (years)",
       ylim=range(lvbdata$Length,625))
  envelope(lvb_jags_out2_multi_t0_free_lnorm$sims.list$mufit[,,j], x=lvb_data2_multi$agefit, add=T, col=6)
  envelope(replicate(2,lvb_jags_out2_multi_t0_free_lnorm$sims.list$L_inf[,j]), col=6, x=c(0,40), add=T)
}
for(j in 1:4) {   # t0 set to zero
  plot(lvbdata$Age, lvbdata$Length, main=paste(modelnames[j], "- t0 set to zero"), 
       ylab="Fork Length (mm)", xlab="Age (years)",
       ylim=range(lvbdata$Length,625))
  envelope(lvb_jags_out2_multi_t0_0_lnorm$sims.list$mufit[,,j], x=lvb_data2_multi$agefit, add=T, col=4)
  envelope(replicate(2,lvb_jags_out2_multi_t0_0_lnorm$sims.list$L_inf[,j]), col=4, x=c(0,40), add=T)
}



#### one more try with k-fold cross validation, using the new kfold function
niter <- 1000*100#0  # 20k in 15 min
{
tstart <- Sys.time()
ncores <- min(parallel::detectCores()-1, 10)
kfold_t0_free <- kfold(model.file="lvb_jags2_multi_t0_free_lnorm", data=lvb_data2_multi,
                       k=10, save_postpred = TRUE, fold_dims = 1, p="L",
                       n.chains=ncores, parallel=T, n.iter=niter,
                       n.burnin=niter/2, n.thin=niter/2000)
rmse_vec <- mae_vec <- NA
rmse_mat <- mae_mat <- matrix(NA, nrow=dim(kfold_t0_free$postpred_y)[1], ncol=8)
for(i in 1:4) {
  rmse_vec[i] <- jagshelper:::rmse(lvb_data2_multi$L[,i], kfold_t0_free$pred_y[,i])
  mae_vec[i] <- jagshelper:::mae(lvb_data2_multi$L[,i], kfold_t0_free$pred_y[,i])
  for(i_mcmc in 1:dim(kfold_t0_free$postpred_y)[1]) {
    rmse_mat[i_mcmc,i] <- jagshelper:::rmse(kfold_t0_free$postpred_y[i_mcmc,,i], 
                                            kfold_t0_free$data_y[,i])
    mae_mat[i_mcmc,i] <- jagshelper:::mae(kfold_t0_free$postpred_y[i_mcmc,,i], 
                                            kfold_t0_free$data_y[,i])
  }
}

kfold_t0_0 <- kfold(model.file="lvb_jags2_multi_t0_0_lnorm", data=lvb_data2_multi,
                    k=10, save_postpred = TRUE, fold_dims = 1, p="L",
                    n.chains=ncores, parallel=T, n.iter=niter,
                    n.burnin=niter/2, n.thin=niter/2000)
for(i in 1:4) {
  rmse_vec[i+4] <- jagshelper:::rmse(lvb_data2_multi$L[,i], kfold_t0_0$pred_y[,i])
  mae_vec[i+4] <- jagshelper:::mae(lvb_data2_multi$L[,i], kfold_t0_0$pred_y[,i])
  for(i_mcmc in 1:dim(kfold_t0_0$postpred_y)[1]) {
    rmse_mat[i_mcmc,i+4] <- jagshelper:::rmse(kfold_t0_0$postpred_y[i_mcmc,,i], 
                                              kfold_t0_0$data_y[,i])
    mae_mat[i_mcmc,i+4] <- jagshelper:::mae(kfold_t0_0$postpred_y[i_mcmc,,i], 
                                            kfold_t0_0$data_y[,i])
  }
}
print(Sys.time() - tstart)
}
caterpillar(rmse_mat)
caterpillar(mae_mat)


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

rownames(length_kspmat) <- rownames(age_kspmat) <- c("All fish", "Males", "Females")
colnames(length_kspmat) <- colnames(age_kspmat) <- c("12/19/23 vs. 12/28/23", 
                                                     "12/28/23 vs. 1/3/24", 
                                                     "12/19/23 vs. 1/3/24")

# still need to figure out how to apply sidak correction
alpha <- .1#0.05  # 0.05?  #0.1?
sidak_alpha <- 1-((1-alpha)^(1/3))
length_kspmat < sidak_alpha
age_kspmat < sidak_alpha

write_output
if(write_output) {
  write.csv(
rbind(c("Length", "", ""),
array(paste0(round(length_kspmat, 2),
             ifelse(length_kspmat < sidak_alpha, "*", "")), 
      dim=dim(length_kspmat),
      dimnames=list(rownames(length_kspmat), colnames(length_kspmat))),
c("", "", ""),c("Age", "", ""),

array(paste0(round(age_kspmat, 2),
             ifelse(age_kspmat < sidak_alpha, "*", "")), 
      dim=dim(age_kspmat),
      dimnames=list(rownames(age_kspmat), colnames(age_kspmat)))),
file="FDS_2024/R_output/Tab8.csv")
}




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



## cohort thing??
asl_2021 <- read_csv("FDS_2024/flat_data/ASL_2021.csv") %>% 
  as.data.frame %>%
  mutate(cohort = 2021 - `Age B`)

cohort21 <- asl_2021$cohort
cohort23 <- 2023 - spawn_sample$Age

thebreaks <- min(cohort21, cohort23, na.rm=TRUE):max(cohort21, cohort23, na.rm=TRUE)
hist(cohort21, breaks=thebreaks)
hist(cohort23, breaks=thebreaks)

count21 <- table(factor(cohort21, levels=thebreaks))
count23 <- table(factor(cohort23, levels=thebreaks))

p21 <- count21/sum(count21)
p23 <- count23/sum(count23)

cilo21 <- qbeta(p=0.025, shape1=count21+0.5, shape2=sum(count21)-count21+0.5)
cihi21 <- qbeta(p=0.975, shape1=count21+0.5, shape2=sum(count21)-count21+0.5)
cilo23 <- qbeta(p=0.025, shape1=count23+0.5, shape2=sum(count23)-count23+0.5)
cihi23 <- qbeta(p=0.975, shape1=count23+0.5, shape2=sum(count23)-count23+0.5)

cols <- c(2,4)
plot(x=thebreaks-.1, y=as.numeric(p21), ylim=c(0,max(cihi21, cihi23)), col=cols[1], pch=16)
segments(x0=thebreaks-.1, y0=cilo21, y1=cihi21, col=cols[1])
points(x=thebreaks+.1, y=as.numeric(p23), col=cols[2], pch=16)
segments(x0=thebreaks+.1, y0=cilo23, y1=cihi23, col=cols[2])

hist(cohort21, breaks=thebreaks, at=2*thebreaks)

aa <- 4
plot(NA, ylim=c(0, max(cihi21, cihi23)), xlim=c(1, aa*length(p21)), xaxt="n", 
     xlab="Cohort", ylab="Proportion")
axis(side=1, at=aa*seq_along(p21)-1, labels=thebreaks, las=2)
rect(xleft=aa*seq_along(p21)-2,
     xright=aa*seq_along(p21)-1,
     ybottom=0*p21, ytop=p21,
     col=adjustcolor(cols[1], alpha.f=.5),
     border=cols[1])
rect(xleft=aa*seq_along(p23)-1,
     xright=aa*seq_along(p23)-0,
     ybottom=0*p23, ytop=p23,
     col=adjustcolor(cols[2], alpha.f=.5),
     border=cols[2])
segments(x0=aa*seq_along(p21)-1.5, y0=cilo21, y1=cihi21, col=cols[1])
segments(x0=aa*seq_along(p21)-0.5, y0=cilo23, y1=cihi23, col=cols[2])
legend("topleft", fill=adjustcolor(cols, alpha.f=.5), border=cols, lty=1, col=cols,
       legend=c("2021 Sample (95% CI)", "2023 Sample (95% CI)"))
