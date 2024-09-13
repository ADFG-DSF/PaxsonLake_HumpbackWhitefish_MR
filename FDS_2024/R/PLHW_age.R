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




## trying multi-model
cat('model {
  for(i in 1:n) {
    for(j in 1:4) {
      L[i,j] ~ dnorm(mu[i,j], tau[j])
      ypp[i,j] ~ dnorm(mu[i,j], tau[j])
    }
    mu[i,1] <- L_inf[1]*(1-exp(-k[1]*(a[i]-t0[1])))      # VBM - VB
    mu[i,2] <- L_inf[2]*(1-exp(-k[2]*(a[i]-t0[2])))^P[2]    # GGM - generalized VB
    mu[i,3] <- L_inf[3]/(1+exp(-k[3]*(a[i]-t0[3])))      # LGM - logistic growth model
    mu[i,4] <- L_inf[4]*exp((-1/k[4])*exp(-k[4]*(a[i]-t0[4]))) # PGM - Gompertz
  }

  for(j in 1:4) {
    tau[j] <- pow(sig[j], -2)
    sig[j] ~ dunif(0, 1000)
    L_inf[j] ~ dnorm(400, .0001)
    k[j] ~ dnorm(0,.1)T(0,)
    t0[j] ~ dnorm(0,.01)T(,6)#<-0#
    P[j] ~ dlnorm(0,1)
  }
  
  for(ifit in 1:nfit) {
    mufit[ifit,1] <- L_inf[1]*(1-exp(-k[1]*(agefit[ifit]-t0[1])))      # VBM - VB
    mufit[ifit,2] <- L_inf[2]*(1-exp(-k[2]*(agefit[ifit]-t0[2])))^P[2]    # GGM - generalized VB
    mufit[ifit,3] <- L_inf[3]/(1+exp(-k[3]*(agefit[ifit]-t0[3])))      # LGM - logistic growth model
    mufit[ifit,4] <- L_inf[4]*exp((-1/k[4])*exp(-k[4]*(agefit[ifit]-t0[4]))) # PGM - Gompertz
    for(j in 1:4) {
      ppfit[ifit,j] ~ dnorm(mufit[ifit,j], tau[j])
    }
  }
  
}', file="lvb_jags2_multi_t0_free")

cat('model {
  for(i in 1:n) {
    for(j in 1:4) {
      L[i,j] ~ dnorm(mu[i,j], tau[j])
      ypp[i,j] ~ dnorm(mu[i,j], tau[j])
    }
    mu[i,1] <- L_inf[1]*(1-exp(-k[1]*(a[i]-t0[1])))      # VBM - VB
    mu[i,2] <- L_inf[2]*(1-exp(-k[2]*(a[i]-t0[2])))^P[2]    # GGM - generalized VB
    mu[i,3] <- L_inf[3]/(1+exp(-k[3]*(a[i]-t0[3])))      # LGM - logistic growth model
    mu[i,4] <- L_inf[4]*exp((-1/k[4])*exp(-k[4]*(a[i]-t0[4]))) # PGM - Gompertz
  }

  for(j in 1:4) {
    tau[j] <- pow(sig[j], -2)
    sig[j] ~ dunif(0, 1000)
    L_inf[j] ~ dnorm(400, .0001)
    k[j] ~ dnorm(0,.1)T(0,)
    t0[j] <- 0 #~ dnorm(0,.001)T(,6)#
    P[j] ~ dlnorm(0,1)
  }
  
  for(ifit in 1:nfit) {
    mufit[ifit,1] <- L_inf[1]*(1-exp(-k[1]*(agefit[ifit]-t0[1])))      # VBM - VB
    mufit[ifit,2] <- L_inf[2]*(1-exp(-k[2]*(agefit[ifit]-t0[2])))^P[2]    # GGM - generalized VB
    mufit[ifit,3] <- L_inf[3]/(1+exp(-k[3]*(agefit[ifit]-t0[3])))      # LGM - logistic growth model
    mufit[ifit,4] <- L_inf[4]*exp((-1/k[4])*exp(-k[4]*(agefit[ifit]-t0[4]))) # PGM - Gompertz
    for(j in 1:4) {
      ppfit[ifit,j] ~ dnorm(mufit[ifit,j], tau[j])
    }
  }
  
}', file="lvb_jags2_multi_t0_0")



# ## trying multi-model  --- LOGNORMAL ERROR
# cat('model {
#   for(i in 1:n) {
#     for(j in 1:4) {
#       L[i,j] ~ dlnorm(logmu[i,j], tau[j])
#       ypp[i,j] ~ dlnorm(logmu[i,j], tau[j])
#       mu[i,j] <- exp(logmu[i,j])
#     }
#     logmu[i,1] <- log(L_inf[1]*(1-exp(-k[1]*(a[i]-t0[1]))))      # VBM - VB
#     logmu[i,2] <- log(L_inf[2]*(1-exp(-k[2]*(a[i]-t0[2])))^P[2])    # GGM - generalized VB
#     logmu[i,3] <- log(L_inf[3]/(1+exp(-k[3]*(a[i]-t0[3]))))      # LGM - logistic growth model
#     logmu[i,4] <-log( L_inf[4]*exp((-1/k[4])*exp(-k[4]*(a[i]-t0[4])))) # PGM - Gompertz
#   }
# 
#   for(j in 1:4) {
#     tau[j] <- pow(sig[j], -2)
#     sig[j] ~ dunif(0, 10)#00)
#     L_inf[j] ~ dnorm(400, .0001)
#     k[j] ~ dnorm(0,.1)T(0,)
#     t0[j] ~ dnorm(0,.001)T(,6)#<-0#
#     P[j] ~ dlnorm(0,1)
#   }
#   
#   for(ifit in 1:nfit) {
#     mufit[ifit,1] <- L_inf[1]*(1-exp(-k[1]*(agefit[ifit]-t0[1])))      # VBM - VB
#     mufit[ifit,2] <- L_inf[2]*(1-exp(-k[2]*(agefit[ifit]-t0[2])))^P[2]    # GGM - generalized VB
#     mufit[ifit,3] <- L_inf[3]/(1+exp(-k[3]*(agefit[ifit]-t0[3])))      # LGM - logistic growth model
#     mufit[ifit,4] <- L_inf[4]*exp((-1/k[4])*exp(-k[4]*(agefit[ifit]-t0[4]))) # PGM - Gompertz
#     for(j in 1:4) {
#       ppfit[ifit,j] ~ dnorm(mufit[ifit,j], tau[j])
#     }
#   }
#   
# }', file="lvb_jags2_multi")


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

# JAGS controls - t0 free
niter <- 1000000 # 56 min at 1000k
ncores <- min(parallel::detectCores()-1, 10)

{
  tstart <- Sys.time()
  print(tstart)
  lvb_jags_out2_multi_t0_free <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free", data=lvb_data2_multi,
                                      parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P","mufit","ppfit","logmu","whichmodel","pp"),
                                      n.chains=ncores, parallel=T, n.iter=niter,
                                      n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
nbyname(lvb_jags_out2_multi_t0_free)
plotRhats(lvb_jags_out2_multi_t0_free)
par(mfrow=c(2,2))
traceworstRhat(lvb_jags_out2_multi_t0_free)

save(lvb_jags_out2_multi_t0_free, file="FDS_2024/posts/lvb_jags_out2_multi_t0_free.Rdata")



# JAGS controls - t0 set to zero
niter <- 500000 # 18 min at 500k
ncores <- min(parallel::detectCores()-1, 10)

{
  tstart <- Sys.time()
  print(tstart)
  lvb_jags_out2_multi_t0_0 <- jagsUI::jags(model.file="lvb_jags2_multi_t0_0", data=lvb_data2_multi,
                                           parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P","mufit","ppfit","logmu","whichmodel","pp"),
                                           n.chains=ncores, parallel=T, n.iter=niter,
                                           n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
nbyname(lvb_jags_out2_multi_t0_0)
plotRhats(lvb_jags_out2_multi_t0_0)
par(mfrow=c(2,2))
traceworstRhat(lvb_jags_out2_multi_t0_0)

save(lvb_jags_out2_multi_t0_0, file="FDS_2024/posts/lvb_jags_out2_multi_t0_0.Rdata")




par(mfrow=c(2,2))
for(j in 1:4) {   # t0 free
  plot(lvbdata$Age, lvbdata$Length)
  envelope(lvb_jags_out2_multi_t0_free$sims.list$mufit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}
for(j in 1:4) {   # t0 set to zero
  plot(lvbdata$Age, lvbdata$Length)
  envelope(lvb_jags_out2_multi_t0_0$sims.list$mufit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}

for(j in 1:4) {   # t0 free
  plot(lvbdata$Age, lvbdata$Length)
  envelope(lvb_jags_out2_multi_t0_free$sims.list$ppfit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}
for(j in 1:4) {   # t0 set to zero
  plot(lvbdata$Age, lvbdata$Length)
  envelope(lvb_jags_out2_multi_t0_0$sims.list$ppfit[,,j], x=lvb_data2_multi$agefit, add=T, col=j+1)
}

for(j in 1:4) {   # t0 free
  qq_postpred(ypp = lvb_jags_out2_multi_t0_free$sims.list$ypp[,,j], y=lvbdata$Length)
}
for(j in 1:4) {   # t0 set to zero
  qq_postpred(ypp = lvb_jags_out2_multi_t0_0$sims.list$ypp[,,j], y=lvbdata$Length)
}

par(mfrow=c(2,2))   
# t0 free
caterpillar(lvb_jags_out2_multi_t0_free, p="sig")
caterpillar(lvb_jags_out2_multi_t0_free, p="k")
caterpillar(lvb_jags_out2_multi_t0_free, p="t0")
caterpillar(lvb_jags_out2_multi_t0_free, p="L_inf")

# t0 set to zero
caterpillar(lvb_jags_out2_multi_t0_0, p="sig")
caterpillar(lvb_jags_out2_multi_t0_0, p="k")
caterpillar(lvb_jags_out2_multi_t0_0, p="t0")
caterpillar(lvb_jags_out2_multi_t0_0, p="L_inf")


for(j in 1:4) {
  par(mfcol = c(3,3))  
  plot_postpred(ypp = lvb_jags_out2_multi_t0_free$sims.list$ypp[,,j], 
                y=lvbdata$Length, x=lvbdata$Age, col = j+1)
  par(mfcol = c(3,3))  
  plot_postpred(ypp = lvb_jags_out2_multi_t0_0$sims.list$ypp[,,j], 
                y=lvbdata$Length, x=lvbdata$Age, col = j+1)
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