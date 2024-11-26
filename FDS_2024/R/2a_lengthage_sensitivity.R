library(tidyverse)
library(dsftools)
library(jagsUI)
library(jagshelper)



# reading data
spawn_sample <- read_csv("FDS_2024/flat_data/spawn_sample.csv")%>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(Date = as.Date(Date, format="%m/%d/%Y"))

mcmc_overlap <- function(x1, x2) {
  # thebreaks <- seq(from=min(x1, x2, na.rm=TRUE), 
  #                  to=max(x1, x2, na.rm=TRUE), 
  #                  length.out=nbins+1)
  # t1 <- as.numeric(table(cut(x1, breaks=thebreaks)))
  # t2 <- as.numeric(table(cut(x2, breaks=thebreaks)))
  dd1 <- density(x1, from=min(x1, x2, na.rm=TRUE), 
                 to=max(x1, x2, na.rm=TRUE),
                 n=1024)
  dd2 <- density(x2, from=min(x1, x2, na.rm=TRUE), 
                 to=max(x1, x2, na.rm=TRUE),
                 n=1024)
  t1 <- dd1$y
  t2 <- dd2$y
  themin <- ifelse(t1<t2, t1, t2)
  themax <- ifelse(t1>t2, t1, t2)
  sum(themin)/sum(themax)
}




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
  # L_inf ~ dnorm(400, .0001)         # weak prior
  L_inf ~ dnorm(381.5, 1/(20.5^2))    # INFORMED PRIOR from prev study
  k ~ dnorm(0,.1)T(0,)
  t0 ~ dnorm(0,.01)T(,6)#<-0#
  P ~ dlnorm(0,1)
  
}', file="lvb_jags2_multi_t0_free_4DIC")


# bundle data to pass into JAGS
lvbdata <- subset(spawn_sample, !is.na(Age) & !is.na(Length))

lvb_data2_multi_4DIC <- list(L=lvbdata$Length,
                        a=lvbdata$`Age`,
                        n=nrow(lvbdata))#,
                        # agefit=seq(5,40,length.out=nfit),
                        # nfit=nfit,
                        # ppp=rep(1,4))

lvb_data2_multi_4DIC$whichmodel <- rep(0,4)

{
################################ starts here #################################

### switch for which model
lvb_data2_multi_4DIC$whichmodel <- rep(0,4)
lvb_data2_multi_4DIC$whichmodel[1] <- 1
# lvb_data2_multi_4DIC$whichmodel[2] <- 1
### switch for which model


## run the null model with all the data
niter <- 100*1000             ##############    2 minutes at 100k, 47 seconds at 40k
ncores <- min(parallel::detectCores()-1, 10)

{
tstart <- Sys.time()
print(tstart)
jags0 <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free_4DIC", data=lvb_data2_multi_4DIC,
                               parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P"),  # 
                               n.chains=ncores, parallel=T, n.iter=niter,
                               n.burnin=niter/2, n.thin=niter/2000)
print(Sys.time() - tstart)
}

plotRhats(jags0)
qq_postpred(jags0, p="ypp", y=lvb_data2_multi_4DIC$L)
par(mfcol=c(3,3))
plot_postpred(jags0, p="ypp", y=lvb_data2_multi_4DIC$L, x=lvb_data2_multi_4DIC$a)


# jags0$q50

# initialize places to put results
jags1_q50 <- jags1_sd <- jags1_overlap <- list()
for(ip in 1:length(jags0$q50)) {
  jags1_q50[[names(jags0$q50)[ip]]] <- rep(NA, nrow(lvbdata))
  jags1_sd[[names(jags0$q50)[ip]]] <- rep(NA, nrow(lvbdata))
  jags1_overlap[[names(jags0$q50)[ip]]] <- rep(NA, nrow(lvbdata))
}


## loop over dataset and NA-out each observation in turn 
# niter <- 20*1000             ##############    2.15 hours at 20k
{
  tstart <- Sys.time()
for(i in 1:nrow(lvbdata)) {
  lvb_data2_multi_4DIC1 <- lvb_data2_multi_4DIC
  lvb_data2_multi_4DIC1$L[i] <- NA
  
  jags1 <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free_4DIC", data=lvb_data2_multi_4DIC1,
                        parameters.to.save=c("mu","sig","L_inf","k","t0","P"),  # "ypp",
                        n.chains=ncores, parallel=T, n.iter=niter,
                        n.burnin=niter/2, n.thin=niter/2000,
                        verbose=FALSE)
  
  for(ip in 1:length(jags0$q50)) {
    jags1_q50[[names(jags0$q50)[ip]]][i] <- jags1$q50[[names(jags0$q50)[ip]]]
    jags1_sd[[names(jags0$q50)[ip]]][i] <- jags1$sd[[names(jags0$q50)[ip]]]
    jags1_overlap[[names(jags0$q50)[ip]]][i] <- mcmc_overlap(jags1$sims.list[[names(jags0$q50)[ip]]],
                                                             jags0$sims.list[[names(jags0$q50)[ip]]])
  }
  print(i)
}
  print(Sys.time() - tstart)
}
plotRhats(jags1)

par(mfrow=c(2,3))
for(ip in 1:length(jags1_q50)) {
  plot(jags1_q50[[ip]], main=names(jags1_q50)[ip], 
       col=adjustcolor("grey",alpha.f=.5),
       ylim=range(jags1_q50[[ip]], jags0$q50[[ip]], na.rm=TRUE), xlim=c(0,310))
  abline(h=jags0$q50[[ip]], lty=2)
  yy0 <- jags1_q50[[ip]]
  xx0 <- seq_along(yy0)
  yy1 <- yy0[rank(yy0)<=3]
  xx1 <- xx0[rank(yy0)<=3]
  yy2 <- yy0[rank(yy0)>sum(!is.na(yy0))-3]
  xx2 <- xx0[rank(yy0)>sum(!is.na(yy0))-3]
  text(x=xx1, y=yy1, pos=4, col=2,
       labels=xx1)
  text(x=xx2, y=yy2, pos=4, col=4,
       labels=xx2)
  points(xx1, yy1, col=2, pch=16)
  points(xx2, yy2, col=4, pch=16)
}
for(ip in 1:length(jags1_q50)) {
  yy0 <- jags1_q50[[ip]]
  these1 <- (rank(yy0)<=3) & !is.na(yy0)
  these2 <- (rank(yy0)>sum(!is.na(yy0))-3) & !is.na(yy0)
  plot(x=lvb_data2_multi_4DIC$a, y=lvb_data2_multi_4DIC$L, 
       col=adjustcolor("grey",alpha.f=.5),
       main=names(jags1_q50)[ip], xlim=c(5,40))
  points(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1],
         pch=16, col=2)
  points(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2],
         pch=16, col=4)
  text(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1], pos=4,
       labels=seq_along(yy0)[these1], col=2)
  text(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2], pos=4,
       labels=seq_along(yy0)[these2], col=4)
}
par(mfrow=c(2,3))
for(ip in 1:length(jags1_q50)) {
  plot(jags1_sd[[ip]], main=names(jags1_q50)[ip], 
       col=adjustcolor("grey",alpha.f=.5),
       ylim=range(jags1_sd[[ip]], jags0$sd[[ip]], na.rm=TRUE), xlim=c(0,310))
  abline(h=jags0$q50[[ip]], lty=2)
  yy0 <- jags1_sd[[ip]]
  xx0 <- seq_along(yy0)
  yy1 <- yy0[rank(yy0)<=3]
  xx1 <- xx0[rank(yy0)<=3]
  yy2 <- yy0[rank(yy0)>sum(!is.na(yy0))-3]
  xx2 <- xx0[rank(yy0)>sum(!is.na(yy0))-3]
  text(x=xx1, y=yy1, pos=4, col=2,
       labels=xx1)
  text(x=xx2, y=yy2, pos=4, col=4,
       labels=xx2)
  points(xx1, yy1, col=2, pch=16)
  points(xx2, yy2, col=4, pch=16)
}
for(ip in 1:length(jags1_q50)) {
  yy0 <- jags1_sd[[ip]]
  these1 <- (rank(yy0)<=3) & !is.na(yy0)
  these2 <- (rank(yy0)>sum(!is.na(yy0))-3) & !is.na(yy0)
  plot(x=lvb_data2_multi_4DIC$a, y=lvb_data2_multi_4DIC$L, 
       col=adjustcolor("grey",alpha.f=.5),
       main=names(jags1_q50)[ip], xlim=c(5,40))
  points(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1],
         pch=16, col=2)
  points(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2],
         pch=16, col=4)
  text(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1], pos=4,
       labels=seq_along(yy0)[these1], col=2)
  text(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2], pos=4,
       labels=seq_along(yy0)[these2], col=4)
}
par(mfrow=c(2,3))
for(ip in 1:length(jags1_q50)) {
  plot(jags1_overlap[[ip]], main=names(jags1_q50)[ip], 
       col=adjustcolor("grey",alpha.f=.5),
       ylim=range(jags1_overlap[[ip]], na.rm=TRUE), xlim=c(0,310), log="y")
  abline(h=jags0$q50[[ip]], lty=2)
  yy0 <- jags1_overlap[[ip]]
  xx0 <- seq_along(yy0)
  yy1 <- yy0[rank(yy0)<=3]
  xx1 <- xx0[rank(yy0)<=3]
  yy2 <- yy0[rank(yy0)>sum(!is.na(yy0))-3]
  xx2 <- xx0[rank(yy0)>sum(!is.na(yy0))-3]
  text(x=xx1, y=yy1, pos=4, col=2,
       labels=xx1)
  text(x=xx2, y=yy2, pos=4, col=4,
       labels=xx2)
  points(xx1, yy1, col=2, pch=16)
  points(xx2, yy2, col=4, pch=16)
}
for(ip in 1:length(jags1_q50)) {
  yy0 <- jags1_overlap[[ip]]
  these1 <- (rank(yy0)<=3) & !is.na(yy0)
  these2 <- (rank(yy0)>sum(!is.na(yy0))-3) & !is.na(yy0)
  plot(x=lvb_data2_multi_4DIC$a, y=lvb_data2_multi_4DIC$L, 
       col=adjustcolor("grey",alpha.f=.5),
       main=names(jags1_q50)[ip], xlim=c(5,40))
  points(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1],
         pch=16, col=2)
  points(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2],
         pch=16, col=4)
  text(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1], pos=4,
       labels=seq_along(yy0)[these1], col=2)
  text(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2], pos=4,
       labels=seq_along(yy0)[these2], col=4)
}




################################ starts here #################################

### switch for which model
lvb_data2_multi_4DIC$whichmodel <- rep(0,4)
# lvb_data2_multi_4DIC$whichmodel[1] <- 1
lvb_data2_multi_4DIC$whichmodel[2] <- 1
### switch for which model


## run the null model with all the data
# niter <- 100*1000             ##############    2 minutes at 100k, 47 seconds at 40k
# ncores <- min(parallel::detectCores()-1, 10)

{
  tstart <- Sys.time()
  print(tstart)
  mod2_jags0 <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free_4DIC", data=lvb_data2_multi_4DIC,
                        parameters.to.save=c("mu","sig","L_inf","k","t0","ypp","P"),  # 
                        n.chains=ncores, parallel=T, n.iter=niter,
                        n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

par(mfrow=c(1,1))
plotRhats(mod2_jags0)
qq_postpred(mod2_jags0, p="ypp", y=lvb_data2_multi_4DIC$L)
par(mfcol=c(3,3))
plot_postpred(mod2_jags0, p="ypp", y=lvb_data2_multi_4DIC$L, x=lvb_data2_multi_4DIC$a)

# mod2_jags0$q50

# initialize places to put results
mod2_jags1_q50 <- mod2_jags1_sd <- mod2_jags1_overlap <- list()
for(ip in 1:length(mod2_jags0$q50)) {
  mod2_jags1_q50[[names(mod2_jags0$q50)[ip]]] <- rep(NA, nrow(lvbdata))
  mod2_jags1_sd[[names(mod2_jags0$q50)[ip]]] <- rep(NA, nrow(lvbdata))
  mod2_jags1_overlap[[names(mod2_jags0$q50)[ip]]] <- rep(NA, nrow(lvbdata))
}


## loop over dataset and NA-out each observation in turn 
# niter <- 20*1000             ##############    2.15 hours at 20k
{
  tstart <- Sys.time()
  for(i in 1:nrow(lvbdata)) {
    lvb_data2_multi_4DIC1 <- lvb_data2_multi_4DIC
    lvb_data2_multi_4DIC1$L[i] <- NA
    
    mod2_jags1 <- jagsUI::jags(model.file="lvb_jags2_multi_t0_free_4DIC", data=lvb_data2_multi_4DIC1,
                          parameters.to.save=c("mu","sig","L_inf","k","t0","P"),  # "ypp",
                          n.chains=ncores, parallel=T, n.iter=niter,
                          n.burnin=niter/2, n.thin=niter/2000,
                          verbose=FALSE)
    
    for(ip in 1:length(mod2_jags0$q50)) {
      mod2_jags1_q50[[names(mod2_jags0$q50)[ip]]][i] <- mod2_jags1$q50[[names(mod2_jags0$q50)[ip]]]
      mod2_jags1_sd[[names(mod2_jags0$q50)[ip]]][i] <- mod2_jags1$sd[[names(mod2_jags0$q50)[ip]]]
      mod2_jags1_overlap[[names(mod2_jags0$q50)[ip]]][i] <- mcmc_overlap(mod2_jags1$sims.list[[names(mod2_jags0$q50)[ip]]],
                                                               mod2_jags0$sims.list[[names(mod2_jags0$q50)[ip]]])
    }
    print(i)
  }
  print(Sys.time() - tstart)
}
plotRhats(mod2_jags1)

par(mfrow=c(2,3))
for(ip in 1:length(mod2_jags1_q50)) {
  plot(mod2_jags1_q50[[ip]], main=names(mod2_jags1_q50)[ip], 
       col=adjustcolor("grey",alpha.f=.5),
       ylim=range(mod2_jags1_q50[[ip]], mod2_jags0$q50[[ip]], na.rm=TRUE), xlim=c(0,310))
  abline(h=mod2_jags0$q50[[ip]], lty=2)
  yy0 <- mod2_jags1_q50[[ip]]
  xx0 <- seq_along(yy0)
  yy1 <- yy0[rank(yy0)<=3]
  xx1 <- xx0[rank(yy0)<=3]
  yy2 <- yy0[rank(yy0)>sum(!is.na(yy0))-3]
  xx2 <- xx0[rank(yy0)>sum(!is.na(yy0))-3]
  text(x=xx1, y=yy1, pos=4, col=2,
       labels=xx1)
  text(x=xx2, y=yy2, pos=4, col=4,
       labels=xx2)
  points(xx1, yy1, col=2, pch=16)
  points(xx2, yy2, col=4, pch=16)
}
for(ip in 1:length(mod2_jags1_q50)) {
  yy0 <- mod2_jags1_q50[[ip]]
  these1 <- (rank(yy0)<=3) & !is.na(yy0)
  these2 <- (rank(yy0)>sum(!is.na(yy0))-3) & !is.na(yy0)
  plot(x=lvb_data2_multi_4DIC$a, y=lvb_data2_multi_4DIC$L, 
       col=adjustcolor("grey",alpha.f=.5),
       main=names(mod2_jags1_q50)[ip], xlim=c(5,40))
  points(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1],
         pch=16, col=2)
  points(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2],
         pch=16, col=4)
  text(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1], pos=4,
       labels=seq_along(yy0)[these1], col=2)
  text(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2], pos=4,
       labels=seq_along(yy0)[these2], col=4)
}
par(mfrow=c(2,3))
for(ip in 1:length(mod2_jags1_q50)) {
  plot(mod2_jags1_sd[[ip]], main=names(mod2_jags1_q50)[ip], 
       col=adjustcolor("grey",alpha.f=.5),
       ylim=range(mod2_jags1_sd[[ip]], mod2_jags0$sd[[ip]], na.rm=TRUE), xlim=c(0,310))
  abline(h=mod2_jags0$q50[[ip]], lty=2)
  yy0 <- mod2_jags1_sd[[ip]]
  xx0 <- seq_along(yy0)
  yy1 <- yy0[rank(yy0)<=3]
  xx1 <- xx0[rank(yy0)<=3]
  yy2 <- yy0[rank(yy0)>sum(!is.na(yy0))-3]
  xx2 <- xx0[rank(yy0)>sum(!is.na(yy0))-3]
  text(x=xx1, y=yy1, pos=4, col=2,
       labels=xx1)
  text(x=xx2, y=yy2, pos=4, col=4,
       labels=xx2)
  points(xx1, yy1, col=2, pch=16)
  points(xx2, yy2, col=4, pch=16)
}
for(ip in 1:length(mod2_jags1_q50)) {
  yy0 <- mod2_jags1_sd[[ip]]
  these1 <- (rank(yy0)<=3) & !is.na(yy0)
  these2 <- (rank(yy0)>sum(!is.na(yy0))-3) & !is.na(yy0)
  plot(x=lvb_data2_multi_4DIC$a, y=lvb_data2_multi_4DIC$L, 
       col=adjustcolor("grey",alpha.f=.5),
       main=names(mod2_jags1_q50)[ip], xlim=c(5,40))
  points(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1],
         pch=16, col=2)
  points(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2],
         pch=16, col=4)
  text(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1], pos=4,
       labels=seq_along(yy0)[these1], col=2)
  text(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2], pos=4,
       labels=seq_along(yy0)[these2], col=4)
}

par(mfrow=c(2,3))
for(ip in 1:length(mod2_jags1_q50)) {
  plot(mod2_jags1_overlap[[ip]], main=names(mod2_jags1_q50)[ip], 
       col=adjustcolor("grey",alpha.f=.5),
       ylim=range(mod2_jags1_overlap[[ip]], na.rm=TRUE), xlim=c(0,310), log="y")
  abline(h=mod2_jags0$q50[[ip]], lty=2)
  yy0 <- mod2_jags1_overlap[[ip]]
  xx0 <- seq_along(yy0)
  yy1 <- yy0[rank(yy0)<=3]
  xx1 <- xx0[rank(yy0)<=3]
  yy2 <- yy0[rank(yy0)>sum(!is.na(yy0))-3]
  xx2 <- xx0[rank(yy0)>sum(!is.na(yy0))-3]
  text(x=xx1, y=yy1, pos=4, col=2,
       labels=xx1)
  text(x=xx2, y=yy2, pos=4, col=4,
       labels=xx2)
  points(xx1, yy1, col=2, pch=16)
  points(xx2, yy2, col=4, pch=16)
}
for(ip in 1:length(mod2_jags1_q50)) {
  yy0 <- mod2_jags1_overlap[[ip]]
  these1 <- (rank(yy0)<=3) & !is.na(yy0)
  these2 <- (rank(yy0)>sum(!is.na(yy0))-3) & !is.na(yy0)
  plot(x=lvb_data2_multi_4DIC$a, y=lvb_data2_multi_4DIC$L, 
       col=adjustcolor("grey",alpha.f=.5),
       main=names(mod2_jags1_q50)[ip], xlim=c(5,40))
  points(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1],
         pch=16, col=2)
  points(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2],
         pch=16, col=4)
  text(x=lvb_data2_multi_4DIC$a[these1], y=lvb_data2_multi_4DIC$L[these1], pos=4,
       labels=seq_along(yy0)[these1], col=2)
  text(x=lvb_data2_multi_4DIC$a[these2], y=lvb_data2_multi_4DIC$L[these2], pos=4,
       labels=seq_along(yy0)[these2], col=4)
}
}







# could also do z test (ish) based on (mn0-mn1)/sqrt(sd0^2+sd1^2)
jags1_pval <- jags1_q50
for(ip in 1:length(jags1_q50)) {
  # for(j in 1:length(jags1_pval[[names(jags1_q50)[ip]]])) {
    jags1_pval[[names(jags1_q50)[ip]]] <- 2*(1-pnorm(abs(
      (jags0$q50[[names(jags1_q50)[ip]]] - jags1_q50[[names(jags1_q50)[ip]]]) / # [j]
      sqrt(jags0$sd[[names(jags1_q50)[ip]]]^2 + jags1_sd[[names(jags1_q50)[ip]]]^2))))
}
for(ip in 1:length(jags1_q50)) {
  plot(jags1_pval[[names(jags1_q50)[ip]]], main=names(jags1_q50)[ip], log="y")
}

# do mcmc_overlap thing for sure, either store all posts or make an object