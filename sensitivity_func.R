library(jagshelper)
library(jagsUI)

# specify model, which is written to a temporary file
NAME_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0 + b1*x[i] + a[grp[i]]
  }

  for(j in 1:ngrp) {
    a[j] ~ dnorm(0, tau_a)
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  tau_a <- pow(sig_a, -2)
  sig_a ~ dunif(0, 10)
}', file=NAME_jags)


# simulate data to go with the example model
n <- 60
x <- rnorm(n, sd=3)
grp <- sample(1:3, n, replace=T)
y <- rnorm(n, mean=grp-x)

# bundle data to pass into JAGS
NAME_data <- list(x=x,
                  y=y,
                  n=length(x),
                  grp=as.numeric(as.factor(grp)),
                  ngrp=length(unique(grp)))

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  NAME_jags_out <- jagsUI::jags(model.file=NAME_jags, data=NAME_data,
                                parameters.to.save=c("b0","b1","sig","a","sig_a"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(NAME_jags_out)
plotRhats(NAME_jags_out)
traceworstRhat(NAME_jags_out, parmfrow = c(3, 3))


str(NAME_data)
str(NAME_jags_out$q50)

sensitivity <- function(model.file, data, 
                        p_data, 
                        p_inference,
                        # save_postpred = FALSE, 
                        ...) {
  
  # fit it once for the thing to compare it to as well as structure of output
  jags0 <- jagsUI::jags(model.file = model.file, data = data, 
                        parameters.to.save = p_inference, 
                        verbose = FALSE, 
                        codaOnly = FALSE, bugs.format = FALSE, 
                        n.chains=n.chains, parallel=parallel, n.iter=n.iter,
                        n.burnin=n.burnin, n.thin=n.thin)
                        # ... = ...)
  
  # read the dimensions of data and parameters to save
  if(is.null(dim(data[[p_data]]))) {
    dimdata <- length(data[[p_data]])
  } else {
    dimdata <- dim(data[[p_data]])
  }
  
  # set up result objects for q50 and sd
  np <- nbyname(jags0)
  q50_list <- sd_list <- list()
  for(i_p in seq_along(p_inference)) {
    q50_list[[p_inference[i_p]]] <- 
      sd_list[[p_inference[i_p]]] <- 
      array(dim=c(dimdata, np[[p_inference[i_p]]]))
  }
  
  # loop over data (maybe different versions for different dim lengths??)
  if(length(dimdata) == 1) {
    data1 <- data
    for(i1 in 1:dimdata[1]) {
      data1[[p_data]][i1] <- NA
      jags1 <- jagsUI::jags(model.file = model.file, data = data1, 
                            parameters.to.save = p_inference, 
                            verbose = FALSE, 
                            codaOnly = FALSE, bugs.format = FALSE, 
                            n.chains=n.chains, parallel=parallel, n.iter=n.iter,
                            n.burnin=n.burnin, n.thin=n.thin)
                            # ... = ...)
      for(i_p in seq_along(p_inference)) {
        q50_list[[p_inference[i_p]]][i1]  ##### uggh make it different for each possible dims of q50_list[[thing]]
      }
    }
  }
  if(length(dimdata) == 2) {
    
  }
  if(length(dimdata) == 3) {
    
  }
  if(length(dimdata) == 4) {
    
  }
  #   make new data object with data NA'd out
  #   run jags
  #   loop over each data & put the stuff in the appropriate place
}

sensitivity(model.file=NAME_jags, 
            data=NAME_data, 
            p_data="y", 
            p_inference=c("b0","b1","sig","a","sig_a"),
            n.chains=ncores, parallel=T, n.iter=niter,
            n.burnin=niter/2, n.thin=niter/2000)


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
mcmc_overlap(x1 = rnorm(10000, mean=0, sd=1), x2 = rnorm(10000, mean=4, sd=1))
