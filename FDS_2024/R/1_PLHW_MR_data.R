library(tidyverse)
library(recapr)
library(dsftools)

### read data
Event1 <- read_csv("FDS_2024/flat_data/Event1.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`))

Event2 <- read_csv("FDS_2024/flat_data/Event2.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`))

Event2_recaps <- subset(Event2, !is.na(`Tag Number`)) %>%
  arrange(`Tag Number`)
Event1_recaps <- subset(Event1, `Tag Number` %in% Event2_recaps$`Tag Number`) %>%
  arrange(`Tag Number`)



### reproduce KS tests
# reproduce growth correction
xreg <- Event2_recaps$`Fork Length (mm)`[Event2_recaps$`Tag Number` != "TL"]
yreg <- Event1_recaps$`Fork Length (mm)`
lm_growth <- lm(yreg ~ xreg)
lm_predict <- predict(lm_growth, newdata = data.frame(xreg=Event2$`Fork Length (mm)`))
plot(Event2$`Fork Length Corrected (mm)`, lm_predict)
abline(0,1)
# I think it makes more sense to predict event 1 size from event 2, but this is 
# essentially fine

# censor < 345mm
table(is.na(Event1$`Fork Length (mm)`))
table(is.na(Event2$`Fork Length Corrected (mm)`))
table(is.na(Event1_recaps$`Fork Length (mm)`))
table(is.na(Event2_recaps$`Fork Length Corrected (mm)`))

Event1 <- filter(Event1, `Fork Length (mm)` >= 345)
Event1_recaps <- filter(Event1_recaps, `Fork Length (mm)` >= 345)
Event2 <- filter(Event2, `Fork Length Corrected (mm)` >= 345)
Event2_recaps <- filter(Event2_recaps, `Fork Length Corrected (mm)` >= 345)

ks.test(Event1$`Fork Length (mm)`, 
        Event1_recaps$`Fork Length (mm)`)
ks.test(Event2$`Fork Length Corrected (mm)`, 
        Event2_recaps$`Fork Length Corrected (mm)`)

ks.test(Event1$`Fork Length (mm)` %s_l% 400, 
        Event1_recaps$`Fork Length (mm)` %s_l% 400)
ks.test(Event2$`Fork Length Corrected (mm)` %s_l% 400, 
        Event2_recaps$`Fork Length Corrected (mm)` %s_l% 400)

ks.test(Event1$`Fork Length (mm)` %s_geq% 400, 
        Event1_recaps$`Fork Length (mm)` %s_geq% 400)
ks.test(Event2$`Fork Length Corrected (mm)` %s_geq% 400, 
        Event2_recaps$`Fork Length Corrected (mm)` %s_geq% 400)

ksplot <- function(x1, x2, legend=c("x1","x2"), main="", col=c(1,1), lty=c(1,1), xlab="") {
  d1 <- density(x1, na.rm=TRUE)
  d2 <- density(x2, na.rm=TRUE)
  plot(d1, main=main, col=col[1], lty=lty[1],
       xlim=range(d1$x, d2$x), ylim=range(d1$y, d2$y), xlab=xlab)
  lines(density(x2, na.rm=TRUE), col=col[2], lty=lty[2])
  legend("topright", lty=lty, col=col,
         legend=paste0(legend, " (n=",c(sum(!is.na(x1)), sum(!is.na(x2))),")"))
  
  ksks <- suppressWarnings(ks.test(x1, x2))
  plot(ecdf(x1), main=c(main,
                        paste0("D=", signif(ksks$statistic, digits=3), ", ",
                               "pval=", signif(ksks$p.value, digits=3), collapse=NULL)),
       col=col[1], lty=lty[1], xlab=xlab)
  plot(ecdf(x2), col=col[2], lty=lty[2], add=TRUE)
  legend("bottomright", lty=lty, col=col,
         legend=paste0(legend, " (n=",c(sum(!is.na(x1)), sum(!is.na(x2))),")"))
}
par(mfrow=c(2,2))
ksplot(Event1$`Fork Length (mm)`, 
        Event1_recaps$`Fork Length (mm)`,
       xlab="mm FL", main="Event 1",
       col=c(1,2), legend=c("All","Recaps"))
ksplot(Event2$`Fork Length Corrected (mm)`, 
        Event2_recaps$`Fork Length Corrected (mm)`,
       xlab="mm FL", main="Event 2",
       col=c(1,4), legend=c("All","Recaps"))

ksplot(Event1$`Fork Length (mm)` %s_l% 400, 
        Event1_recaps$`Fork Length (mm)` %s_l% 400,
       xlab="mm FL", main="Event 1 - Small (<400mm)",
       col=c(1,2), legend=c("All","Recaps"))
ksplot(Event2$`Fork Length Corrected (mm)` %s_l% 400, 
        Event2_recaps$`Fork Length Corrected (mm)` %s_l% 400,
       xlab="mm FL", main="Event 2 - Small (<400mm)",
       col=c(1,4), legend=c("All","Recaps"))

ksplot(Event1$`Fork Length (mm)` %s_geq% 400, 
        Event1_recaps$`Fork Length (mm)` %s_geq% 400,
       xlab="mm FL", main="Event 1 - Large (>=400mm)",
       col=c(1,2), legend=c("All","Recaps"))
ksplot(Event2$`Fork Length Corrected (mm)` %s_geq% 400, 
        Event2_recaps$`Fork Length Corrected (mm)` %s_geq% 400,
       xlab="mm FL", main="Event 2 - Large (>=400mm)",
       col=c(1,4), legend=c("All","Recaps"))
## I get different numbers in my KS tests, but the inferences are the same


### reproduce X2 tests
consistencytest(n1 = table(Event1$Area),
                n2 = c(0, 0, 0, table(Event2$Area)),
                m2strata1 = as.numeric(as.factor(Event1_recaps$Area)),
                m2strata2 = rep(4, nrow(Event1_recaps)))

consistencytest(n1 = table(Event1$Area[Event1$`Fork Length (mm)` < 400]),
                n2 = c(0, 0, 0, table(Event2$Area[Event2$`Fork Length Corrected (mm)` < 400])),
                m2strata1 = as.numeric(as.factor(Event1_recaps$Area[Event1_recaps$`Fork Length (mm)` < 400])),
                m2strata2 = rep(4, sum(Event1_recaps$`Fork Length (mm)` < 400)))

consistencytest(n1 = table(Event1$Area[Event1$`Fork Length (mm)` >= 400]),
                n2 = c(0, 0, 0, table(Event2$Area[Event2$`Fork Length Corrected (mm)` >= 400])),
                m2strata1 = as.numeric(as.factor(Event1_recaps$Area[Event1_recaps$`Fork Length (mm)` >= 400]))+1,
                m2strata2 = rep(4, sum(Event1_recaps$`Fork Length (mm)` >= 400)))



### reproduce abundance estimation
NChapman(n1 = table(Event1$`Fork Length (mm)` >= 400),
         n2 = table(Event2$`Fork Length Corrected (mm)` >= 400),
         m2 = table(Event2_recaps$`Fork Length Corrected (mm)` >= 400))
seChapman(n1 = table(Event1$`Fork Length (mm)` >= 400),
         n2 = table(Event2$`Fork Length Corrected (mm)` >= 400),
         m2 = table(Event2_recaps$`Fork Length Corrected (mm)` >= 400))

Nstrat(n1 = table(Event1$`Fork Length (mm)` >= 400),
       n2 = table(Event2$`Fork Length Corrected (mm)` >= 400),
       m2 = table(Event2_recaps$`Fork Length Corrected (mm)` >= 400),
       estimator="Chapman")
sestrat(n1 = table(Event1$`Fork Length (mm)` >= 400),
       n2 = table(Event2$`Fork Length Corrected (mm)` >= 400),
       m2 = table(Event2_recaps$`Fork Length Corrected (mm)` >= 400),
       estimator="Chapman")

## I get slightly different abundance estimates: this is because in tab Chi Square and Abundance,
## the cells counting n1, n2 corrected, and m2 corrected have a "<" in the formula and the upper bound
## is set to 399, so fish measuring exactly 399mm are omitted.



