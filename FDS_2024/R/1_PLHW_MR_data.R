## The primary purpose of this script is to validate the Excel-based analysis
## conducted by Corey Schwanke and April Behr.
##
## There are no analysis outputs produced in this script that are directly used,
## however, the results should match (at least inferentially) what is reported.
##
## The Excel-based analysis consisted of a few iterations, as some discrepancies
## were encountered.
## * Growth was detected between the two capture events, and the second-event lengths
##   were adjusted to be comparable to the first.  Originally, a linear regression
##   was used to model the second-event lengths as a function of the first-event
##   lengths, and the regression formula was inverted to "correct" the second-
##   event lengths.  This was changed to a regression of the first-event lengths
##   as a function of the second, which could be directly applied to correct the 
##   second-event lengths.
## * Counts of fish per size stratum did not match exactly; this was a result of
##   the Excel formula missing a few fish with corrected lengths falling between
##   the 345-399mm bin and the 400+mm bin.  The Excel sheet was modified to 
##   calculate the FLOOR() of the corrected lengths.
##
## KS and Chi^2 tests all result in the same inferences as the Excel sheet, though
## there are some differences in the exact numbers.  It is likely that the KS 
## tests are calculated slightly differently in R vs Excel.




### load packages
library(tidyverse)
library(recapr)
library(dsftools)  # devtools::install_github("ADFG-DSF/dsftools")


### read data
Event1 <- Event1_raw <- read_csv("FDS_2024/flat_data/Event1.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`))

Event2 <- Event2_raw <- read_csv("FDS_2024/flat_data/Event2.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`))

Event2_recaps <- Event2_recaps_raw <- subset(Event2, !is.na(`Tag Number`)) %>%
  arrange(`Tag Number`)
Event1_recaps <- Event1_recaps_raw <- subset(Event1, `Tag Number` %in% Event2_recaps$`Tag Number`) %>%
  arrange(`Tag Number`)



### reproduce KS tests
# reproduce growth correction
xreg <- Event2_recaps$`Fork Length (mm)`[Event2_recaps$`Tag Number` != "TL"]
yreg <- Event1_recaps$`Fork Length (mm)`
lm_growth <- lm(yreg ~ xreg)
summary(lm_growth)
lm_predict <- predict(lm_growth, newdata = data.frame(xreg=Event2$`Fork Length (mm)`))

t.test(yreg-xreg)


####### addition 11/22
# impute first event lengths for recaps
for(i in 1:nrow(Event2)) {
  if(!is.na(Event2$`Tag Number`[i])) {
    if(Event2$`Tag Number`[i] %in% Event1_recaps$`Tag Number`) {
      lm_predict[i] <- Event1_recaps$`Fork Length (mm)`[Event1_recaps$`Tag Number` == Event2$`Tag Number`[i]]
    }
  }
}

plot(Event2$`Fork Length Corrected (mm)`, lm_predict)
abline(0,1)

####### addition 11/22
Event2$`Fork Length Corrected (mm)` <- lm_predict

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

## creating a function to streamline display of KS test output
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
       # Event1_recaps$`Fork Length (mm)`,
       Event2_recaps$`Fork Length Corrected (mm)`,
       xlab="mm FL", main="Event 1",
       col=c(1,2), legend=c("All","Recaps"))
ksplot(Event2$`Fork Length Corrected (mm)`, 
        Event2_recaps$`Fork Length Corrected (mm)`,
       xlab="mm FL", main="Event 2",
       col=c(1,4), legend=c("All","Recaps"))

ksplot(Event1$`Fork Length (mm)` %s_l% 400, 
       # Event1_recaps$`Fork Length (mm)` %s_l% 400,
       Event2_recaps$`Fork Length Corrected (mm)` %s_l% 400,
       xlab="mm FL", main="Event 1 - Small (<400mm)",
       col=c(1,2), legend=c("All","Recaps"))
ksplot(Event2$`Fork Length Corrected (mm)` %s_l% 400, 
        Event2_recaps$`Fork Length Corrected (mm)` %s_l% 400,
       xlab="mm FL", main="Event 2 - Small (<400mm)",
       col=c(1,4), legend=c("All","Recaps"))

ksplot(Event1$`Fork Length (mm)` %s_geq% 400, 
       # Event1_recaps$`Fork Length (mm)` %s_geq% 400,
       Event2_recaps$`Fork Length Corrected (mm)` %s_geq% 400,
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
## seems to match (at least for the "small" stratum)
####### though does not match tallies as of 11/22



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

## tallies do not match 11/22 (missing 399-400cm corrected lengths)

## 5/20/25 verify cap probabilities
# first event
table(Event2_recaps$`Fork Length Corrected (mm)` >= 400) /
  table(Event2$`Fork Length Corrected (mm)` >= 400)

# second event
table(Event2_recaps$`Fork Length Corrected (mm)` >= 400) /
  table(Event1$`Fork Length (mm)` >= 400)


## update 5/19/25 to verify stratified proportions
lengthinput <- c(Event1$`Fork Length (mm)`,
                 Event2$`Fork Length Corrected (mm)`)
lengthstrata <- cut(lengthinput, breaks=c(0,400, 1000), right=FALSE)
lengthbin <- cut(lengthinput, breaks=c(seq(345, 515, by=10)), right=FALSE)

Nhat <- NChapman(n1 = table(Event1$`Fork Length (mm)` >= 400),
                 n2 = table(Event2$`Fork Length Corrected (mm)` >= 400),
                 m2 = table(Event2_recaps$`Fork Length Corrected (mm)` >= 400))
seNhat <- seChapman(n1 = table(Event1$`Fork Length (mm)` >= 400),
          n2 = table(Event2$`Fork Length Corrected (mm)` >= 400),
          m2 = table(Event2_recaps$`Fork Length Corrected (mm)` >= 400))

the_tbl <- ASL_table(age=lengthbin,
          stratum=as.numeric(lengthstrata),
          Nhat=as.numeric(Nhat),
          se_Nhat=as.numeric(seNhat))

par(mfrow=c(1,1))
plot(the_tbl$phat)
segments(x0=1:nrow(the_tbl),
         y0=the_tbl$phat - qnorm(0.95)*the_tbl$se_phat,
         y1=the_tbl$phat + qnorm(0.95)*the_tbl$se_phat)



# As of 11/25, April attempted to reconcile the methods by ROUNDing the corrected lengths
# This is an investigation of whether that worked
# IT DID NOT FULLY WORK, recommended switching to FLOOR for full consistency
lc <- data.frame(
          V1 = c(339,331,362,355,360,357,364,355,351,
                 380,393,352,344,340,390,408,393,324,408,372,359,
                 357,433,400,369,388,391,457,424,376,372,393,394,415,
                 383,380,417,371,366,395,385,418,402,363,377,342,
                 441,379,365,391,373,370,365,359,375,379,404,376,344,
                 364,400,364,381,375,453,361,339,346,365,381,405,331,
                 381,353,360,355,335,377,344,375,405,399,382,373,
                 378,364,376,367,366,359,388,408,381,437,374,400,373,
                 384,411,362,385,383,405,375,384,391,365,390,366,
                 380,389,399,375,407,366,350,378,400,350,385,351,329,
                 353,338,368,361,386,358,347,343,369,377,354,384,420,
                 380,364,412,431,413,434,393,388,395,397,379,367,
                 373,360,355,383,412,376,382,409,435,371,366,405,371,
                 377,344,400,420,368,399,415,373,394,393,446,389,
                 385,357,393,378,385,395,363,371,380,364,378,389,381,
                 372,409,366,372,359,355,369,365,380,350,354,342,340,
                 385,383,372,339,386,389,372,342,334,365,350,368,
                 383,376,407,366,424,391,394,450,356,337,368,342,369,
                 373,371,404,398,340,346,337,335,347,355,329,412,
                 370,339,360,359,331,329,342,415,408,351,365,372,356,
                 362,342,359,340,345,340,336,346,368,364,359,366,352,
                 381,364,345,352,366,346,361,386,385,366,347,338,
                 340,346,342,353,359,358,390,373,347,390,422,373,389,
                 385,362,356,342,368,362,360,375,355,353,415,262,
                 364,371,359,341,347,353,402,339,389,359,341,340,388,
                 368,391,395,398,374,385,377,364,376,376,317,347,408,
                 352,355,353,355,343,328,356,348,410,354,339,387,
                 360,346,328,325,333,351,342,381,365,325,343,371,350,
                 345,341,346,395,364,351,348,399,347,345,356,342,
                 350,347,341,351,365,348,322,359,341,340,345,353,365,
                 371,357,341,327,346,348,347,375,370,358,356,363,371,
                 375,425,405,346,375,349,350,346,369,346,351,408,
                 370,374,367,361,376,373,407,394,445,369,392,404,414,
                 358,390,404,374,378,413,365,386,402,392,366,368,
                 413,371,395,359,366,382,373,374,350,370,432,381,362,
                 381,362,375,361,443,353,355,371,360,352,402,357,361,
                 365,351,379,370,345,328,325,358,346,338,346,364,
                 339,340,332,395,356,361,362,364,333,358,346,410,343,
                 346,408,342,341,313,347,358,399,349,346,384,346,
                 342,354,335,395,366,356,384,376,317,336,412,392,375,
                 356,370,352,351,347,355,340,354,360,323,340,341,355,
                 356,338,347,346,452,355,404,335,345,356,354,361,
                 343,364,395,397,369,439,362,368,380,361,376,364,353,
                 342,384,356,398,400,438,367,373,413,383,412,388,
                 375,382,385,365,396,377,385,424,365,405,417,365,365,
                 414,424,417,371,395,389,364,416,383,426,354,384,375,
                 412,365,375,395,424,375,380,371,352,370,375,438,
                 383,402,412,375,361,375,395,384,387,363,365,370,373,
                 390,408,385,375,383,385,412,370,364,391,387,385,
                 380,365,369,365,367,370,364,415,386,383,414,370,395,
                 386,360,382,383,380,395,384,375,356,393,392,375,386,
                 385,355,392,336,383,412,381,395,365,338,356,349,
                 354,335,364,356,361,365,360,334,362,346,371,363,358,
                 353,367,400,336,346,385,346,356,372,353,346,346,
                 347,358,351,338,360,356,338,362,347,361,376,361,364,
                 356,356,358,328,350,338,355,340,344,360,351,348,366,
                 371,347,325,353,380,336,350,372,371,376,348,327,
                 328,352,430,347,350,327,328,315,336,357,336,350,353,
                 389,375,424,356,385,375,383,353,404,443,365,383,
                 414,377,361,353,398,360,402,375,400,404,374,375,358,
                 384,412,364,390,373,373,422,395,381,395,371,384,404,
                 361,380,419,395,381,365,376,369,409,366,374,395,
                 373,375,412,377,354,364,385,379,369,439,360,367,395,
                 415,458,367,385,382,362,406,402,378,465,389,375,
                 373,383,390,373,368,405,366,384,369,379,435,389,393,
                 379,359,373,418,361,364,366,385,412,384,394,383,373,
                 464,360,428,418,412,392,365,422,389,375,369,378,
                 386,422,496,386,462,374,442,375,358,389,380,465,379,
                 413,407,362,391,387,384,390,399,453,442,393,379,
                 374,411,383,364,375,391,365,385,416,373,399,384,371,
                 396,382,393,375,379,359,394,376,365,379,388,426,372,
                 365,447,389,388,357,370,373,383,380,356,369,386,
                 377,420,395,391,385,389,372,397,374,408,398,373,381,
                 406,390,397,393,395,379,394,369,372,382,379,410,
                 380,385,385,392,356,386,408,375,412,400,387,381,375,
                 369,380,400,423,388,381,422,366,375,442,422,372,368,
                 390,418,376,388,382,414,367,400,376,426,425,452,
                 372,377,383,399,362,408,400,400,346,382,391,384,364,
                 423,397,375,405,371,380,395,385,384,387,385,373,
                 410,384,399,382,388,435,401,404,399,400,387,384,394,
                 339,391,409,368,373,390,399,423,394,379,370,396,405,
                 372,382,375,420,368,404,386,396,375,388,443,367,
                 408,391,388,369,415,373,387,390,404,393,392,413,389,
                 496,397,374,383,388,400,402,398,370,393,372,396,
                 368,418,430,404,382,378,420,394,387,373,408,420,411,
                 394,415,384,417,412,419,400,425,358,352,377,374,410,
                 402,367,410,379,368,415,398,404,386,372,390,374,
                 375,408,390,405,408,403,416,394,367,381,380,378,406,
                 375,379,434,399,389,368,397,391,430,374,420,386,
                 362,377,390,386,405,374,372,383,381,468,435,399,400,
                 419,406,390,381,447,382,461,428,407,379,396,373,377,
                 398,405,356,380,415,384,387,381,382,413,382,386,
                 408,413,414,389,404,395,423,364,397,405,389,370,412,
                 398,378,420,419,383,454,361,399,366,442,405,428,
                 391,404,435,472,376,391,403,426,391,360,412,455,408,
                 412,371,437,436,377,429,353,372,455,387,401,423,396,
                 373,398,374,409,426,386,408,368,400,404,384,359,
                 355,328,340,332,328,339,331,340,316,343,354,358,347,
                 305,342,354,349,360,345,351,372,337,320,342,344,
                 345,360,356,315,358,376,346,336,339,360,351,354,335,
                 333,354,340,336,335,321,323,340,352,325,335,338,335,
                 333,343,365,325,340,352,350,337,335,364,321,356,
                 348,349,357,368,358,327,345,330,425,329,327,325,321,
                 352,340,344,348,322,329,345,356,340,344,345,333,
                 338,331,376,374,361,342,359,334,326,393,341,355,339,
                 341,344,335,334,368,379,386,344,336,357,334,367,332,
                 372,365,339,364,318,330,358,355,353,369,345,354,
                 323,356,332,326,316,348,308,352,307,347,383,320,346,
                 356,334,340,344,372,348,351,335,347,318,330,329,
                 354,391,352,343,336,369,330,375,335,340,345,343,359,
                 357,337,355,360,329,353,358,337,333,323,370,358,327,
                 337,328,325,369,354,343,364,358,332,350,324,360,
                 337,341,359,344,352,356,348,357,325,365,340,369,359,
                 371,340,363,350,336,338,347,352,356,369,349,359,
                 369,361,376,371,374,379,356,366,336,342,348,337,350,
                 359,338,330,337,353,335,384,323,430,360,360,362,370,
                 354,321,335,351,344,355,330,372,337,341,351,365,
                 326,340,344,356,376,360,338,332,330,356,329,385,363,
                 336,367,345,355,331,353,350,346,341,345,328,344,
                 365,362,347,342,341,360,373,349,352,353,323,356,347,
                 252,368,348,346,338,362,353,286,347,340,445,373,344,
                 345,354,362,346,382,337,349,336,375,362,366,353,
                 344,368,369,357,347,336,403,347,334,345,358,328,348,
                 360,353,344,362,334,365,364,368,367,332,335,333,
                 368,346,351,343,376,329,362,345,331,380,343,328,329,
                 354,387,367,341,362,328,396,378,338,335,362,330,365,
                 359,339,325,343,399,356,376,367,410,340,336,344,
                 347,352,356,354,339,337,341,316,353,348,368,346,340,
                 368,356,375,330,342,328,361,337,348,330,326,353,
                 298,331,334,356,350,461,371,353,356,318,334,346,393,
                 348,333,342,347,345,357,344,358,320,340,368,373,337,
                 351,327,315,323,355,354,339,375,373,388,341,340,
                 363,430,336,347,337,434,344,358,344,340,326,343,373,
                 405,313,369,354,352,350,378,372,332,375,335,332,
                 339,354,342,361,359,364,328,320,346,361,372,357,323,
                 358,365,380,385,395,417,418,419,441,398,364,416,387,
                 388,396,361,381,364,393,377,420,366,388,427,454,
                 385,367,362,363,395,367,346,414,373,400,368,360,409,
                 364,359,403,416,341,362,354,374,394,347,360,352,
                 394,363,284,362,358,366,343,386,340,346,351,359,362,
                 378,357,346,354,338,349,365,329,352,328,343,337,358,
                 350,366,317,352,345,318,332,328,347,343,325,348,
                 353,347,357,342,364,338,366,316,329,348,357,328,366,
                 342,335,339,334,327,363,327,340,376,381,346,318,
                 330,350,343,366,369,385,383,388,367,386,420,363,359,
                 377,386,383,372,381,370,368,387,391,389,372,382,356,
                 438,377,367,371,400,378,392,371,367,369,362,380,
                 386,359,388,399,384,397,363,360,372,368,378,375,373,
                 361,351,402,372,360,374,360,350,331,323,378,450,
                 358,368,342,331,353,329,340,378,358,323,360,328,315,
                 336,373,357,370,332,339,348,337,334,370,331,333,340,
                 403,347,378,413,359,380,399,383,371,435,358,412,
                 392,467,405,400,390,383,381,343,356,386,398,379,406,
                 370,404,376,348,357,364,364,373,329,333,328,325,
                 355,347,350,346,349,351,355,302,374,336,334,352,327,
                 425,374,400,402,370,373,388,383,393,380,434,359,379,
                 396,386,390,392,405,389,379,370,433,383,412,373,
                 389,363,369,401,366,379,435,350,354,359,352,326,339,
                 360,325,327,368,332,350,353,345,367,372,339,346,
                 357,396,350,352,331,351,373,356,350,353,346,406,402,
                 432,375,395,429,379,386,372,362,391,405,389,404,432,
                 396,365,381,411,379,389,374,356,381,348,330,392,
                 348,359,337,328,331,337,359,352,324,337,400,357,367,
                 359,364,372,356,370,415,378,405,385,420,381,378,
                 369,433,473,356,412,386,352,372,349,338,359,372,321,
                 334,357,350,326,359,327,347,342,372,333,370,349,350,
                 377,366,368,364,408,422,383,391,380,376,381,394,
                 392,387,396,415,347,374,359,333,349,362,350,336,361,
                 358,346,359,457,345,332,344,323,329,326,355,347,
                 346,349,333,374,392,420,367,361,372,354,375,364,386,
                 375,388,380,398,369,410,412,383,374,366,408,390,365,
                 425,351,368,380,382,372,365,337,397,355,444,303,
                 356,352,354,350,348,349,337,375,354,337,389,344,335,
                 352,347,334,325,332,362,347,361,358,336,354,336,
                 346,357,354,364,414,334,337,345,368,409,420,374,369,
                 394,400,359,391,372,430,435,374,449,365,376,366,359,
                 403,431,404,370,387,389,408,319,333,338,309,322,
                 340,328,344,339,357,343,329,433,341,324,325,387,328,
                 330,360,366,326,344,384,370,384,348,328,324,336,
                 339,397,343,336,325,342,337,329,335,362,340,324,314,
                 335,343,345,343,339,340,351,381,301,353,348,339,336,
                 339,312,346,334,321,322,369,332,330,323,347,372,
                 356,338,343,319,342,413,471,435,322,339,330,316,331,
                 374,336,329,338,360,338,336,346,349,338,334,337,
                 356,340,417,308,328,324,390,330,364,382,359,367,310,
                 364,364,328,364,319,350,413,320,326,318,347,344,342,
                 328,312,333,367,328,358,339,310,346,388,336,324,
                 336,360,339,327,361,353,343,345,339,350,340,328,321,
                 308,349,368,337,342,348,332,333,342,344,326,336,
                 329,337,355,343,356,353,345,323,336,348,325,364,302,
                 328,393,327,323,388,337,338,333,345,332,337,400,341,
                 329,304,340,328,314,326,353,347,328,319,331,309,
                 330,399,338,323,337,325,411,321,328,343,325,345,334,
                 367,348,345,378,342,318,320,315,335,339,249,350,
                 390,340,358,325,366,321,357,337,373,345,365,364,345,
                 342,332,343,338,329,379,351,374,347,343,316,364,371,
                 344,436,399,365,360,364,379,343,339,322,353,410,
                 347,345,348,364,343,337,362,336,339,346,366,364,352,
                 337,340,344,358,362,359,354,325,335,346,356,342,
                 351,327,333,347,268,337,364,369,349,376,333,245,336,
                 370,335,347,360,348,345,337,315,349,351,345,339,361,
                 391,341,349,377,345,399,350,366,405,365,342,393,
                 360,357,340,356,349,335,346,348,340,371,354,356,372,
                 344,376,346,342,349,339,356,342,348,345,363,344,
                 345,347,349,340,380,399,355,350,356,344,348,346,328,
                 342,349,361,348,459,362,354,342,324,343,386,346,352,
                 355,349,337,335,345,370,347,366,364,334,363,365,
                 366,365,360,351,341,376,375,354,345,340,350,351,368,
                 336,362,325,323,383,360,331,341,353,349,370,379,
                 400,360,368,334,354,326,356,352,358,337,350,332,354,
                 411,352,383,408,345,349,351,361,347,363,342,341,352,
                 322,356,321,327,356,338,393,364,438,359,361,346,
                 329,337,375,342,345,332,355,340,326,405,347,331,334,
                 384,349,354,332,395,400,379,364,366,370,336,384,
                 358,353,357,400,369,375,364,357,362,425,358,369,359,
                 471,427,354,379,344,356,362,370,355,359,372,344,368,
                 374,341,337,351,350,337,339,360,349,340,345,328,
                 352,348,367,356,362,341,354,343,357,408,364,369,332,
                 356,365,363,354,332,328,350,340,449,326,354,343,
                 405,347,346,357,335,333,340,344,395,379,349,374,352,
                 340,347,350,325,371,318,336,342,335,346,343,360,362,
                 345,351,344,339,342,335,347,369,338,358,343,390,
                 371,342,343,362,344,340,335,347,323,365,355,363,344,
                 370,356,336,348,330,373,358,349,326,336,362,369,
                 340,367,320,349,352,357,351,343,334,358,364,432,342,
                 366,362,417,339,345,371,394,346,328,383,343,331,333,
                 359,348,380,343,370,378,335,351,337,369,339,328,
                 355,347,354,343,365,344,349,377,340,357,345,342,358,
                 361,348,371,337,340,335,344,331,360,323,313,364,
                 339,342,328,336,322,413,344,328,333,348,320,337,373,
                 367,349,344,328,327,338,343,329,325,339,328,319,350,
                 328,339,305,325,331,339,343,329,330,326,352,319,
                 324,328,349,336,332,327,318,354,384,348,386,350,323,
                 347,308,314,343,335,346,325,334,325,331,335,337,
                 332,311,330,322,320,330,341,359,325,350,339,309,336,
                 347,333,316,313,322,339,318,348,327)
) %>% filter(V1>=345)

table(cut(lc$V1, breaks=c(345, 400, 1000), right=F))
table(Event2$`Fork Length Corrected (mm)` >= 400)

table(lc$V1 >= 400, Event2$`Fork Length Corrected (mm)` >= 400)
subset(data.frame(lc=lc$V1, E2=Event2$`Fork Length Corrected (mm)`),
       lc >= 400 & E2 < 400)
