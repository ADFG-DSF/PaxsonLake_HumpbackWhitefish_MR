library(tidyverse)
library(recapr)
library(dsftools)  # devtools::install_github("ADFG-DSF/dsftools")

source("FDS_2024/R/recapr_prep.R")


### read data
Event1 <- read_csv("FDS_2024/flat_data/Event1.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`)) %>%
  filter(`Fork Length (mm)` >= 345)

Event2 <- read_csv("FDS_2024/flat_data/Event2.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`)) %>%
  filter(`Fork Length (mm)` >= 345)

# Event2_recaps <- subset(Event2, !is.na(`Tag Number`)) %>%
#   arrange(`Tag Number`)
# Event1_recaps <- subset(Event1, `Tag Number` %in% Event2_recaps$`Tag Number`) %>%
#   arrange(`Tag Number`)


allMR <- recapr_prep(ID="Tag Number", event1=Event1, event2=Event2, recap_codes="TL")

paired_recaps <- allMR$recaps$matched

# reproduce growth correction
# yreg <- paired_recaps$`Fork Length (mm)_event1`
# xreg <- paired_recaps$`Fork Length (mm)_event2`
lm1from2 <- with(paired_recaps, lm(`Fork Length (mm)_event1` ~ `Fork Length (mm)_event2`))
# lm1from2 <- lm(yreg ~ xreg)
summary(lm1from2)

Event2$fl_corr_mt <-
  predict(lm1from2, newdata = tibble(`Fork Length (mm)_event2`=
                                        Event2$`Fork Length (mm)`))
# Event2$fl_corr_mt <- 
#   predict(lm1from2, newdata = data.frame(xreg =
#                                         Event2$`Fork Length (mm)`))

# kind of want to impute actual event 1 lengths when they exist
for(i in 1:nrow(Event2)) {
  if(!is.na(Event2$`Tag Number`[i])) {
    if(Event2$`Tag Number`[i] %in% paired_recaps$`Tag Number_event1`) {
      Event2$fl_corr_mt[i] <- 
      paired_recaps$`Fork Length (mm)_event1`[paired_recaps$`Tag Number_event1`==Event2$`Tag Number`[i]]
    }
  }
}

Event1$length_stratum <- as.numeric(cut(Event1$`Fork Length (mm)`, 
                                        breaks=c(300,400,1000), right=FALSE)) 
Event2$length_stratum <- as.numeric(cut(Event2$fl_corr_mt, 
                                        breaks=c(300,400,1000), right=FALSE)) 


allMR2 <- recapr_prep(ID="Tag Number", event1=Event1, event2=Event2, recap_codes="TL")

# event 1 marks vs recaps
ks.test(allMR2$input_data$event1$`Fork Length (mm)`,
        allMR2$recaps$all$event1$`Fork Length (mm)`)

# event 2 caps vs recaps
ks.test(allMR2$input_data$event2$`Fork Length (mm)`,
        allMR2$recaps$all$event2$`Fork Length (mm)`)
ks.test(allMR2$input_data$event2$fl_corr_mt,
        allMR2$recaps$all$event2$fl_corr_mt)

## smalls
# event 1 marks vs recaps
ks.test(allMR2$input_data$event1$`Fork Length (mm)`[allMR2$input_data$event1$length_stratum==1],
        allMR2$recaps$all$event1$`Fork Length (mm)`[allMR2$recaps$all$event1$length_stratum==1])

# event 2 caps vs recaps
ks.test(allMR2$input_data$event2$`Fork Length (mm)`[allMR2$input_data$event2$length_stratum==1],
        allMR2$recaps$all$event2$`Fork Length (mm)`[allMR2$recaps$all$event2$length_stratum==1])
ks.test(allMR2$input_data$event2$fl_corr_mt[allMR2$input_data$event2$length_stratum==1],
        allMR2$recaps$all$event2$fl_corr_mt[allMR2$recaps$all$event2$length_stratum==1])

## bigs
# event 1 marks vs recaps
ks.test(allMR2$input_data$event1$`Fork Length (mm)`[allMR2$input_data$event1$length_stratum==2],
        allMR2$recaps$all$event1$`Fork Length (mm)`[allMR2$recaps$all$event1$length_stratum==2])

# event 2 caps vs recaps
ks.test(allMR2$input_data$event2$`Fork Length (mm)`[allMR2$input_data$event2$length_stratum==2],
        allMR2$recaps$all$event2$`Fork Length (mm)`[allMR2$recaps$all$event2$length_stratum==2])
ks.test(allMR2$input_data$event2$fl_corr_mt[allMR2$input_data$event2$length_stratum==2],
        allMR2$recaps$all$event2$fl_corr_mt[allMR2$recaps$all$event2$length_stratum==2])
