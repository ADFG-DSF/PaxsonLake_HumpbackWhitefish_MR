orderlikeanumber <- function(x, stopiferror=TRUE) {
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  if(any(x_num < 0, na.rm=TRUE)) {
    if(stopiferror) {
      stop("All numeric values must be positive")
    } else {
      return(order(x))
    }
  }
  
  num_digits <- floor(log10(x_num[!is.na(x_num)])) + 1
  num_digits[num_digits < 1] <- 1
  if(length(num_digits) > 0) {
    num_zeroes <- max(num_digits, na.rm=TRUE) - num_digits
  } else {
    num_zeroes <- NULL
  }
  
  zeroes <- sapply(num_zeroes, \(xx) paste0(rep(0, times=xx), collapse=""))
  
  x1 <- rep(NA, length(x))
  x1[is.na(x_num)] <- x[is.na(x_num)]
  x1[!is.na(x_num)] <- paste0(zeroes, x[!is.na(x_num)])
  
  return(order(x1))
}
reorderlikeanumber <- function(x, ...) x[orderlikeanumber(x, ...=...)]


recapr_prep <- function(ID, event=NULL, recap_codes=NULL, ...) {
  dots <- list(...)
  if(length(dots) == 0) stop("Need to add input data")
  
  # check to make sure all the data objects are data.frame or similar
  if(all(sapply(dots, inherits, c("data.frame", "matrix")))) {
    for(idots in 1:length(dots)) {
      if(!inherits(dots[[idots]], "data.frame")) {
        dots[[idots]] <- as.data.frame(dots[[idots]])
      }
    }
    if(length(dots) > 2) {
      stop("More than two data tables detected")
    }
    if(length(dots) == 2) {
      if(!is.null(event)) {
        warning("Two data tables detected, argument event= will be ignored")
      }
      if(!ID[1] %in% colnames(dots[[1]])) {
        stop("Specified ID= column not detected in the first data table")
      }
      if(!ID[length(ID)] %in% colnames(dots[[2]])) {
        stop("Specified ID= column not detected in the second data table")
      }
      out <- list()
      out$input_data <- dots
      the_events <- names(out$input_data)
    }
    if(length(dots) == 1) {
      if(!event %in% colnames(dots[[1]])) {
        stop("Specified event= column not detected in the data table")
      }
      if(length(ID) > 1) {
        stop("Only one ID= column may be specified if a single data table is used")
      }
      if(!ID %in% colnames(dots[[1]])) {
        stop("Specified ID= column not detected in the data table")
      }
      events_vec <- dots[[1]][[event]]
      the_events <- unique(events_vec)
      if(length(the_events) < 2) stop("Fewer than two events detected")
      if(length(the_events) > 2) stop("More than two events detected")
      
      out <- list()
      out$input_data <- list(subset(dots[[1]], events_vec==the_events[1]),
                  subset(dots[[1]], events_vec==the_events[2]))
      names(out$input_data) <- the_events
    }
  } else {
    stop("All data inputs must be data.frames or similar")
  }
  
  # this only works if ID is column names
  if(!length(ID) %in% 1:2) stop("Argument ID= may only have one or two elements")
  
  # full vectors of tag ID for both events
  recap_tags1 <- out$input_data[[1]][[ID[1]]]
  recap_tags2 <- out$input_data[[2]][[ID[length(ID)]]]
  
  # tabulate to see if there are multiple records for individuals
  t1 <- table(recap_tags1, useNA="no")
  t2 <- table(recap_tags2, useNA="no")
  problems1 <- t1[t1>1]
  problems2 <- t2[t2>1]
  if(length(problems1) > 0) {
    warning(c("Multiple records exist for individuals in event ", the_events[1], ": ",
           paste0(reorderlikeanumber(names(problems1), stopiferror = FALSE), 
           c(rep(", ", length(problems1)-1), ""))))
  }
  if(length(problems2) > 0) {
    warning(c("Multiple records exist for individuals in event ", the_events[2], ": ",
           paste0(reorderlikeanumber(names(problems2), stopiferror = FALSE), 
                  c(rep(", ", length(problems2)-1), ""))))
  }
  
  # a vector of JUST ID for recaptured individuals (duplicates are ignored)
  recaps_vec <- base::intersect(recap_tags1, recap_tags2)
  recaps_vec <- c(recaps_vec[!is.na(recaps_vec)], recap_codes)
  
  # data tables of recaps for both events
  recaps1 <- out$input_data[[1]][recap_tags1 %in% recaps_vec, ]
  recaps2 <- out$input_data[[2]][recap_tags2 %in% recaps_vec, ]
  
  recaps1 <- recaps1[orderlikeanumber(recaps1[[ID[1]]], stopiferror = FALSE), ]
  recaps2 <- recaps2[orderlikeanumber(recaps2[[ID[length(ID)]]], stopiferror = FALSE), ]
  
  # tabulate recaps for both events to check for multiple entries per individual
  t1 <- table(recaps1[[ID[1]]], useNA="no")
  t2 <- table(recaps2[[ID[length(ID)]]], useNA="no")
  
  # separate matched and non-matched subsets of recaps1 and recaps2
  recaps_vec_matched <- base::intersect(names(t1[t1==1]), names(t2[t2==1]))
  recaps1_matched <- recaps1[recaps1[[ID[1]]] %in% recaps_vec_matched, ]
  recaps2_matched <- recaps2[recaps2[[ID[length(ID)]]] %in% recaps_vec_matched, ]
  recaps1_unmatched <- recaps1[!recaps1[[ID[1]]] %in% recaps_vec_matched, ]
  recaps2_unmatched <- recaps2[!recaps2[[ID[length(ID)]]] %in% recaps_vec_matched, ]
  
  # throw a warning if there are unmatched individuals in recaps 
  if(nrow(recaps1_unmatched) > 0) {
    problems <- unique(recaps1_unmatched[[ID[1]]])   # maybe not unique
    warning(c("Unmatched records exist for recaptured individuals in event ", the_events[1], ": ",
              paste0(problems, c(rep(", ", length(problems)-1), ""))))
  }
  if(nrow(recaps2_unmatched) > 0) {
    problems <- unique(recaps2_unmatched[[ID[length(ID)]]])   # maybe not unique
    warning(c("Unmatched records exist for recaptured individuals in event ", the_events[2], ": ",
              paste0(problems, c(rep(", ", length(problems)-1), ""))))
  }
    
  # interleave matched
  # names(recaps1_matched) <- paste(names(recaps1), names(out)[1], sep="_")
  # names(recaps2_matched) <- paste(names(recaps2), names(out)[2], sep="_")
  # recaps_matched <- cbind(recaps1_matched, recaps2_matched)
  recaps_matched <- interleave(recaps1_matched, recaps2_matched, thenames=the_events)
  
  # # initialize sub-list
  out$recaps <- list()
  
  out$recaps$matched <- recaps_matched # [, order(names(recaps_matched))]
  out$recaps$unmatched <- list(recaps1_unmatched, recaps2_unmatched)
  names(out$recaps$unmatched) <- the_events
  
  out$recaps$all <- list(recaps1, recaps2)
  names(out$recaps$all) <- the_events
  
  # if((length(recaps_vec)==nrow(recaps1)) & (length(recaps_vec)==nrow(recaps2))) {
  #   # interleave columns
  #   recaps <- cbind(recaps1, recaps2)
  #   out$recaps <- recaps[, order(names(recaps))]
  # } else {
  #   # print(recaps1[[ID[1]]])
  #   # print(t2)
  #   problems1 <- t1[t1>1]
  #   problems2 <- t2[t2>1]
  #   if(length(problems1) > 0) {
  #     warning(c("Multiple records exist for RECAPTURED individuals in event ", names(out)[1], ": ",
  #               paste0(names(problems1), c(rep(", ", length(problems1)-1), ""))))
  #   }
  #   if(length(problems2) > 0) {
  #     warning(c("Multiple records exist for RECAPTURED individuals in event ", names(out)[2], ": ",
  #               paste0(names(problems2), c(rep(", ", length(problems2)-1), ""))))
  #   }
  # 
  #   out$recaps <- list(recaps1, recaps2)
  #   names(out$recaps) <- names(out)[1:2]
  # }
  
  
  return(out)
}

# propmatch <- function(x1,x2) {
#   maxn <- max(nchar(x1), nchar(x2))
#   minn <- min(nchar(x1), nchar(x2))
#   return(sum(unlist(strsplit(substr(x1, 1, minn), split="")) ==
#         unlist(strsplit(substr(x2, 1, minn), split="")))/maxn)
# }

interleave <- function(x1, x2, thenames=NULL) {
  
  names1 <- colnames(x1)
  names2 <- colnames(x2)
  
  if(is.null(thenames)) thenames <- 1:2
  colnames(x1) <- paste(colnames(x1), thenames[1], sep="_")
  colnames(x2) <- paste(colnames(x2), thenames[2], sep="_")
  
  if(length(intersect(names1, names2)) > 0) {
    colnum <- 1
    unused1 <- rep(TRUE, ncol(x1))
    unused2 <- rep(TRUE, ncol(x2))
    for(ix1 in 1:ncol(x1)) {
      if(any(names1[ix1] %in% names2)) {
        if(colnum == 1) {
          outdf <- as.data.frame(x1[ix1])
        } else {
          outdf[colnum] <- x1[ix1]   ######
        }
        outdf[colnum+1] <- x2[which(names2==names1[ix1])]   ######
        unused1[ix1] <- FALSE
        unused2[which(names2==names1[ix1])] <- FALSE
        
        # print(ix1)
        # print(outdf)
        colnum <- colnum+2
      }
    }
    unused <- cbind(x1[unused1], x2[unused2])#)as.data.frame(
    unused <- unused[, order(names(unused))]
    out <- cbind(outdf, unused)
  } else {
    unused <- cbind(x1, x2)
    out <- unused[order(names(unused))]   ######
  }
  
  return(out)
  
  
  # propmat <- outer(names1, names2, FUN=Vectorize(propmatch)) # rows from names1, columns from names2
  # thedim <- ifelse(length(names1) <= length(names2), 1, 2)
  # 
  # whichmaxes <- apply(propmat, thedim, which.max)
  
  # if(method=="first" | (method=="smaller" & ncol(x1) <= ncol(x2)) | (method=="larger" & ncol(x1) > ncol(x2))) {
  #   exactmatch <- sapply(sapply(names1, \(x) which(names2==x)), \(x) ifelse(is.null(x), NA, x))
  #   
  # }
}
# interleave(x1=recaps1_matched[,1:2], x2=recaps2_matched[,2:3])
# interleave(x1=as.data.frame(recaps1_matched)[,1:2], 
           # x2=as.data.frame(recaps2_matched)[3])


# simplify recap warning - maybe IDs aren't needed
# need more test cases!!
# - length(ID)==2
# - duplicates in recaps
# - TRY THINGS THAT SHOULD THROW AN ERROR
# make out$recaps$all an rbindish thing (single appended data table) --- EACH ROW SHOULD BE A UNIQUE FISH
# - make new sub-function append() to go with interleave() ???  -- NO

# ! done ! make sure the name $recaps isn't problematic!! maybe add a $input_data
# ! done ! add an error message when ... is empty! AND IF WE GIVE UP ON VECTOR INPUT
# ! done ! make interleaving smarter - try to match the ordering for event 1 (!!!)
# ! done ! - maybe make interleave() function
# ! done ! coerce matrix input to data.frame
# ! not going to ! handle vector input???

# maybe rename as recapr_data, and write another function recapr_tabulate??
# - would need columns for stratum - ANY OTHERS??


# aa <- recapr_prep(ID="Tag Number", event1=Event1, event2=Event2, recap_codes="TL")
# str(aa)
# 
# bothevents <- rbind(select(Event1, c("Event", "Tag Number", "Fork Length (mm)")),
#                     select(Event2, c("Event", "Tag Number", "Fork Length (mm)"))) #%>% as.data.frame
# aa <- recapr_prep(ID="Tag Number", data=bothevents, event="Event", recap_codes="TL")
# str(aa)
# 
# aa <- recapr_prep(ID="Tag Number", data=(as.matrix(bothevents)), event="Event", recap_codes="TL")
# str(aa) 
# 
# 
# bothevents1 <- bothevents
# bothevents1$`Tag Number`[1] <- 770
# aa <- recapr_prep(ID="Tag Number", data=(as.matrix(bothevents1)), event="Event", recap_codes="TL")
# str(aa)
# 
# 
# Event11 <- Event1
# names(Event11)[7] <- "TagNumber"
# aa <- recapr_prep(ID=c("TagNumber","Tag Number"), event1=Event1, event2=Event2, recap_codes="TL")
# str(aa)
# 
# bothevents1 <- bothevents
# bothevents1$Event <- 3
# aa <- recapr_prep(ID=c("Tag Number","Tag Number","steve"),event1=Event1, event2=Event2)   # needed data=bothevents
# str(aa)

correct_growth <- function(x, 
                           event_keep, event_adjust,
                           column_keep, column_adjust,
                           ID_keep, ID_adjust) {
  # insert error checking
  
  # make a copy to modify
  x1 <- x
  
  # regression bit
  yreg <- x$recaps$matched[[paste(column_keep, event_keep, sep="_")]]
  xreg <- x$recaps$matched[[paste(column_adjust, event_adjust, sep="_")]]
  lm1 <- lm(yreg ~ xreg)
  
  # predict from regression
  ypred <- predict(lm1, newdata = data.frame(xreg=x$input_data[[event_adjust]][[column_adjust]]))
  
  # fill in individuals as available
  ytag <- x$input_data[[event_adjust]][[ID_adjust]]
  keeptag <- x$recaps$matched[[paste(ID_keep, event_keep, sep="_")]]
  for(iy in seq_along(ytag)) {
    if(!is.na(ytag[iy])) {
      if(ytag[iy] %in% keeptag) {
        ypred[iy] <- yreg[keeptag==ytag[iy]]
      }
    }
  }
  
  ## actually should make new columns for these:
  # need to change it in x1$input_data[[event_adjust]][[column_adjust]]
  x1$input_data[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]] <- unname(ypred)
  
  # change it in matched
  x1$recaps$matched[[paste(column_adjust, event_adjust, "adjusted", sep="_")]] <- yreg
  
  # need to change it in x1$recaps$unmatched[[event_adjust]][[column_adjust]]
  x1$recaps$unmatched[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]] <-
    unname(predict(lm1, newdata = data.frame(xreg=x1$recaps$unmatched[[event_adjust]][[column_adjust]])))
  
  # need to change it in x1$recaps$all[[event_adjust]][[column_adjust]]
  x1$recaps$all[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]] <-
    unname(predict(lm1, newdata = data.frame(xreg=x1$recaps$all[[event_adjust]][[column_adjust]])))
  ytag <- x1$recaps$all[[event_adjust]][[ID_adjust]]
  keeptag <- x$recaps$matched[[paste(ID_keep, event_keep, sep="_")]]
  for(iy in seq_along(ytag)) {
    if(!is.na(ytag[iy])) {
      if(ytag[iy] %in% keeptag) {
        x1$recaps$all[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]][iy] <- 
          yreg[keeptag==ytag[iy]]  
      }
    }
  }
  
  return(x1)
}
aa1 <- correct_growth(x=aa, 
               event_keep="event1", 
               event_adjust="event2", 
               column_keep="Fork Length (mm)", 
               column_adjust="Fork Length (mm)",
               ID_keep="Tag Number",
               ID_adjust="Tag Number")

lm1 <- lm(aa$recaps$matched$`Fork Length (mm)_event1` ~ aa$recaps$matched$`Fork Length (mm)_event2`)

plot(aa1$input_data$event2$`Fork Length (mm)`, aa1$input_data$event2$`Fork Length (mm)_adjusted`,
     pch=ifelse(aa1$input_data$event2$`Tag Number` %in% aa$recaps$matched$`Tag Number_event1`, 16, 1))
abline(lm1)
abline(0, 1, lty=3)

plot(aa1$recaps$matched$`Fork Length (mm)_event2`, aa1$recaps$matched$`Fork Length (mm)_event2_adjusted`)
plot(aa1$recaps$matched$`Fork Length (mm)_event1`, aa1$recaps$matched$`Fork Length (mm)_event2_adjusted`)
abline(0, 1, lty=3)

plot(aa1$recaps$unmatched$event2$`Fork Length (mm)`, aa1$recaps$unmatched$event2$`Fork Length (mm)_adjusted`)
abline(lm1)

plot(aa1$recaps$all$event2$`Fork Length (mm)`, aa1$recaps$all$event2$`Fork Length (mm)_adjusted`,
     pch=ifelse(aa1$recaps$all$event2$`Tag Number` %in% aa$recaps$matched$`Tag Number_event1`, 16, 1))
abline(lm1)

all.equal(aa1$input_data$event1, aa$input_data$event1)
all.equal(aa1$recaps$unmatched$event1, aa$recaps$unmatched$event1)
all.equal(aa1$recaps$all$event1, aa$recaps$all$event1)
