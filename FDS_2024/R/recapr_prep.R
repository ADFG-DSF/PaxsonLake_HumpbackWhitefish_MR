# a helper function that returns the order of a character vector, but preserves
# numeric values correctly (i.e. 9 comes before 10)
orderlikeanumber <- function(x, stopiferror=TRUE) {
  
  # error handling: 
  # - check if input is a vector
  if(!(is.vector(x) & is.atomic(x))) stop("Input must be a vector")
  
  # - what happens with NA?  It puts them at the end by default.  We can live with that.
  
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

# a wrapper of orderlikeanumber that returns a re-ordered vector
reorderlikeanumber <- function(x, ...) x[orderlikeanumber(x, ...=...)]


## this is the big one: it reorganizes mark-recapture datasets
recapr_prep <- function(ID, event=NULL, recap_codes=NULL, ...) {
  dots <- list(...)
  if(length(dots) == 0) stop("Need to add input data")
  
  # don't know why anyone would use NA in recap_codes, but not allowing this
  if(any(is.na(recap_codes))) {
    warning("NA not supported in recap_codes, removing these")
    recap_codes <- recap_codes[!is.na(recap_codes)]
  }
  
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
      if(is.null(names(dots))) {
        names(dots) <- 1:2
      }
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
  recaps_vec_matched <- base::intersect(names(t1[t1==1]), names(t2[t2==1]))   # need to exclude recap_codes
  recaps_vec_matched <- recaps_vec_matched[!(recaps_vec_matched %in% recap_codes)]
  
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
  
  ### trying appending unmatched to matched   ------- see if i can make this NOT have to be a data.frame
  recaps1_unmatched_rbind1 <- recaps1_unmatched
  names(recaps1_unmatched_rbind1) <- paste(names(recaps1_unmatched), the_events[1], sep="_")
  recaps1_unmatched_rbind <- data.frame(row.names=rownames(recaps1_unmatched_rbind1))
  for(j in 1:ncol(recaps_matched)) {
    if(names(recaps_matched)[j] %in% names(recaps1_unmatched_rbind1)) {
      addthis <- recaps1_unmatched_rbind1[[names(recaps_matched)[j]]]
    } else {
      addthis <- rep(NA, nrow(recaps1_unmatched_rbind1))
    }
    recaps1_unmatched_rbind[[names(recaps_matched)[j]]] <- addthis
  }

  recaps2_unmatched_rbind2 <- recaps2_unmatched
  names(recaps2_unmatched_rbind2) <- paste(names(recaps2_unmatched), the_events[2], sep="_")
  recaps2_unmatched_rbind <- data.frame(row.names=rownames(recaps2_unmatched_rbind2))
  for(j in 1:ncol(recaps_matched)) {
    if(names(recaps_matched)[j] %in% names(recaps2_unmatched_rbind2)) {
      addthis <- recaps2_unmatched_rbind2[[names(recaps_matched)[j]]]
    } else {
      addthis <- rep(NA, nrow(recaps2_unmatched_rbind2))
    }
    recaps2_unmatched_rbind[[names(recaps_matched)[j]]] <- addthis
  }
  
  recaps_all <- rbind(recaps_matched, recaps1_unmatched_rbind, recaps2_unmatched_rbind)
  ###
  
  # # initialize sub-list
  out$recaps <- list()
  
  out$recaps$matched <- recaps_matched # [, order(names(recaps_matched))]
  out$recaps$unmatched <- list(recaps1_unmatched, recaps2_unmatched)
  names(out$recaps$unmatched) <- the_events
  
  ## previous version
  # out$recaps$all <- list(recaps1, recaps2)
  # names(out$recaps$all) <- the_events
  
  out$recaps$all <- recaps_all
  
  class(out) <- "MR_data"
  return(out)
}

# propmatch <- function(x1,x2) {
#   maxn <- max(nchar(x1), nchar(x2))
#   minn <- min(nchar(x1), nchar(x2))
#   return(sum(unlist(strsplit(substr(x1, 1, minn), split="")) ==
#         unlist(strsplit(substr(x2, 1, minn), split="")))/maxn)
# }



# a helper function that splices two data.frames by column name
interleave <- function(x1, x2, thenames=NULL) {
  
  # need error checking: x1 and x2 need to be data.frames or similar, and of the same length
  if(!inherits(x1, c("data.frame", "matrix")) | !inherits(x2, c("data.frame", "matrix"))) {
    stop("Inputs must be data.frames or similar")
  }
  if(nrow(x1) != nrow(x2)) stop("Inputs must have the same number of rows")
  
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



## the next big one - this takes an input from recapr_prep and uses linear regression to correct for growth
correct_growth <- function(x, 
                           event_keep, event_adjust,
                           column_keep, column_adjust=column_keep,
                           ID_keep, ID_adjust=ID_keep,
                           impute = TRUE) {
  # insert error checking
  # - needs to be an object returned by recapr_prep
  if(!inherits(x, "MR_data")) stop("Argument x= must be an object returned from recapr_prep()")
  
  # - event names must exist
  if(!(event_keep %in% names(x$input_data))) {
    stop("Argument event_keep= not found in event names")
  }
  if(!(event_adjust %in% names(x$input_data))) {
    stop("Argument event_adjust= not found in event names")
  }
  
  # - column names must exist in both events
  if(!(column_keep %in% names(x$input_data[[event_keep]]))) {
    stop(paste("Argument column_keep= not found in", event_keep, "data"))
  }
  if(!(column_adjust %in% names(x$input_data[[event_adjust]]))) {
    stop(paste("Argument column_adjust= not found in", event_adjust, "data"))
  }
  
  # - ID names must exist (maybe make ID_adjust=ID_keep by default)
  if(!(ID_keep %in% names(x$input_data[[event_keep]]))) {
    stop(paste("Argument ID_keep= not found in", event_keep, "data"))
  }
  if(!(ID_adjust %in% names(x$input_data[[event_adjust]]))) {
    stop(paste("Argument ID_adjust= not found in", event_adjust, "data"))
  }
  
  # make a copy to modify
  x1 <- x
  
  # regression bit
  yreg <- x$recaps$matched[[paste(column_keep, event_keep, sep="_")]]
  xreg <- x$recaps$matched[[paste(column_adjust, event_adjust, sep="_")]]
  lm1 <- lm(yreg ~ xreg)
  
  # predict from regression
  ypred <- predict(lm1, newdata = data.frame(xreg=x$input_data[[event_adjust]][[column_adjust]]))
  
  # fill in individuals as available
  if(impute) {
    ytag <- x$input_data[[event_adjust]][[ID_adjust]]
    keeptag <- x$recaps$matched[[paste(ID_keep, event_keep, sep="_")]]
    for(iy in seq_along(ytag)) {
      if(!is.na(ytag[iy])) {
        if(ytag[iy] %in% keeptag) {
          ypred[iy] <- yreg[keeptag==ytag[iy]]
        }
      }
    }
  }
  
  ## actually should make new columns for these:
  # need to change it in x1$input_data[[event_adjust]][[column_adjust]]
  x1$input_data[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]] <- unname(ypred)
  
  # change it in matched
  if(impute) {
    # x1$recaps$matched[[paste(column_adjust, event_adjust, "adjusted", sep="_")]] <- yreg
    x1$recaps$matched[[paste(column_adjust, "adjusted", event_adjust, sep="_")]] <- yreg    # reordered
  } else {
    x1$recaps$matched[[paste(column_adjust, "adjusted", event_adjust, sep="_")]] <- unname(predict(lm1))
  }
  
  # need to change it in x1$recaps$unmatched[[event_adjust]][[column_adjust]]
  x1$recaps$unmatched[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]] <-
    unname(predict(lm1, newdata = data.frame(xreg=x1$recaps$unmatched[[event_adjust]][[column_adjust]])))
  
  # need to change it in x1$recaps$all[[event_adjust]][[column_adjust]]
  # x1$recaps$all[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]] <-
  #   unname(predict(lm1, newdata = data.frame(xreg=x1$recaps$all[[event_adjust]][[column_adjust]])))
  # if(impute) {
  #   ytag <- x1$recaps$all[[event_adjust]][[ID_adjust]]
  #   keeptag <- x$recaps$matched[[paste(ID_keep, event_keep, sep="_")]]
  #   for(iy in seq_along(ytag)) {
  #     if(!is.na(ytag[iy])) {
  #       if(ytag[iy] %in% keeptag) {
  #         x1$recaps$all[[event_adjust]][[paste(column_adjust, "adjusted", sep="_")]][iy] <- 
  #           yreg[keeptag==ytag[iy]]  
  #       }
  #     }
  #   }
  # }
  x1$recaps$all[[paste(column_adjust, "adjusted", event_adjust, sep="_")]] <-
    unname(predict(lm1, newdata = data.frame(xreg=x1$recaps$all[[paste(column_adjust, event_adjust, sep="_")]])))
  if(impute) {
    ytag <- x1$recaps$all[[paste(ID_adjust, event_adjust, sep="_")]]
    keeptag <- x$recaps$matched[[paste(ID_keep, event_keep, sep="_")]]
    for(iy in seq_along(ytag)) {
      if(!is.na(ytag[iy])) {
        if(ytag[iy] %in% keeptag) {
          x1$recaps$all[[paste(column_adjust, "adjusted", event_adjust, sep="_")]][iy] <- 
            yreg[keeptag==ytag[iy]]  
        }
      }
    }
  }
  
  return(x1)
}



# pulling out the common pieces of error checking
recapr_errorcheck <- function(x, event_names, column_names) {
  # error checking
  # - needs to be an object returned by recapr_prep
  if(!inherits(x, "MR_data")) stop("Argument x= must be an object returned from recapr_prep()")
  
  if(length(event_names) > 2) stop("No more than two event names can be used")
  if(length(column_names) > 2) stop("No more than two column names can be used")
  
  ## only require event names when multiple column names are used
  ## allow length-1 column name
  if(length(column_names) == 2) {
    if(length(event_names) == 0) {
      if(column_names[1] != column_names[2]) {
        stop("Must supply two event names if two column names are used")
      } else {
        event_names <- names(x$input_data)
      }
    }
    if(length(event_names) == 1) stop("Must supply two event names if two column names are used")
  } 
  if(length(column_names) == 1) {
    column_names <- rep(column_names, 2)
    if(is.null(event_names)) {
      event_names <- names(x$input_data)
    }
    if(length(event_names) == 1) stop("Invalid input to event_names")
  }
  
  # # must be two event names and two column names
  # if(length(event_names) != 2) stop("Must supply two event names")
  # if(length(column_names) != 2) stop("Must supply two column names")
  
  # - event names must exist
  if(!(all(event_names %in% names(x$input_data)))) {
    stop("Supplied event names not found in data event names")
  }
  
  # - column names must exist in both events
  if(!(column_names[1] %in% names(x$input_data[[event_names[1]]]))) {
    stop(paste("Column name not found in", event_names[1], "data"))
  }
  if(!(column_names[2] %in% names(x$input_data[[event_names[2]]]))) {
    stop(paste("Column name not found in", event_names[2], "data"))
  }
  
  return(list(event_names=event_names, column_names=column_names))
}



# big function to truncate objects according to a min and max value
truncate <- function(x, event_names, column_names, min=NULL, max=NULL) {
  
  # error checking
  # - needs to be an object returned by recapr_prep
  if(!inherits(x, "MR_data")) stop("Argument x= must be an object returned from recapr_prep()")

  if(length(event_names) > 2) stop("No more than two event names can be used")
  if(length(column_names) > 2) stop("No more than two column names can be used")

  ## only require event names when multiple column names are used
  ## allow length-1 column name
  if(length(column_names) == 2) {
    if(length(event_names) == 0) {
      if(column_names[1] != column_names[2]) {
        stop("Must supply two event names if two column names are used")
      } else {
        event_names <- names(x$input_data)
      }
    }
    if(length(event_names) == 1) stop("Must supply two event names if two column names are used")
  }
  if(length(column_names) == 1) {
    column_names <- rep(column_names, 2)
    if(is.null(event_names)) {
      event_names <- names(x$input_data)
    }
    if(length(event_names) == 1) stop("Invalid input to event_names")
  }

  # # must be two event names and two column names
  # if(length(event_names) != 2) stop("Must supply two event names")
  # if(length(column_names) != 2) stop("Must supply two column names")

  # - event names must exist
  if(!(all(event_names %in% names(x$input_data)))) {
    stop("Supplied event names not found in data event names")
  }

  # - column names must exist in both events
  if(!(column_names[1] %in% names(x$input_data[[event_names[1]]]))) {
    stop(paste("Column name not found in", event_names[1], "data"))
  }
  if(!(column_names[2] %in% names(x$input_data[[event_names[2]]]))) {
    stop(paste("Column name not found in", event_names[2], "data"))
  }
  # ec <- recapr_errorcheck(x=x, event_names=event_names, column_names=column_names)
  # event_names <- ec$event_names
  # column_names <- ec$column_names
  
  # - columns must be numeric in both events
  if(!is.numeric(x$input_data[[event_names[1]]][[column_names[1]]]) | 
     !is.numeric(x$input_data[[event_names[2]]][[column_names[2]]])) {
    stop("Non-numeric columns detected")
  }
  
  # - min and max must be numeric if they are not null
  if(!is.null(min) & !is.numeric(min)) stop("Argument min= must be numeric")
  if(!is.null(max) & !is.numeric(max)) stop("Argument max= must be numeric")
  
  # create a copy to modify
  x1 <- x
  
  # x1$input_data[[event_names[1]]] %>% str
  # x1$input_data[[event_names[2]]] %>% str
  # x1$recaps$unmatched[[event_names[1]]] %>% str
  # x1$recaps$unmatched[[event_names[2]]] %>% str
  # x1$recaps$all[[event_names[1]]] %>% str
  # x1$recaps$all[[event_names[2]]] %>% str
  for(ii in 1:2) {
    if(!is.null(min)) {
      x1$input_data[[event_names[ii]]] <- subset(x1$input_data[[event_names[ii]]], 
                                                 x1$input_data[[event_names[ii]]][[column_names[ii]]] >= min)
      x1$recaps$unmatched[[event_names[ii]]] <- subset(x1$recaps$unmatched[[event_names[ii]]],
                                                       x1$recaps$unmatched[[event_names[ii]]][[column_names[ii]]] >= min)
      # x1$recaps$all[[event_names[ii]]] <- subset(x1$recaps$all[[event_names[ii]]],
      #                                            x1$recaps$all[[event_names[ii]]][[column_names[ii]]] >= min)
    }
    if(!is.null(max)) {
      x1$input_data[[event_names[ii]]] <- subset(x1$input_data[[event_names[ii]]], 
                                                 x1$input_data[[event_names[ii]]][[column_names[ii]]] <= max)
      x1$recaps$unmatched[[event_names[ii]]] <- subset(x1$recaps$unmatched[[event_names[ii]]],
                                                       x1$recaps$unmatched[[event_names[ii]]][[column_names[ii]]] <= max)
      # x1$recaps$all[[event_names[ii]]] <- subset(x1$recaps$all[[event_names[ii]]],
      #                                            x1$recaps$all[[event_names[ii]]][[column_names[ii]]] <= max)
    }
  }
  
  # check to make sure subsetting agrees for matched fish
  logi1 <- logi2 <- rep(TRUE, nrow(x1$recaps$matched))
  if(!is.null(min)) {
    logi1[x1$recaps$matched[[paste(column_names[1], event_names[1], sep="_")]] < min] <- FALSE
    logi2[x1$recaps$matched[[paste(column_names[2], event_names[2], sep="_")]] < min] <- FALSE
  }
  if(!is.null(max)) {
    logi1[x1$recaps$matched[[paste(column_names[1], event_names[1], sep="_")]] > max] <- FALSE
    logi2[x1$recaps$matched[[paste(column_names[2], event_names[2], sep="_")]] > max] <- FALSE
  }
  if(!all(logi1 == logi2)) {
    warning(paste("Disagreement in truncation for", sum(logi1 != logi2), "recaptured individuals"))
  }
  
  # x1$recaps$matched %>% str
  if(!is.null(min)) {
    x1$recaps$matched <- subset(x1$recaps$matched,
                                x1$recaps$matched[[paste(column_names[1], event_names[1], sep="_")]] >= min & 
                                  x1$recaps$matched[[paste(column_names[2], event_names[2], sep="_")]] >= min)
    x1$recaps$all <- subset(x1$recaps$all,
                                (x1$recaps$all[[paste(column_names[1], event_names[1], sep="_")]] >= min & 
                                  x1$recaps$all[[paste(column_names[2], event_names[2], sep="_")]] >= min) |
                              (x1$recaps$all[[paste(column_names[1], event_names[1], sep="_")]] >= min & 
                                 is.na(x1$recaps$all[[paste(column_names[2], event_names[2], sep="_")]])) |
                              (is.na(x1$recaps$all[[paste(column_names[1], event_names[1], sep="_")]]) & 
                                 x1$recaps$all[[paste(column_names[2], event_names[2], sep="_")]] >= min))
      
  }
  if(!is.null(max)) {
    x1$recaps$matched <- subset(x1$recaps$matched,
                                x1$recaps$matched[[paste(column_names[1], event_names[1], sep="_")]] <= max & 
                                  x1$recaps$matched[[paste(column_names[2], event_names[2], sep="_")]] <= max)
    x1$recaps$all <- subset(x1$recaps$all,
                            (x1$recaps$all[[paste(column_names[1], event_names[1], sep="_")]] <= max & 
                               x1$recaps$all[[paste(column_names[2], event_names[2], sep="_")]] <= max) |
                              (x1$recaps$all[[paste(column_names[1], event_names[1], sep="_")]] <= max & 
                                 is.na(x1$recaps$all[[paste(column_names[2], event_names[2], sep="_")]])) |
                              (is.na(x1$recaps$all[[paste(column_names[1], event_names[1], sep="_")]]) & 
                                 x1$recaps$all[[paste(column_names[2], event_names[2], sep="_")]] <= max))
  }
  
  return(x1)
}



stratify <- function(x, event_names=NULL, column_names, breaks, right=FALSE, dig.lab=6) {
  
  # error checking
  # - needs to be an object returned by recapr_prep
  if(!inherits(x, "MR_data")) stop("Argument x= must be an object returned from recapr_prep()")

  if(length(event_names) > 2) stop("No more than two event names can be used")
  if(length(column_names) > 2) stop("No more than two column names can be used")

  ## only require event names when multiple column names are used
  ## allow length-1 column name
  if(length(column_names) == 2) {
    if(length(event_names) == 0) {
      if(column_names[1] != column_names[2]) {
        stop("Must supply two event names if two column names are used")
      } else {
        event_names <- names(x$input_data)
      }
    }
    if(length(event_names) == 1) stop("Must supply two event names if two column names are used")
  }
  if(length(column_names) == 1) {
    column_names <- rep(column_names, 2)
    if(is.null(event_names)) {
      event_names <- names(x$input_data)
    }
    if(length(event_names) == 1) stop("Invalid input to event_names")
  }

  # # must be two event names and two column names
  # if(length(event_names) != 2) stop("Must supply two event names")
  # if(length(column_names) != 2) stop("Must supply two column names")

  # - event names must exist
  if(!(all(event_names %in% names(x$input_data)))) {
    stop("Supplied event names not found in data event names")
  }

  # - column names must exist in both events
  if(!(column_names[1] %in% names(x$input_data[[event_names[1]]]))) {
    stop(paste("Column name not found in", event_names[1], "data"))
  }
  if(!(column_names[2] %in% names(x$input_data[[event_names[2]]]))) {
    stop(paste("Column name not found in", event_names[2], "data"))
  }
  # ec <- recapr_errorcheck(x=x, event_names=event_names, column_names=column_names)
  # event_names <- ec$event_names
  # column_names <- ec$column_names
  
  # - columns must be numeric in both events
  if(!is.numeric(x$input_data[[event_names[1]]][[column_names[1]]]) | 
     !is.numeric(x$input_data[[event_names[2]]][[column_names[2]]])) {
    stop("Non-numeric columns detected")
  }
  
  # - breaks must be numeric if they are not null
  if(!is.null(breaks) & !is.numeric(breaks)) stop("Argument breaks= must be numeric")
  
  # create a copy to modify
  x1 <- x
  
  # create a dummy factor to make sure levels are consistent and check range
  levelmaker <- c(x1$input_data[[event_names[1]]][[column_names[1]]],
                  x1$input_data[[event_names[2]]][[column_names[2]]])
  testcut <- cut(levelmaker, breaks=breaks, right=right, dig.lab=dig.lab)
  thelevels <- levels(testcut)
  
  # checking range
  if(right) {
    if(min(levelmaker, na.rm=TRUE) <= min(breaks) | max(levelmaker, na.rm=TRUE) > max(breaks)) {
      warning("Data exists beyond range of stratum breaks: min and max must be included")
    }
  } else {
    if(min(levelmaker, na.rm=TRUE) < min(breaks) | max(levelmaker, na.rm=TRUE) >= max(breaks)) {
      warning("Data exists beyond range of stratum breaks: min and max must be included")
    }
  }
  
  # creating factor columns
  for(ii in 1:2) {
    x1$input_data[[event_names[ii]]][[paste(column_names[ii], "strat", sep="_")]] <- 
      factor(cut(x1$input_data[[event_names[ii]]][[column_names[ii]]], 
                 breaks=breaks, right=right, dig.lab=dig.lab), levels=thelevels)
    
    x1$recaps$unmatched[[event_names[ii]]][[paste(column_names[ii], "strat", sep="_")]] <- 
      factor(cut(x1$recaps$unmatched[[event_names[ii]]][[column_names[ii]]], 
                 breaks=breaks, right=right, dig.lab=dig.lab), levels=thelevels)
    
    # x1$recaps$all[[event_names[ii]]][[paste(column_names[ii], "strat", sep="_")]] <- 
    #   factor(cut(x1$recaps$all[[event_names[ii]]][[column_names[ii]]], 
    #              breaks=breaks, right=right, dig.lab=dig.lab), levels=thelevels)
    
    x1$recaps$all[[paste(column_names[ii], "strat", event_names[ii], sep="_")]] <- 
      factor(cut(x1$recaps$all[[paste(column_names[ii], event_names[ii], sep="_")]], 
                 breaks=breaks, right=right, dig.lab=dig.lab), levels=thelevels)
    
    x1$recaps$matched[[paste(column_names[ii], "strat", event_names[ii], sep="_")]] <- 
      factor(cut(x1$recaps$matched[[paste(column_names[ii], event_names[ii], sep="_")]], 
                 breaks=breaks, right=right, dig.lab=dig.lab), levels=thelevels)
  }
  
  # checking if there is disagreement in $recaps$matched factor columns
  s1 <- x1$recaps$matched[[paste(column_names[1], "strat", event_names[1], sep="_")]]
  s2 <- x1$recaps$matched[[paste(column_names[2], "strat", event_names[2], sep="_")]]
  if(!isTRUE(all.equal(s1, s2))) {
    warning(paste("Disagreement in stratum assignment for", 
                  sum(s1 != s2, na.rm=TRUE) + sum(is.na(s1) != is.na(s2)), 
                  "recaptured individuals"))
  }
  
  return(x1)
}  



tabulate_samples <- function(x,
                             column_names=NULL,
                             event_names=NULL) {
  # error checking
  # - needs to be an object returned by recapr_prep
  if(!inherits(x, "MR_data")) stop("Argument x= must be an object returned from recapr_prep()")
  
  if(length(event_names) > 2) stop("No more than two event names can be used")
  if(length(column_names) > 2) stop("No more than two column names can be used")
  
  ## only require event names when multiple column names are used
  ## allow length-1 column name
  if(length(column_names) == 2) {
    if(length(event_names) == 0) {
      if(column_names[1] != column_names[2]) {
        stop("Must supply two event names if two column names are used")
      } else {
        event_names <- names(x$input_data)
      }
    }
    if(length(event_names) == 1) stop("Must supply two event names if two column names are used")
  }
  if(length(column_names) == 1) {
    column_names <- rep(column_names, 2)
    if(is.null(event_names)) {
      event_names <- names(x$input_data)
    }
    if(length(event_names) == 1) stop("Invalid input to event_names")
  }
  if(is.null(column_names) & is.null(event_names)) {
    event_names <- names(x$input_data)
  }
  
  # # must be two event names and two column names
  # if(length(event_names) != 2) stop("Must supply two event names")
  # if(length(column_names) != 2) stop("Must supply two column names")
  
  # - event names must exist
  if(!(all(event_names %in% names(x$input_data)))) {
    stop("Supplied event names not found in data event names")
  }
  
  if(!is.null(column_names)) {  # this condition does not exist in other funcs
    # - column names must exist in both events
    if(!(column_names[1] %in% names(x$input_data[[event_names[1]]]))) {
      stop(paste("Column name not found in", event_names[1], "data"))
    }
    if(!(column_names[2] %in% names(x$input_data[[event_names[2]]]))) {
      stop(paste("Column name not found in", event_names[2], "data"))
    }
  }
  # ec <- recapr_errorcheck(x=x, event_names=event_names, column_names=column_names)
  # event_names <- ec$event_names
  # column_names <- ec$column_names
  
  ### START TABULATING!
  out <- list(captures=list(), recaps=list())
  if(!is.null(column_names)) {  # if there are strata
    out$captures[[event_names[1]]] <- table(x$input_data[[event_names[1]]][[column_names[1]]], useNA = "ifany")
    out$captures[[event_names[2]]] <- table(x$input_data[[event_names[2]]][[column_names[2]]], useNA = "ifany")
    out$recaps$matched <- table(factor(x$recaps$matched[[paste(column_names[1], event_names[1], sep="_")]],
                                       levels = names(out$captures[[event_names[1]]])),
                                factor(x$recaps$matched[[paste(column_names[2], event_names[2], sep="_")]],
                                       levels = names(out$captures[[event_names[2]]])),
                                useNA = "ifany", dnn=paste(column_names, event_names, sep="_"))
    out$recaps$unmatched <- list()
    out$recaps$unmatched[[event_names[1]]] <- table(factor(x$recaps$unmatched[[event_names[1]]][[column_names[1]]],
                                                           levels = names(out$captures[[event_names[1]]])), useNA="ifany")
    out$recaps$unmatched[[event_names[2]]] <- table(factor(x$recaps$unmatched[[event_names[2]]][[column_names[2]]],
                                                           levels = names(out$captures[[event_names[2]]])), useNA="ifany")
    out$recaps$all <- table(factor(x$recaps$all[[paste(column_names[1], event_names[1], sep="_")]],
                                   levels = names(out$captures[[event_names[1]]])),
                            factor(x$recaps$all[[paste(column_names[2], event_names[2], sep="_")]],
                                   levels = names(out$captures[[event_names[2]]])),
                            useNA = "ifany", dnn=paste(column_names, event_names, sep="_"))
  } else { # if there are no strata
    out$captures[[event_names[1]]] <- nrow(x$input_data[[event_names[1]]])
    out$captures[[event_names[2]]] <- nrow(x$input_data[[event_names[2]]])
    out$recaps$matched <- nrow(x$recaps$matched)
    out$recaps$unmatched <- list()
    out$recaps$unmatched[[event_names[1]]] <- nrow(x$recaps$unmatched[[event_names[1]]])
    out$recaps$unmatched[[event_names[2]]] <- nrow(x$recaps$unmatched[[event_names[2]]])
    out$recaps$all <- nrow(x$recaps$all)
  }
  
  return(out)
}



############ ---------- testing zone ------------ #############

### load packages
library(tidyverse)
library(recapr)
library(dsftools)  # devtools::install_github("ADFG-DSF/dsftools")


### read data
Event1 <- read_csv("FDS_2024/flat_data/Event1.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`))

Event2 <- read_csv("FDS_2024/flat_data/Event2.csv", skip = 1) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::remove_empty(which = "rows") %>%
  mutate(`Tag Number` = as.character(`Tag Number`))



aa <- recapr_prep(ID="Tag Number", event1=Event1, event2=Event2, recap_codes="TL")
# aa <- recapr_prep(ID="Tag Number", event=NULL, recap_codes="TL", Event1, Event2)


aa1 <- correct_growth(x=aa, 
                      # impute=FALSE, 
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

# plot(aa1$recaps$matched$`Fork Length (mm)_event2`, aa1$recaps$matched$`Fork Length (mm)_event2_adjusted`)
# plot(aa1$recaps$matched$`Fork Length (mm)_event1`, aa1$recaps$matched$`Fork Length (mm)_event2_adjusted`)
plot(aa1$recaps$matched$`Fork Length (mm)_event2`, aa1$recaps$matched$`Fork Length (mm)_adjusted_event2`)
plot(aa1$recaps$matched$`Fork Length (mm)_event1`, aa1$recaps$matched$`Fork Length (mm)_adjusted_event2`)
abline(0, 1, lty=3)

plot(aa1$recaps$unmatched$event2$`Fork Length (mm)`, aa1$recaps$unmatched$event2$`Fork Length (mm)_adjusted`)
abline(lm1)

plot(aa1$recaps$all$event2$`Fork Length (mm)`, aa1$recaps$all$event2$`Fork Length (mm)_adjusted`,
     pch=ifelse(aa1$recaps$all$event2$`Tag Number` %in% aa$recaps$matched$`Tag Number_event1`, 16, 1))
plot(aa1$recaps$all$`Fork Length (mm)_event2`, aa1$recaps$all$`Fork Length (mm)_adjusted_event2`,
     pch=ifelse(aa1$recaps$all$`Tag Number_event2` %in% aa$recaps$matched$`Tag Number_event1`, 16, 1))
abline(lm1)

all.equal(aa1$input_data$event1, aa$input_data$event1)
all.equal(aa1$recaps$unmatched$event1, aa$recaps$unmatched$event1)
all.equal(aa1$recaps$all$event1, aa$recaps$all$event1)



aa2 <- truncate(x = aa1,
                event_names = c("event1", "event2"),
                column_names = c("Fork Length (mm)", "Fork Length (mm)_adjusted"),
                min = 345, max = 400)
str(aa2)
range(aa1$input_data$event1$`Fork Length (mm)`)
range(aa2$input_data$event1$`Fork Length (mm)`)
range(aa1$input_data$event2$`Fork Length (mm)_adjusted`)
range(aa2$input_data$event2$`Fork Length (mm)_adjusted`)
range(aa1$recaps$matched$`Fork Length (mm)_event1`)
range(aa2$recaps$matched$`Fork Length (mm)_event1`)
range(aa1$recaps$matched$`Fork Length (mm)_adjusted_event2`)
range(aa2$recaps$matched$`Fork Length (mm)_adjusted_event2`)
range(aa1$recaps$unmatched$event1$`Fork Length (mm)`)
range(aa2$recaps$unmatched$event1$`Fork Length (mm)`)
range(aa1$recaps$unmatched$event2$`Fork Length (mm)`)
range(aa2$recaps$unmatched$event2$`Fork Length (mm)_adjusted`)
# range(aa1$recaps$all$event1$`Fork Length (mm)`)
# range(aa2$recaps$all$event1$`Fork Length (mm)`)
# range(aa1$recaps$all$event2$`Fork Length (mm)_adjusted`)
# range(aa2$recaps$all$event2$`Fork Length (mm)_adjusted`)
range(aa1$recaps$all$`Fork Length (mm)_event1`, na.rm=TRUE)
range(aa2$recaps$all$`Fork Length (mm)_event1`, na.rm=TRUE)
range(aa1$recaps$all$`Fork Length (mm)_adjusted_event2`, na.rm=TRUE)
range(aa2$recaps$all$`Fork Length (mm)_adjusted_event2`, na.rm=TRUE)


aa3 <- stratify(x = aa1,
                event_names = c("event1", "event2"),
                column_names = c("Fork Length (mm)", "Fork Length (mm)_adjusted"),#
                breaks = c(200, 345, 400, 1000))
str(aa3)
plot(aa3$input_data$event1$`Fork Length (mm)` ~ aa3$input_data$event1$`Fork Length (mm)_strat`)
plot(aa3$input_data$event2$`Fork Length (mm)_adjusted` ~ aa3$input_data$event2$`Fork Length (mm)_adjusted_strat`)
plot(aa3$recaps$matched$`Fork Length (mm)_event1` ~ aa3$recaps$matched$`Fork Length (mm)_strat_event1`)
plot(aa3$recaps$matched$`Fork Length (mm)_adjusted_event2` ~ aa3$recaps$matched$`Fork Length (mm)_adjusted_strat_event2`)
plot(aa3$recaps$unmatched$event2$`Fork Length (mm)_adjusted` ~ aa3$recaps$unmatched$event2$`Fork Length (mm)_adjusted_strat`)
# plot(aa3$recaps$all$event1$`Fork Length (mm)` ~ aa3$recaps$all$event1$`Fork Length (mm)_strat`)
# plot(aa3$recaps$all$event2$`Fork Length (mm)_adjusted` ~ aa3$recaps$all$event2$`Fork Length (mm)_adjusted_strat`)
plot(aa3$recaps$all$`Fork Length (mm)_event1` ~ aa3$recaps$all$`Fork Length (mm)_strat_event1`)
plot(aa3$recaps$all$`Fork Length (mm)_adjusted_event2` ~ aa3$recaps$all$`Fork Length (mm)_adjusted_strat_event2`)



Event2na <- Event2
Event2na$`Fork Length (mm)`[Event2na$`Tag Number` == 770] <- NA

aana <- recapr_prep(ID="Tag Number", event1=Event1, event2=Event2na, recap_codes="TL")
all.equal(aa, aana)

aa1na <- correct_growth(x=aana, 
                        event_keep="event1", 
                        event_adjust="event2", 
                        column_keep="Fork Length (mm)", 
                        column_adjust="Fork Length (mm)",
                        ID_keep="Tag Number",
                        ID_adjust="Tag Number")
all.equal(aa1, aa1na)


## should do these without growth correcting ##

aa2na <- truncate(x = aana,
                  event_names = c("event1", "event2"),
                  column_names = c("Fork Length (mm)", "Fork Length (mm)"),
                  min = 345, max = 400)
aa2b <- truncate(x = aa,
                  event_names = c("event1", "event2"),
                  column_names = c("Fork Length (mm)", "Fork Length (mm)"),
                  min = 345, max = 400)
all.equal(aa2b, aa2na)
dim(aa2na$input_data$event1)
dim(aa2b$input_data$event1)
dim(aa2na$input_data$event2)
dim(aa2b$input_data$event2)
dim(aa2na$recaps$matched)
dim(aa2b$recaps$matched)
dim(aa2na$recaps$unmatched$event1)
dim(aa2b$recaps$unmatched$event1)
# dim(aa2na$recaps$all$event1)
# dim(aa2b$recaps$all$event1)
dim(aa2na$recaps$all)
dim(aa2b$recaps$all)
# dim(aa2na$recaps$all$event2)
# dim(aa2b$recaps$all$event2)



aa3na <- stratify(x = aana,
                  event_names = c("event1", "event2"),
                  column_names = c("Fork Length (mm)", "Fork Length (mm)"),#
                  breaks = c(200, 345, 400, 1000))
aa3b <- stratify(x = aa,
                  event_names = c("event1", "event2"),
                  column_names = c("Fork Length (mm)", "Fork Length (mm)"),#
                  breaks = c(200, 345, 400, 1000))
all.equal(aa3b, aa3na)
dim(aa3na$input_data$event1)
dim(aa3b$input_data$event1)
dim(aa3na$input_data$event2)
dim(aa3b$input_data$event2)
dim(aa3na$recaps$matched)
dim(aa3b$recaps$matched)
dim(aa3na$recaps$unmatched$event1)
dim(aa3b$recaps$unmatched$event1)
# dim(aa3na$recaps$all$event1)
# dim(aa3b$recaps$all$event1)
# dim(aa3na$recaps$all$event2)
# dim(aa3b$recaps$all$event2)
dim(aa3na$recaps$all)
dim(aa3b$recaps$all)

table(aa3na$input_data$event1$`Fork Length (mm)_strat`, useNA = "always")
table(aa3b$input_data$event1$`Fork Length (mm)_strat`, useNA = "always")
table(aa3na$input_data$event2$`Fork Length (mm)_strat`, useNA = "always")
table(aa3b$input_data$event2$`Fork Length (mm)_strat`, useNA = "always")
table(aa3na$recaps$matched$`Fork Length (mm)_strat_event1`, useNA = "always")
table(aa3b$recaps$matched$`Fork Length (mm)_strat_event1`, useNA = "always")
table(aa3na$recaps$matched$`Fork Length (mm)_strat_event2`, useNA = "always")
table(aa3b$recaps$matched$`Fork Length (mm)_strat_event2`, useNA = "always")
table(aa3na$recaps$unmatched$event1$`Fork Length (mm)_strat`, useNA = "always")
table(aa3b$recaps$unmatched$event1$`Fork Length (mm)_strat`, useNA = "always")
table(aa3na$recaps$unmatched$event2$`Fork Length (mm)_strat`, useNA = "always")
table(aa3b$recaps$unmatched$event2$`Fork Length (mm)_strat`, useNA = "always")
# table(aa3na$recaps$all$event1$`Fork Length (mm)_strat`, useNA = "always")
# table(aa3b$recaps$all$event1$`Fork Length (mm)_strat`, useNA = "always")
# table(aa3na$recaps$all$event2$`Fork Length (mm)_strat`, useNA = "always")
# table(aa3b$recaps$all$event2$`Fork Length (mm)_strat`, useNA = "always")
table(aa3na$recaps$all$`Fork Length (mm)_strat_event1`, useNA="always")
table(aa3b$recaps$all$`Fork Length (mm)_strat_event1`, useNA="always")
table(aa3na$recaps$all$`Fork Length (mm)_strat_event2`, useNA="always")
table(aa3b$recaps$all$`Fork Length (mm)_strat_event2`, useNA="always")


Event1TL <- Event1
Event1TL$`Tag Number`[1] <- "TL"
aaTL <- recapr_prep(ID="Tag Number", event1=Event1TL, event2=Event2, recap_codes="TL")
str(aaTL$recaps)
aaTL$recaps$matched$`Tag Number_event1`   # seems to not include TL??
aaTL$recaps$unmatched$event1$`Tag Number`
aaTL$recaps$all$`Tag Number_event1`


# could theoretically
# - automate ks tests
# - automate chisq tests -> make inputs for recapr::consistencytest
# DONE - df of all unique fish ***************************************
#   * or make master length column somehow
#   DONE * REWORK $recaps$all to semi-interleaved **************** !!!!!!!!!!
#   DONE   - and of course change truncate/stratify/correct_growth
# - tabulate stratificationses

# edge cases i can think of
# DONE - will $matched ever be constructed of things with different lengths? NO
# DONE - what happens when there are NA values in (length)  ALL FINE
#   * truncate  REMOVES NA
#   * stratify  MAKES NA
#   * correct_growth MAKES NA
# - (maybe change adjusted to corrected, or else correct_growth to adjust_growth or growth_correction)
# DONE - could there be NA in ID column? NO
# DONE - what happens when data is not named in recapr_prep? MAKES NEW NAMES
# DONE - what happens when min and max are left null?  NOTHING
# - what happens when elements of recap_codes are repeated in both events?

# idea: add length and stratification column names for each event to MR_data object
# print or summary method: tabulate counts by stratum?? range of numeric values? pop estimates??
# DONE in correct_growth: add a logical argument for whether to directly impute matched 

# problem: how to make lengths/stratum/counts consistent for $recaps

## DONE robustify default column names etc (allow length-1 if column is present)




# # make a toy dataset!
# n1 <- 30
# n2 <- 30
# m2 <- 10
# cap1 <- data.frame(tag = 1:n1,
#                    FL = round(rnorm(n1, mean=400, sd=50)),
#                    area = sample(LETTERS[1:3], n1, replace=TRUE))
# cap2 <- data.frame(tag = NA,
#                    FL = round(rnorm(n1, mean=400, sd=50)),
#                    area = sample(LETTERS[1:3], n1, replace=TRUE))
# cap1_recaps <- sample(n1, m2)
# cap2_recaps <- sample(n2, m2)
# cap2$tag[cap2_recaps] <- cap1$tag[cap1_recaps]
# cap2$FL[cap2_recaps] <- cap1$FL[cap1_recaps] + round(rnorm(m2, mean=10, sd=4))
# 
# cap1$tag[sample(n1, 1)] <- "TL"
# cap1$FL[sample(n1, 1)] <- NA
# cap1$area[sample(n1, 1)] <- NA
# 
# cap2$tag[sample(cap2_recaps, 2)] <- "TL"
# cap2$FL[sample(n2, 1)] <- NA



cap1 <-
  structure(list(tag = c(1L, "TL", 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 
                         11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 
                         24L, 25L, 26L, 27L, 28L, 29L, 30L), 
                 FL = c(320, 463, 397, 447, 
                        422, 476, 342, 416, 395, NA, 461, 377, 414, 387, 310, 313, 396, 
                        432, 374, 487, 361, 426, 369, 360, 406, 367, 413, 329, 400, 340), 
                 area = c("B", "B", "A", "B", "A", "A", "B", "A", NA, "A", 
                          "C", "C", "C", "B", "C", "C", "A", "B", "A", "C", "B", "A", "A", 
                          "A", "A", "C", "C", "C", "C", "B")), 
            row.names = c(NA, -30L), class = "data.frame")
cap2 <-
  structure(list(tag = c("25", NA, NA, NA, NA, "17", "14", NA, 
                         NA, "2", NA, NA, NA, NA, "TL", NA, NA, "29", NA, NA, "20", NA, 
                         "9", NA, NA, "13", NA, NA, "TL", NA), 
                 FL = c(411, 405, 425, 436, 
                        395, 408, 399, 441, 309, 478, 341, 443, 368, 305, 361, NA, 317, 
                        411, 398, 387, 501, 390, 402, 372, 367, 424, 346, 360, 420, 411), 
                 area = c("B", "B", "C", "A", "C", "A", "B", "A", "A", "A", 
                          "B", "A", "A", "B", "A", "B", "B", "A", "C", "B", "A", "A", "A", 
                          "B", "C", "B", "A", "B", "C", "B")), 
            row.names = c(NA, -30L), class = "data.frame")


x <- recapr_prep(ID="tag", recap_codes = "TL", cap1=cap1, cap2=cap2)
x
truncate(x=x, event_names=c("cap1","cap2"), column_names = c("FL","FL"), min=400, max=450)
stratify(x=x, event_names=c("cap1","cap2"), column_names = c("FL","FL"), breaks=c(300, 400, 500))
correct_growth(x=x, event_keep="cap1", event_adjust="cap2", column_keep="FL", column_adjust="FL", ID_keep="tag")#, impute=F)

stratify(x=x, column_names = c("FL", "FL"), breaks=c(300, 400, 550), event_names=c("cap1","cap2"))
stratify(x=x, column_names = c("FL"), breaks=c(300, 400, 550), event_names=c("cap1","cap2"))
stratify(x=x, column_names = c("FL", "FL"), breaks=c(300, 400, 550))   
stratify(x=x, column_names = c("FL"), breaks=c(300, 400, 550))
stratify(x=x, column_names = c("FL", "FL"), breaks=c(300, 400, 550), event_names=c("cap1"))
stratify(x=x, column_names = c("FL"), breaks=c(300, 400, 550), event_names=c("cap1"))
stratify(x=x, column_names = c("FL", "FL", "FL"), breaks=c(300, 400, 550))
stratify(x=x, breaks=c(300, 400, 550))
stratify(x=x, column_names = c("FL", "FL1"), breaks=c(300, 400, 550)) 
stratify(x=x, column_names = c("FL", "FL1"), breaks=c(300, 400, 550), event_names=c("cap1","cap2")) 
tabulate_samples(x=x, column_names = "area")
tabulate_samples(x)
