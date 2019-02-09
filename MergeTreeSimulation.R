# We are trying to merging multiple trees into one single cohort definition
library(MASS)

rm(list = ls())      # Clear all variables
graphics.off()    # Close graphics windows

########################### Helper Functions ##########################
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

calculatePvalue <- function(df, metric.idx, treatment.idx) {
  treatment.ids <- which(df[[treatment.idx]]  == 1)
  control.ids <- which(df[[treatment.idx]] == 0)
  y.1 <- mean(df[, metric.idx][treatment.ids])
  y.0 <- mean(df[, metric.idx][control.ids])
  n.1 <- as.numeric(nrow(df[treatment.ids, ]))
  n.0 <- as.numeric(nrow(df[control.ids, ]))
  var.y1 <- var(df[, metric.idx][treatment.ids]) / n.1
  var.y0 <- var(df[, metric.idx][control.ids]) / n.0
  delta <- y.1 - y.0

  variance.delta <- var.y1 + var.y0
  tstats.delta <- delta / sqrt(variance.delta)
  pValue.delta <- 2 * pnorm(- abs(tstats.delta))
  pValue.delta
}

extractSet <- function(str){
  if (str == "NULL") set <- NULL
  else set <- as.numeric(trim(unlist(strsplit(str, ","))))
  set
}

#take intersection of two Segments by evaluate the sets per feature
intersectSegment <- function(S.1, S.2){
  Segment.out <- list()
  for (i in 1:length(S.1)){
    set.1 <- extractSet(S.1[[i]])
    set.2 <- extractSet(S.2[[i]])

    set.intersect <- intersect(set.1, set.2)
    if (is.null(set.intersect) || length(set.intersect) == 0) Segment.out[[i]] <- "NULL"
    else Segment.out[[i]] <- paste(noquote(set.intersect), collapse=",")
  }
  Segment.out
}

#Get (S.1 - S.2) by evaluate the sets per feature
complementSegment <- function(S.1, S.2){
  Segment.out <- list()
  for (i in 1:length(S.1)){
    set.1 <- extractSet(S.1[[i]])
    set.2 <- extractSet(S.2[[i]])

    set.complement <- setdiff(set.1, set.2)
    if (is.null(set.complement) || length(set.complement) == 0) Segment.out[[i]] <- "NULL"
    else Segment.out[[i]] <- paste(noquote(set.complement), collapse=",")
  }
  Segment.out
}

buildSubtableSet <- function(Segments, Metrics){
  subtable.set <- list()
  for (j in 1:length(Segments)){
    segment <- Segments[[j]]
    subtable <- Metrics
    for (i in 1: length(segment)){
      set.1 <- extractSet(segment[[i]])
      subtable <- subtable[which(subtable[,i] %in% set.1), ]
    }
    subtable.set[[j]] <- subtable
  }
  subtable.set
}

isStatSig <- function(Segments, Metrics, metrics.idx.list){
  subtable.set <- buildSubtableSet(Segments, Metrics)
  for (k in metrics.idx.list){
    flag <- TRUE
    for (subtable in subtable.set){
      if (nrow(subtable) == 0 || calculatePvalue(subtable, k, 6) > 0.1){
        flag <- FALSE
      }
    }
    if (flag == TRUE) return(flag)
  }
  FALSE
}

equalSegment <- function(S.1, S.2){
  flag <- TRUE
  for (i in 1:length(S.1)){
    set.1 <- extractSet(S.1[[i]])
    set.2 <- extractSet(S.2[[i]])
    if (!setequal(set.1,set.2)) flag <- FALSE
  }
  flag
}

replaceSegments <- function(C, Segments.remove, Segments.add){
  C.out <- list()
  cnt <- 1
  for (Seg in C){
    if (!equalSegment(Seg, Segments.remove)){
      C.out[[cnt]] <- Seg
      cnt <- cnt + 1
    }
    else{
      for (seg in Segments.add){
        C.out[[cnt]] <- unlist(seg)
        cnt <- cnt + 1
      }
    }
  }
  C.out
}

mergeTrees <- function(Cohorts, Metrics){
  C.out <- Cohorts[[1]]
  for (l in 2:length(Cohorts)){
    C.base <- C.out
    for (Segment.1 in C.base){
      for (Segment.2 in Cohorts[[l]]){
        Segment.int <- intersectSegment(Segment.1, Segment.2)
        Segment.complement <- complementSegment(Segment.1, Segment.int)
        if (isStatSig(list(Segment.int, Segment.complement), Metrics, c(4,5)) == TRUE){
          C.out <- replaceSegments(C.out, Segment.1, list(Segment.int,Segment.complement))
          Segment.1 <- Segment.complement
        }
      }
    }
  }
  C.out
}

########################### Preparing Example Input Data ##########################
#Tree models in the forms of Cohorts
C.1 <- list(c("0,1", "2,3,4,5,6,7,8,9,10", "1,2,3"),
           c("-5,-4,-3,-2,-1,2,3,4,5,6,7", "0,1", "4,5"))

C.2 <- list(c("0", "2,3,4,5,6", "1,2"),
           c("-5,-4,-3,-2,-1", "0,1", "5"),
           c("2,3,4,5,6,7", "7,8,9,10","3,4"))

C.3 <- list(c("-1,0,1,2,3,4", "5,6", "1,2"),
           c("-5,-4,-3-2", "0,1,2,3,4", "3"),
           c("5,6", "7,8,9", "4"),
           c("7", "10", "5"))


Cohorts <- list(C.1, C.2, C.3)

set.seed(1)
#Metrics Data
size <- 100000
treatment.assignment <- sample(0:1,size, replace = TRUE)
feature.1 <- sample(-5:7,size, replace = TRUE)
feature.2 <- sample(0:10,size, replace = TRUE)
feature.3 <- sample(1:5,size, replace = TRUE)

Metrics <- data.frame(
  feature.1 <- feature.1,
  feature.2 <- feature.2,
  feature.3 <- feature.3,
  metric.1 <- feature.1 * rnorm(size) + feature.3 * rnorm(size) + treatment.assignment * 0.1 + rnorm(size),
  metric.2 <- feature.2 * rnorm(size) + treatment.assignment * 5 + rnorm(size),
  treatment.var <- treatment.assignment
)


########################### Running the Tree Merge Algorithm ##########################

finalCohort <- mergeTrees(Cohorts, Metrics)