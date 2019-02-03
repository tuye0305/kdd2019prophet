# We are trying to merging multiple trees into one single cohort definition
library(MASS)

########################### Helper Functions ##########################
calculatePvalue <- function(df, metricIdx) {
  treatmentIds <- which(df[[treatmentVar]] == 1)
  controlIds <- which(df[[treatmentVar]] == 0)
  y1 <- mean(df[,metricIdx][treatmentIds])
  y0 <- mean(df[,metricIdx][controlIds])
  n1 <- as.numeric(nrow(df[treatmentIds, ]))
  n0 <- as.numeric(nrow(df[controlIds, ]))
  varY1 <- var(df[,metricIdx][treatmentIds]) / n1
  varY0 <- var(df[,metricIdx][controlIds]) / n0
  delta <- y1 - y0

  varianceDelta <- varY1 + varY0
  tStatsDelta <- delta / sqrt(varianceDelta)
  pValueDelta <- 2 * pnorm(- abs(tStatsDelta))
  pValueDelta
}

extractSet <- function(str){
  if (str == "NULL") set <- NULL
  else set <- as.numeric(trim(unlist(strsplit(str, ","))))
  set
}

#take intersection of two Segments by evaluate the sets per feature
intersectSegment <- function(S1, S2){
  Segment_out <- list()
  for (i in 1:length(S1)){
    s1 <- S1[[i]]
    s2 <- S2[[i]]
    set1 <- extractSet(s1)
    set2 <- extractSet(s2)

    set_intersect <- intersect(set1, set2)
    if (is.null(set_intersect) || length(set_intersect) == 0) Segment_out[[i]] <- "NULL"
    else Segment_out[[i]] <- paste(noquote(set_intersect), collapse=",")
  }
  Segment_out
}

#Get (S1 - S2) by evaluate the sets per feature
complementSegment <- function(S1, S2){
  Segment_out <- list()
  for (i in 1:length(S1)){
    s1 <- S1[[i]]
    s2 <- S2[[i]]
    set1 <- extractSet(s1)
    set2 <- extractSet(s2)

    set_complement <- setdiff(set1, set2)
    if (is.null(set_complement) || length(set_complement) == 0) Segment_out[[i]] <- "NULL"
    else Segment_out[[i]] <- paste(noquote(set_complement), collapse=",")
  }
  Segment_out
}

buildSubtableSet <- function(Segments, Metrics){
  subtableSet <- c()
  for (segment in Segments){
    subtable <- Metrics
    for (i in 1: length(segment)){
      s <- segment[[i]]
      set1 <- extractSet(s)
      subtable <- subtable[which(subtable[,i] %in% set1), ]
    }
    subtableSet <- append(subtableSet, subtable)
  }
}

isStatSig <- function(Segments, Metrics, metricsIdxList){
  subtableSet <- buildSubtableSet(Segments, Metrics)
  for (k in metricsIdxList){
    flag <- TRUE
    for (subtable in subtableSet){
      print(nrow(subtable) == 0)
      if (nrow(subtable) == 0 || calculatePvalue(subtable, k) > 0.5){
        flag <- FALSE
      }
    }
    if (flag) flag
  }
  FALSE
}

equalSegment <- function(S1, S2){
  flag <- TRUE
  for (i in 1:length(S1)){
    s1 <- S1[[i]]
    s2 <- S2[[i]]
    set1 <- extractSet(s1)
    set2 <- extractSet(s2)
    if (!setequal(set1,set2)) flag <- FALSE
  }
}

replaceSegments <- function(C, Segments_remove, Segments_add){
  C_out <- list()
  cnt <- 1
  for (Seg in C){
    if (!equalSegment(Seg, Segments_remove)){
      C_out[[cnt]] <- Seg
      cnt <- cnt + 1
    }
    else{
      C_out <- append(C_out, Segments_add)
      cnt <- cnt + length(Segments_add)
    }
  }
  C_out
}

mergeTrees <- function(Cohorts, Metrics){
  C_out <- Cohorts[[1]]
  for (l in 2:length(Cohorts)){
    C_base <- C_out
    for (Segment1 in C_base){
      for (Segment2 in Cohorts[[l]]){
        Segment_int <- intersectSegment(Segment1, Segment2)
        Segment_complement <- complementSegment(Segment1, Segment_int)
        print(isStatSig(c(Segment_int, Segment_complement), Metrics, c(4,5)))
        if (isStatSig(c(Segment_int, Segment_complement), Metrics, c(4,5))){
          C_out <- replaceSegments(C_out, Segment1, c(Segment_int,Segment_complement))
          Segment1 <- Segment_complement
        }
      }
    }
  }
  C_out
}

########################### Preparing Example Input Data ##########################
#Tree models in the forms of Cohorts
C1 <- list(c("NULL", "0,1,2,3,4,5,6,7,8,9,10", "1,2,3"),
           c("-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7", "NULL", "4,5"))

C2 <- list(c("0,1", "2,3,4,5,6,7", "5"),
           c("-5,-4,-3,-2,-1", "0,1", "1,2"),
           c("2,3,4,5,6,7", "8,9,10","3,4"))

C3 <- list(c("-1,0,1,2,3,4", "5,6", "1,2"),
           c("-5,-4,-3-2", "0,1,2,3,4", "3"),
           c("5,6", "7,8,9", "4"),
           c("7", "10", "5"))


Cohorts <- list(C1, C2, C3)

set.seed(1)
#Metrics Data
size <- 10000
treatmentAssignment <- sample(0:1,size, replace = TRUE)
feature1 <- sample(-5:7,size, replace = TRUE)
feature2 <- sample(0:10,size, replace = TRUE)
feature3 <- sample(1:5,size, replace = TRUE)

Metrics <- data.frame(
  feature1 <- feature1,
  feature2 <- feature2,
  feature3 <- feature3,
  metric1 <- feature1 * rnorm(size) + feature3 * rnorm(size) + treatmentAssignment * 0.1 + rnorm(size),
  metric2 <- feature2 * rnorm(size) + treatmentAssignment * 5 + rnorm(size),
  treatmentVar <- treatmentAssignment
)


########################### Running the Tree Merge Algorithm ##########################

finalCohort <- mergeTrees(Cohorts, Metrics)