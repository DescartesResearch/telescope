#' @author Marwin Zuefle, Andre Bauer

#' @description Removes the anomalies of the time series
#'
#' @title Remove the Anomalies
#' @param rawValues The raw values as vector of the time series without timestamp
#' @param frequency The frequency of the time series
#' @param replace.zeros If TRUE, all zeros will be replaced by the mean of the non-zero neighbors
#' @return The vector of observations without anomalies (anomalies replaced by mean of normal values around)
removeAnomalies <- function(rawValues,frequency,replace.zeros) {

  vals <- rawValues
  # Gets the anomalies
  anom <- AnomalyDetectionVec(x = vals, period = frequency, direction = "pos",plot = FALSE)
  anomPos <- anom$anoms$index
  anomVals <- anom$anoms$anoms
  allPos <- c(1:length(vals))
  if(!is.null(anomPos)) {
    normPos <- allPos[-anomPos]

    for(i in 1:length(anomPos)) {
      # find the next lower and upper non-anomaly neighbor of each anomaly
      if(length(which(normPos>anomPos[i]))==0) {
        lb <- max(normPos[which(normPos<anomPos[i])])
        ub <- max(normPos[which(normPos<anomPos[i])])
      } else if(length(which(normPos<anomPos[i]))==0) {
        lb <- min(normPos[which(normPos>anomPos[i])])
        ub <- min(normPos[which(normPos>anomPos[i])])
      } else {
        lb <- max(normPos[which(normPos<anomPos[i])])
        ub <- min(normPos[which(normPos>anomPos[i])])
      }
      # Interpolates the values to replace the anomaly
      vals[anomPos[i]] <- (vals[lb]+vals[ub])/2
    }
  }

  if(replace.zeros) {
    anomPos.null <- which(vals==0)
    if(length(anomPos.null)>0) {
      normPos.null <- allPos[-anomPos.null]
      # find the next lower and upper non-zero neighbor for each zero value
      for(i in 1:length(anomPos.null)) {
        if(length(which(normPos.null>anomPos.null[i]))==0) {
          lb <- max(normPos.null[which(normPos.null<anomPos.null[i])])
          ub <- max(normPos.null[which(normPos.null<anomPos.null[i])])
        } else if(length(which(normPos.null<anomPos.null[i]))==0) {
          lb <- min(normPos.null[which(normPos.null>anomPos.null[i])])
          ub <- min(normPos.null[which(normPos.null>anomPos.null[i])])
        } else {
          lb <- max(normPos.null[which(normPos.null<anomPos.null[i])])
          ub <- min(normPos.null[which(normPos.null>anomPos.null[i])])
        }
        # Interpolates the values to replace the zeros
        vals[anomPos.null[i]] <- (vals[lb]+vals[ub])/2
      }
    }
  }

  return(vals)
}
