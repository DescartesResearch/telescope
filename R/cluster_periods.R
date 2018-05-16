#' @author Marwin Zuefle, Andre Bauer

#' @description Calculates the clusters for the periods of the time series regarding mean, variance and range
#'
#' @title Calculates the Clusters for Periods
#' @param timeseries The time series to split
#' @param frequency The determined frequency of the time series
#' @param doAnomDet Boolean whether anomaly detection shall be used for clustering
#' @param replace.zeros If TRUE, all zeros will be replaced by the mean of the non-zero neighbors
#' @param debug Optional parameter: If TRUE, debugging information will be displayed. FALSE by default
#' @return The cluster label for each period
calcClustersForPeriods <- function(timeseries,frequency,doAnomDet,replace.zeros,debug=FALSE) {

  min.sil <- 0.66

  if(doAnomDet) {
    timeseries <- removeAnomalies(as.vector(timeseries),frequency,replace.zeros = replace.zeros)
  }

  set.seed(750)

  # Calculates the number of full periods
  lenSeq <- length(timeseries) + frequency - 1
  times <- seq(frequency,lenSeq)/frequency
  amountFullPeriods <- as.integer(length(timeseries)/frequency)

  means <- c()
  vars <- c()
  ranges <- c()

  # Calculates mean, variance and range over each period
  i <- 0
  while(i<amountFullPeriods) {
    start <- i*frequency+1
    end <- (i+1)*frequency
    means[(i+1)] <- mean(timeseries[start:end])
    vars[(i+1)] <- var(timeseries[start:end])
    ranges[(i+1)] <- max(timeseries[start:end])-min(timeseries[start:end])
    i <- i + 1
  }

  # perform k-means clustering based on mean, variance and range of the data
  ranges <- ranges/max(ranges)
  means <- means/max(means)
  vars <- vars/max(vars)
  features <- cbind(vars,ranges)
  if(debug){
    par(mfrow=c(3,1))
    print(features)
  }

  # Perform kmeans
  clk <- kmeans(features, centers = 2)


  # determine dissimilarity
  dissE <- daisy(features)
  dissE2 <- dissE^2

  # 2 centers
  sil <- silhouette(clk$cluster,dissE2)
  sum.sil <- summary(sil)
  vals.sil <- sum.sil$clus.avg.widths

  # if the silhouette of at least one cluster is below the min.sil threshold, only a vector of ones is created
  # otherwise, a vector containing the cluster labels is created
  if(length(which(vals.sil<min.sil))>0) {
    clusters <- rep(1,length(clk$cluster))
  } else {
    clusters <- clk$cluster
  }

  ret <- cbind(means,vars,ranges,clusters)

  return(ret)
}

#' @description Forecasts the clusters for the forecast periods
#'
#' @title Forecast the Clusters
#' @param clusters The clusters found for the periods
#' @param freuqency The determined frequency
#' @param timeseries The time series
#' @param reps The amount of repeats for ANN
#' @param debug Optional parameter: If TRUE, debugging information will be displayed. FALSE by default
#' @return All cluster labels (history and forecast)
forecastClusters <- function(clusters,frequency,timeseries,reps,horizon,debug=FALSE) {

  # split time series in test and training
  clusterTrain <- ts(clusters)
  total.length <- length(timeseries) + horizon

  fullper <- as.integer(total.length/frequency)

  # Required number of new clusters
  reqCluster <- fullper - length(clusters) + 1

  # perform ANN forecast
  fAnn <- doANN(clusterTrain,reps)
  fcAnn <- forecast(fAnn, h = reqCluster)$mean

  fcCluster <- round(fcAnn)

  # Builds vector of all cluster labels
  forecastCluster <- c(clusters,fcCluster)
  if(debug){
    print(forecastCluster)
  }

  # Labels each datapoint with associated cluster
  clusterLabel <- c()
  for(i in 1:(length(forecastCluster)-1)) {
    clusterLabel <- c(clusterLabel,rep(forecastCluster[i],frequency))
  }

  rest <- total.length-(fullper*frequency)
  clusterLabel <- c(clusterLabel,rep(forecastCluster[length(forecastCluster)],rest))

  return(clusterLabel)
}


