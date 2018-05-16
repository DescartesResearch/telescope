#' @author Marwin Zuefle, Andre Bauer

#' @description Trains XGBoost regarding the time series and its covariates
#'
#' @title XGBoost Model Training
#' @param myts The training part of the time series
#' @param cov The covariates for the training
#' @param booster Select a boosting method: gblinear or gbtree
#' @param verbose Set to 1 to print information of xgboost performance; 0 for no prints
#' @return model The XGBoost model
doXGB.train <- function(myts, cov, booster, verbose) {
  combined <- cbind(myts,cov)
  feature.columns <- c(2:ncol(combined))

  set.seed(200)

  nsample <- as.integer(0.2 * nrow(combined))

  # Create a vector of nsample sample indices between 1 and the amount of rows of the matrix combined
  h <- sample(nrow(combined),nsample)

  # In case of only one feature, convert the resulting vector to a matrix
  # Create the special kind of matrices for xgboost (training and validation)
  if(length(feature.columns)==1) {
    dtrain <- xgb.DMatrix(data = as.matrix(combined[h,feature.columns]), label = myts[h])
    dtest <- xgb.DMatrix(data = as.matrix(combined[-h,feature.columns]), label = myts[-h])
  } else {
    dtrain <- xgb.DMatrix(data = combined[h,feature.columns], label = myts[h])
    dtest <- xgb.DMatrix(data = combined[-h,feature.columns], label = myts[-h])
  }

  # Create a watchlist using a validation and a training set to prevent xgboost from overfitting
  watchlist <- list(val=dtest, train=dtrain)

  param <- list(  objective           = "reg:linear",
                  booster             = booster,
                  eta                 = 0.1,
                  max_depth           = 5,
                  subsample           = 0.8,
                  min_child_weight    = 1,
                  num_parallel_tree   = 2
  )

  # Builds model
  model <- xgb.train(params = param, data = dtrain, nthread = 2, nrounds = 500,
                     early_stop_rounds = 50, verbose = verbose, watchlist = watchlist,
                     maximize = FALSE
  )
  return(model)
}

#' @description Estimates the booster for XGBoost based on if the time series has a significant trend
#'
#' @title Estimating the boosting method
#' @param stl.decomp The STL decomposition of a time series
#' @param lower.percentile Optional parameter: The lower percentile for the IPR. 0.05 by default
#' @param upper.percentile Optional parameter: The upper percentile for the IPR. 0.95 by default
#' @param threshold Optional parameter: Threshold what the propotion of trend component has to exceed to be trendy. 0.33 by default
#' @return True if the time series has a significant trend
estimateBooster <- function(stl.decomp, lower.percentile = 0.05, upper.percentile = 0.95, threshold = 0.33) {

  # Calculates the IPR of the noise as noise may have outliers
  range.noise <- IPR(stl.decomp$time.series[,3], lower.percentile, upper.percentile)
  # Calculates the range
  range.trend <- range(stl.decomp$time.series[,2])[2] - range(stl.decomp$time.series[,2])[1]
  range.season <- range(stl.decomp$time.series[,1])[2] - range(stl.decomp$time.series[,1])[1]

  # Calculates the propotion of each component
  range.vec <- c(range.noise, range.trend, range.season)
  range.aggr <- sum(range.vec)
  range.ratios <- range.vec / range.aggr

  # If trend component has a higher propotion than the threshold
  if(range.ratios[2] >= threshold) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}

#' @description computes the interpercentile range of vector x from percentile lower.percentile to percentile upper.percentile
#'
#' @title Inter Percentile Range
#' @param x The vector to compute interpercentile range of
#' @param lower.percentile The lower percentile
#' @param upper.percentile The upper percentile
#' @return The range between the lower and upper percentile
IPR <- function(x, lower.percentile, upper.percentile) {
  q <- quantile(x, c(lower.percentile, upper.percentile))
  y <- q[2] - q[1]
  return(as.double(y))
}
