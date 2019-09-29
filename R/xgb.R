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

  nsample <- as.integer(0.8 * nrow(combined))

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
  watchlist <- list(train=dtrain, val=dtest)

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


