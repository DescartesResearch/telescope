#' @author Andre Bauer, Marwin Zuefle

#' @description Forecasts a given univariate time series in a hybrid manner and based on time series decomposition
#'
#' @title Perform the Forecast
#' @param tvp The time value pair: either vector of raw values or n-by-2 matrix (raw values in second column), or time series
#' @param horizon The number of values that should be forecast
#' @param rec_model Optional parameter: 
#' @param natural Optional parameter: A flag indicating wheter only natural frequencies (e.g., daily, hourly, ...) or all found frequencies shall be considered.
#' @param boxcox Optional parameter: A flag indicating if the Box-Cox transofrmation should be performed. It is not recommend to disable the transformation. TRUE by default.
#' @param doAnomDet  Optional parameter: Boolean whether anomaly detection shall be used. FALSE by default
#' @param replace.zeros  Optional parameter: If TRUE, all zeros will be replaced by the mean of the non-zero neighbors. TRUE by default
#' @param use.indicators  Optional parameter: If TRUE, additional information (e.g. a flag wheter there is a high remainder) will be returned. TRUE by default
#' @param save_fc  Optional parameter: Boolean wheter the forecast shall be saved as csv. FALSE by default
#' @param csv.path Optional parameter: The path for the saved csv-file. The current workspace by default.
#' @param csv.name Optional parameter: The name of the saved csvfile. Telescope by default.
#' @param debug Optional parameter: If TRUE, debugging information will be displayed. FALSE by default
#' @return The forecast of the input data
#' @examples
#' telescope.forecast(forecast::taylor, horizon=10)
#' @export
telescope.forecast <- function(tvp, horizon, rec_model=NULL, natural=TRUE, boxcox = TRUE, doAnomDet = FALSE, replace.zeros = TRUE, use.indicators = TRUE, save_fc = FALSE, csv.path = '', csv.name = "Telescope", debug = FALSE, plot = TRUE) {
  
    if(anyNA(tvp)) {
      stop("Telescope does not support NA values, only numeric.")
    }
  
    if(!is.null(rec_model) && !is(rec_model, 'telescope.recommendation')){
      stop('The model must be of class: telescope.recommendation!')
    }
  
    # If the time value pair is not a time series, estimate frequency and extract values
    if(!is.ts(tvp)){
      tvp <- extract.ts(tvp, natural, debug)
      print(paste("Found frequency:", frequency(tvp)))
    }
  
    # Alternative for non-seasonal time series 
    # STL requires at least two full periods
    if(frequency(tvp) < 2 || length(tvp) <= 2*frequency(tvp)){
      
      if(boxcox){
        lambda <- BoxCox.lambda(tvp, lower = 0, upper = 1)
        print(paste("Found Lambda for BoxCox:", lambda))
        tvp <- BoxCox(tvp, lambda)
      }  
      
      arima.model <- doArima(tvp,frequency(tvp)>1)
      values <- forecast(arima.model,h=horizon)$mean
      arima.fit <- ts(arima.model$fitted, frequency = frequency(tvp))
      
      if(boxcox){
        values <- InvBoxCox(values, lambda)
        tvp <- InvBoxCox(tvp, lambda)
        arima.fit <- InvBoxCox(arima.fit, lambda)
      }
      
      if(frequency(tvp) < 2){
        print("Switch to fallback as time series has a frequency of 1")
      } else {
        print("Switch to fallback as STL requires at least 2 full periods")
      }
      
      output.mean <- values
      output.x <- as.vector(tvp)
      output.residuals <-
        output.x - arima.fit
      output.method <- "Telescope"
      output.accuracy <- accuracy(arima.fit, tvp)
      output.fitted <- arima.fit
      
      
      output <- list(mean=output.mean, x=output.x, residuals=output.residuals, method=output.method, 
                     fitted=output.fitted)
      
      print(output.accuracy)
      
  
      return(structure(output, class = 'forecast'))
      
    }
  

    startTime <- Sys.time()
    
    
    # Remove all Anomalies on the raw time series first
    if(frequency(tvp) > 10 && doAnomDet) {
      tvp <- ts(removeAnomalies(as.vector(tvp),frequency = frequency(tvp), replace.zeros = replace.zeros),frequency = frequency(tvp))
    }
    
    
    prep <- preprocessing(tvp,natural,boxcox)
    
    tvp <- prep$tvp
    
    tvp.mul <- msts(tvp, seasonal.periods = prep$freqs)
    
    tvp.stl <- stl(tvp,s.window = "periodic",t.window=length(tvp)/2)
    fourier.terms <- fourier(tvp.mul, K = rep(1,length(prep$freqs)))
    
    # check if the time series has a high remainder
    high_remainder <- has.highRemainder(tvp, tvp.stl$time.series[,3], sig.dif.factor=0.5)
    if (high_remainder) {
      print("-------------- ATTENTION: High remainder in STL --------------")
    }
    
    train <- cbind(tvp.stl$time.series[,1:2], fourier.terms)
    colnames(train) <- c('Season', 'Trend', colnames(fourier.terms))
    
    model <- fittingModels(tvp.stl,frequency = frequency(tvp),difFactor = 1.5, debug = FALSE)
    
    if (model$risky_trend_model) {
      print("-------------- ATTENTION: risky trend estimation --------------")
    }
    
    
    # Forecast the features
    fc.trend <- forecast.trend(model$trendmodel,train[,2], frequency(tvp), horizon)
    fc.season <- forecast.season(length(tvp)+horizon, train[,1], frequency(tvp), length(tvp))
    fc.fourier.terms <- fourier(tvp.mul, K = rep(1,length(prep$freqs)), h = horizon)
    
    fc.features <- cbind(fc.season, fc.trend, fc.fourier.terms)
    colnames(fc.features) <- c('Season', 'Trend', colnames(fc.fourier.terms))
   
    xgbcov <- as.matrix(train[,-2])
    xgblabel <- as.vector(tvp - train[,2])
    booster <- "gblinear"
    testcov <- as.matrix(fc.features[,-2])
    
   
    
    if(is.null(rec_model)){
      fXGB <- doXGB.train(myts = xgblabel, cov = xgbcov, booster = booster, verbose = 0)
    } else {
      
      # gets the best machine learning method for the time series
      method <- consultrecommender(tvp=tvp,tvp.stl=tvp.stl,model=rec_model)
      
      
      
      print(paste(method, "is selected."))
      switch(method,
             # "Catboost"= {},
             "Cubist" = {
               
               fXGB <- cubist(x=xgbcov, y=xgblabel)
             },
             "Evtree" = {
               
               data <- as.data.frame(cbind(xgbcov,xgblabel))
               colnames(data) <- c(colnames(xgbcov), 'Target')
               if(nrow(data) <= 20){
                 fXGB <- evtree(Target ~ ., data = data, control = evtree.control(minsplit = 2L, minbucket = 1L))
               } else {
                 fXGB <- evtree(Target ~ ., data = data)
               }
             },
             "Nnetar"={
               
               fXGB <- nnetar(y=ts(xgblabel, frequency = frequency(tvp)),xreg=xgbcov, MaxNWts=2000)
             },
             "RF"= {
               
               fXGB <- randomForest(y = xgblabel, x = xgbcov)
             },
             "Rpart"={
               
               data <- as.data.frame(cbind(xgbcov,xgblabel))
               colnames(data) <- c(colnames(xgbcov), 'Target')
               fXGB <- rpart(Target ~ ., data = data, method="anova",control = rpart.control(minsplit = 2, maxdepth = 30, cp = 0.000001))
             },
             "SVR"= {
               
               fXGB <- svm(y = xgblabel, x = xgbcov)
             },
             "XGBoost"= {
               
               fXGB <- doXGB.train(myts = xgblabel, cov = xgbcov, booster = booster, verbose = 0)
             }
      )
      
      
    }
    
    
   
    # Unify names
    fXGB$feature_names <- colnames(xgbcov)
    
    
    if(horizon == 1){
      testcov <- t(testcov)
    }
    
    
    if(is.null(rec_model)){
      predXGB <- predict(fXGB, testcov)
    } else {
      switch(method,
             # "Catboost"= {},
             "Cubist" = {
               
               predXGB <- predict(fXGB, testcov)
             },
             "Evtree" = {
               
               predXGB <- predict(fXGB, as.data.frame(testcov))
             },
             "Nnetar"={
               
               predXGB <- as.vector(forecast(fXGB, h=nrow(testcov), xreg=testcov)$mean)
             },
             "RF"= {
               
               predXGB <- predict(fXGB, testcov)
             },
             "Rpart"={
               
               predXGB <- predict(fXGB, as.data.frame(testcov))
             },
             "SVR"= {
               
               predXGB <- predict(fXGB, testcov)
             },
             "XGBoost"= {
               
               predXGB <- predict(fXGB, testcov)
             }
      )
    }
    
    predXGB <- predXGB + fc.trend
    
    
    if (prep$minValue2 <= 0) {
      predXGB <- predXGB - abs(prep$minValue2) - 1
    }
    
    if(boxcox){
      # Invert BoxCox transformation
      predXGB <- InvBoxCox(predXGB, prep$lambda)
    }
    
    # Undo adjustment to positive values
    if (prep$minValue <= 0) {
      predXGB <- predXGB - abs(prep$minValue) - 1
    }
    
    if (save_fc) {
      save.csv(values = predXGB, name = csv.name, path = csv.path)
    }
    
    endTime <- Sys.time()
    if(debug){
      plot(tvp.stl)
      print(paste("Time elapsed for the whole forecast:", difftime(endTime, startTime, units = "secs")))
    }
    
    par(mfrow = c(2, 1))
    
    # Get model of the history
    if(!is.null(rec_model) && method=="Nnetar"){
      xgb.model <- as.vector(fXGB$fitted)
      xgb.model[which(is.na(xgb.model))] <- 0
    } else {
      xgb.model <- predict(fXGB, xgbcov)
    }
    
    xgb.model <- xgb.model + train[,2]

    
    
    if (prep$minValue2 <= 0) {
      xgb.model <- xgb.model - abs(prep$minValue2) - 1
      tvp <- tvp - abs(prep$minValue2) - 1
    }
    
    if(boxcox){
      # Invert BoxCox transformation
      xgb.model <- InvBoxCox(xgb.model, prep$lambda)
      tvp <- InvBoxCox(tvp, prep$lambda)
    }  
    
    
    # Undo adjustment to positive values
    if (prep$minValue <= 0) {
      xgb.model <- xgb.model - abs(prep$minValue) - 1
      tvp <- tvp - abs(prep$minValue) - 1
    }
    
    # Calculates the accuracies of the trained model
    accuracyXGB <- accuracy(xgb.model, tvp)
    print(accuracyXGB)
    
    # Build the time series with history and forecast
    fcOnly <- ts(predXGB, frequency = frequency(tvp))
    
    # Plot the model and the time series
    if(plot) {
      y.min <- min(min(tvp[-1]),min(xgb.model[-1]))
      y.max <- max(max(tvp[-1]),max(xgb.model[-1]))
      plot(1:length(tvp[-1]), tvp[-1],type="l",col="black", main = 'History (black) and Model (red)', xlab = 'Index', ylab = 'Observation', xlim = c(0, length(tvp)+horizon), ylim = c(y.min, y.max) )
      lines(1:length(xgb.model[-1]), xgb.model[-1], type = "l", col="red")
    }
    
    
    # Plot the forecasted time series and the original time series
    if(plot) {
      y.min <- min(min(fcOnly),min(tvp))
      y.max <- max(max(fcOnly),max(tvp))
      plot(length(tvp):(length(tvp)+horizon), c(tvp[length(tvp)],as.vector(fcOnly)),type = 'l',col="red",xlab = 'Index', ylab = 'Observation', main = 'History (black) and Forecast (red)', xlim = c(0, length(tvp)+horizon), ylim = c(y.min, y.max))
      lines(1:length(tvp), tvp)
    }
    
    # Collect information for output
    output.mean <- fcOnly
    output.x <- as.vector(tvp)
    output.residuals <-
      output.x - ts(xgb.model, frequency = frequency(tvp))
    output.method <- "Telescope"
    output.accuracy <- accuracyXGB
    output.fitted <- xgb.model
    
    
    if(use.indicators) {
      output.risky.trend.model <- model$risky_trend_model
      output.high.stl.remainder <- high_remainder
      output <- list(mean=output.mean, x=output.x, residuals=output.residuals, method=output.method, 
                     fitted=output.fitted, riskytrend=output.risky.trend.model, highresiduals=output.high.stl.remainder)
    } else {
      output <- list(mean=output.mean, x=output.x, residuals=output.residuals, method=output.method, 
                     fitted=output.fitted)
    }
    
    return(structure(output, class = 'forecast'))
    
}



#' @description Preprocesses the time series and retrieves dominant frequncies
#'
#' @title Perform the Forecast
#' @param tvp The time value pair: a time series
#' @param natural Optional parameter: A flag indicating wheter only natural frequencies (e.g., daily, hourly, ...) or all found frequencies shall be considered.
#' @param boxcox Optional parameter: A flag indicating if the Box-Cox transofrmation should be performed. It is not recommend to disable the transformation. TRUE by default.
#' @return The adjusted time series, the dominant frequencies, and adjustment parameters
preprocessing <- function(tvp,natural=TRUE,boxcox=TRUE){
  # Calculating the most dominant frequencies
  freq <- calcFrequencyPeriodogram(timeValuePair = as.vector(tvp), asInteger = TRUE, difFactor = 0.5, debug = FALSE)
  spec <- freq$pgram$spec[order(freq$pgram$spec, decreasing = TRUE)]/max(freq$pgram$spec)
  
  
  if(natural){
    freqs <- c()
    
    # Get frequencies that are significant
    for(i in 1:(length(freq$pgram$freq))){
      freqs <- c(freqs ,calcFrequencyPeriodogram(timeValuePair = as.vector(tvp), asInteger = TRUE, difFactor = 0.5,maxIters = 10,ithBest = i, PGramTvp = freq$pgram,debug=FALSE)$frequency)
      freqs <- unique(freqs)
      if(spec[i] < 0.5){
        break
      }
    }
  } else {
    freqs <- 1/freq$pgram$freq[order(freq$pgram$spec, decreasing = TRUE)]
    freqs <- freqs[1:(which(spec<0.5)[1]-1)]
  }
  
  
  freqs <- unique(c(freqs, frequency(tvp)))
  if(1 %in% freqs){
    freqs <- freqs[-which(freqs==1)] 
  }
  
  # get the minimum value to shift all observations to real positive values
  minValue <- min(tvp)
  if (minValue <= 0) {
    tvp <- tvp + abs(minValue) + 1
  }
  
  lambda <- NULL
  
  if(boxcox){
    # Calculating lambda for BoxCox
    lambda <- BoxCox.lambda(tvp, lower = 0, upper = 1)
    if(lambda < 0.1) lambda = 0
    print(paste("Found Lambda for BoxCox:", lambda))
    
    # BoxCox Transformation
    tvp <- BoxCox(tvp, lambda)
  }
  
  
  minValue2 <- min(tvp)
  if (minValue2 <= 0) {
    tvp <- tvp + abs(minValue2) + 1
  }
  
  return(list('tvp'=tvp,'lambda'=lambda, 'freqs'=freqs, 'minValue'=minValue, 'minValue2'=minValue2))
}

