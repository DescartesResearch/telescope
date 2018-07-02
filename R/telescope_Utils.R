#' @author Marwin Zuefle, Andre Bauer

#' @description Saves the vales to a csv-file
#'
#' @title Save as CSV
#' @param values The values to be stored in the csv file
#' @param name The name of the csv file to be created
#' @param csv.path Optional parameter: The path for the saved csv-file. The current workspace by default.
save.csv <- function(values, name, path = '') {
  save_here <- paste(path, Sys.time(), "_", name, ".csv", sep = "")
  save_here <- gsub(" ", "_", save_here)
  save_here <- gsub(":", "-", save_here)
  save_here <- paste("C:", save_here, sep = "")
  write.csv(values, file = save_here)
}


#' @description Extract the requird information (e.g., frequency, values, etc.) of the time series
#'
#' @title Extract time series information
#' @param tvp The time value pair: either vector of raw values or n-by-2 matrix (raw values in second column), or time series
#' @param use.sec.freq Determines if a second frequency shall be used
#' @return the time value pair as vector, the frequency of the data, use.sec.freq and if the data was no time series, the last iteration and the periodigram of the frequncy estimation is also returned
extract.info <- function(tvp, use.sec.freq, debug=FALSE) {

  # If the time value pair is a time series, get the values and the frequency
  if (is.ts(tvp)) {
    frequency <- frequency(tvp)
    use.second.freq <- FALSE
    tvp <- tvp[1:length(tvp)]

    return(
      list(
        "values" = as.vector(tvp),
        "frequency" = frequency,
        "use.second.freq" = use.second.freq
      )
    )

  } else  {
    # If the time value pair is not a time series, estimate frequency and extract values
    tvp <- as.matrix(tvp)
    if (ncol(tvp) == 2) {
      tvp <- tvp[, 2]
    } else if (ncol(tvp) > 2) {
      stop(
        "Input time series has to many columns. Either single column with only raw values or two columns with raw values in second column!"
      )
    }
    freq <-
      calcFrequencyPeriodogram(
        timeValuePair = tvp,
        asInteger = TRUE,
        difFactor = 0.5,
        debug = debug
      )

    frequency <- freq$frequency[1]
    use.second.freq <- use.sec.freq

    return(
      list(
        "values" = as.vector(tvp),
        "frequency" = frequency,
        "use.second.freq" = use.second.freq,
        "lastIterfreq" = freq$lastIterfreq,
        "pgram" = freq$pgram
      )
    )
  }



}

#' @description Forecasts the season part of the time series
#'
#' @title Forecasting Season
#' @param tvp The time value pair as vector
#' @param stlTrain The decomposition of tvp
#' @param horizon The forecast horizon
#' @return The seasonal pattern for the forecast horizon
forecast.season <- function(tvp, stlTrain, horizon) {
  # Total length = history + forecast
  total.length <- length(tvp$values) + horizon
  # As the season is per defintion recurring, repeat season over the horizon
  fullper <- as.integer(total.length / tvp$frequency)
  rest <- total.length - (fullper * tvp$frequency)
  fcSeason <- rep(stlTrain$time.series[1:tvp$frequency, 1], fullper)
  if (rest > 0) {
    fcSeason <- c(fcSeason, stlTrain$time.series[1:rest, 1])
  }
  return(fcSeason)
}

#' @description Forecats the trend part of the time series
#'
#' @title Forecasting Trend
#' @param model The model of the trend, either exponential or linear
#' @param tsTrainTrend The trend component of the time series
#' @param frequency The frequency of the time series
#' @param horizon The forecast horizon
#' @return The trend pattern for the forecast horizon
forecast.trend <- function(model, tsTrainTrend, frequency, horizon) {
  # If trend is exponential, log time series
  if (model == "exp") {
      tsTrainLOG <- ts(log(tsTrainTrend), frequency = frequency)
      fArimaLOG <- doArima(tsTrainLOG, FALSE)
      fcArimaLOG <- forecast(fArimaLOG, h = horizon)$mean
      print("exponential Trend detected!")
      fcArima <- exp(fcArimaLOG)
    } else {
      fArima <- doArima(tsTrainTrend, FALSE)
      fcArima <- forecast(fArima, h = horizon)$mean
    }
  return(fcArima)
}

#' @description Checks if the time series has a significant high remainder
#'
#' @title Checking Remainder
#' @param tvp tvp The time value pair as vector
#' @param stlRemainder The remainder part of the deccomposition of the tvp
#' @param use.log A flag if log was used for tvp
#' @param sig.dif.factor The threshold that is to exceed for having a high remainder
#' @return If the the time series has a high remainder compared to the threshold
has.highRemainder <- function(tvp, stlRemainder, use.log, sig.dif.factor) {
  # Calculates the IQR of the remainder and original time series
  if (use.log) {
    remainder.quantiles.log <- quantile(stlRemainder)
    # range from 25% to 75% quantile
    remainder.quantiles.range <- exp(remainder.quantiles.log[4] - remainder.quantiles.log[2])
    tvp.quantiles <- quantile(tvp)
    tvp.quantiles.range <- exp(tvp.quantiles[4] - tvp.quantiles[2])
  } else {
    remainder.quantiles <- quantile(stlRemainder)
    # range from 25% to 75% quantile
    remainder.quantiles.range <-
      remainder.quantiles[4] - remainder.quantiles[2]
    tvp.quantiles <- quantile(tvp)
    tvp.quantiles.range <- tvp.quantiles[4] - tvp.quantiles[2]
  }
  # If the remainder IQR has a higher propotion than the threshold
  if (sig.dif.factor * tvp.quantiles.range < remainder.quantiles.range) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @description Performs the ARIMA forecast of the timeseries.
#'
#' @title Apply Arima
#' @param ts The timeseries.
#' @param season If false, only non-seasonal models will be fitted
#' @return The Arima Forecast of \code{ts}.
doArima <- function(ts, season = TRUE){
  bool <- is.ts(ts)
  if (is.ts(ts)){
    fc <- auto.arima(ts, stepwise = TRUE, seasonal = season)
    return(fc)
  }
  return(NULL)
}


#' @description Performs the ANN forecast of the timeseries.
#'
#' @title Apply ANN
#' @param ts The timeseries.
#' @param rep The amount of repeats
#' @return The ANN Forecast
doANN <- function(myts,rep = 20){
  result <- tryCatch({
    nnetar(myts, repeats = rep)
  }, error = function(e){
    print(paste("Some error doing ANN, probably stl: ", e, sep="" ))
    tryCatch({
      nnetar(myts, repeats = rep, p=1)
    }, error = function(e){
      print(paste("Other unknown error: ", e, sep = ""))
      snaive(x = myts, h = length(myts)/4)
    })
  })
  return(result)
}

#' @description Creates a tag containing the MASE of the forecasting methods for the timeseries.
#'
#' @title Compute MASE
#' @param forecast The forecasted values.
#' @param train The 'historical' data used for forecasting.
#' @param test The 'future' data used for finding MASE.
#' @param plot Boolean indicating whether the forecast should be plotted.
#' @return MASE between forecast and real data
computeMASE <- function(forecast, train, test, plot){
  if(plot){
    plot(1:length(test), test,type="l",col="black", main = 'History (black) and Model (red)', xlab = 'Index', ylab = 'Observation')
    lines(1:length(forecast), forecast, type = "l", col="red")
  }

  forecast <- as.vector(forecast)
  train <- as.vector(train)
  test <- as.vector(test)

  # calculate scaling factor
  test <- append(test, train[length(train)], after = 0)
  n <- length(test)
  if(n == 1){
    stop('Computing MASE: Test vector of length 0 is invalid.')
  }
  scalingFactor <- sum(abs(test[2:n] - test[1:(n-1)])) / (n-1)

  # Avoding to divide by zero
  if(scalingFactor==0) {
    scalingFactor<-0.00001
  }

  # calculate MASE
  et <- abs(test[2:length(test)]-forecast)
  qt <- et/scalingFactor
  meanMASE <- mean(abs(qt))

  return(meanMASE)
}

#' @description MASE for a naive forecast taking the last observation for the whole forecast
#'
#' @title Compute MASE for same value
#' @param forecast The forecasted values.
#' @param train The 'historical' data used for forecasting.
#' @param test The 'future' data used for finding MASE.
#' @param plot Boolean indicating whether the forecast should be plotted.
#' @return MASE between forecast and real data
computeMASEsameValue <- function(forecast, train, test, plot){
  if(plot){
    plot(1:length(test), test,type="l",col="black", main = deparse(substitute(forecast)))
    lines(1:length(forecast), forecast, type = "l", col="red")
  }

  forecast <- as.vector(forecast)
  train <- as.vector(train)
  test <- as.vector(test)

  # calculate scaling factor
  n <- length(test)
  if(n == 1){
    stop('Computing MASE: Test vector of length 0 is invalid.')
  }
  lastObservation <- tail(train,1)
  maseForecast <- rep(lastObservation,n)
  scalingFactor <- sum(abs(maseForecast - test[1:n])) / n

  # Avoding to divide by zero
  if(scalingFactor==0) {
    scalingFactor<-0.00001
    print("scaling factor is 0")
  }

  # calculate MASE
  et <- abs(test-forecast)
  qt <- et/scalingFactor
  meanMASE <- mean(abs(qt))

  return(meanMASE)
}
