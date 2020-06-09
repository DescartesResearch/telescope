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
  write.table(values, file = save_here, col.names = FALSE, row.names = FALSE)
}


#' @description Extract the requird information (e.g., frequency, values) of the time series
#' 
#' @title Extract time series information
#' @param tvp The time value pair: either vector of raw values or n-by-2 matrix (raw values in second column), or time series
#' @param natural Optional parameter: A flag indicating wheter only natural frequencies (e.g., daily, hourly, ...) or all found frequencies shall be considered.
#' @param debug Optional parameter: If TRUE, debugging information will be displayed. FALSE by default
#' @return the time value pair as vector, the frequency of the data, use.sec.freq and if the data was no time series, the last iteration and the periodigram of the frequncy estimation is also returned
extract.ts <- function(tvp, natural=TRUE, debug=FALSE) {

    if (is.matrix(tvp) && ncol(tvp) == 2) {
      tvp <- tvp[, 2]
    } else if (is.matrix(tvp) && ncol(tvp) > 2) {
      stop("Input time series has to many columns. Either single column with only raw values or two columns with raw values in second column!")
    }
    freq <- calcFrequencyPeriodogram(timeValuePair = tvp, asInteger = TRUE, difFactor = 0.5, debug = debug)
    if(natural){
      tvp <- ts(tvp, frequency = freq$frequency[1])
    } else {
      tvp <- ts(tvp, frequency = 1/freq$pgram$freq[order(freq$pgram$spec, decreasing = TRUE)][1])
    }
    

    return(tvp)
    
}

#' @description Forecasts the season part of the time series
#'
#' @title Forecasting Season
#' @param total.length history + horizon
#' @param season The seasonal pattern from stl
#' @param frequency The frequency of the timeseries
#' @param hist.length The length of the history
#' @return The seasonal pattern for the forecast horizon
forecast.season <- function(total.length, season, frequency, hist.length) {
  # As the season is per defintion recurring, repeat season over the horizon
  fullper <- as.integer(total.length / frequency)
  rest <- total.length - (fullper * frequency)
  fcSeason <- rep(season[1:frequency], fullper)
  if (rest > 0) {
    fcSeason <- c(fcSeason, season[1:rest])
  }
  return(fcSeason[(hist.length+1):total.length])
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
    minValue <- min(tsTrainTrend)
    if (minValue <= 0) {
      tsTrainTrend <- tsTrainTrend + abs(minValue) + 1
    }
    tsTrainLOG <- log(tsTrainTrend)
    fcArimaLOG <- tryCatch({
      fArimaLOG <- doArima(tsTrainLOG, FALSE)
      forecast(fArimaLOG, h = horizon)$mean
    }, error = function(e){
      etsmodel <- ets(ts(tsTrainLOG))
      ts(forecast(etsmodel, h = horizon)$mean,frequency)
    }
    )
    print("exponential Trend detected!")
    fcArima <- exp(fcArimaLOG)
    if (minValue <= 0) {
      fcArima <- fcArima - abs(minValue) - 1
    }
  } else {
    fcArima <- tryCatch({
      fArima <- doArima(tsTrainTrend, FALSE)
      forecast(fArima, h = horizon)$mean
    }, error = function(e){
      etsmodel <- ets(ts(tsTrainTrend))
      ts(forecast(etsmodel, h = horizon)$mean,frequency)
    }
    )
  }
  return(fcArima)
}

#' @description Checks if the time series has a significant high remainder
#'
#' @title Checking Remainder
#' @param tvp tvp The time value pair as vector
#' @param stlRemainder The remainder part of the deccomposition of the tvp
#' @param sig.dif.factor The threshold that is to exceed for having a high remainder
#' @return If the the time series has a high remainder compared to the threshold
has.highRemainder <- function(tvp, stlRemainder, sig.dif.factor) {
  # Calculates the IQR of the remainder and original time series
  
  remainder.quantiles <- quantile(stlRemainder)
  # range from 25% to 75% quantile
  remainder.quantiles.range <-
    remainder.quantiles[4] - remainder.quantiles[2]
  tvp.quantiles <- quantile(tvp)
  tvp.quantiles.range <- tvp.quantiles[4] - tvp.quantiles[2]
  
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