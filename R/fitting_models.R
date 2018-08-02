#' @author Marwin Zuefle, Andre Bauer

#' @description Checks the model of the trend: either linear or exponential
#'
#' @title Fitting the Model of the Trend
#' @param tvp The time-value pair
#' @param frequency The frequency of the time-value pair
#' @param difFactor Optional parameter: The factor how much exp needs to be better than linear. 1.5 by default
#' @param debug Optional parameter: If TRUE, debugging information will be displayed. FALSE by default
#' @return The model, either linear or exp and if the estimation is riskys
fittingModels <- function(stl, frequency, difFactor = 1.5, debug = FALSE) {

  stlTrend <- stl$time.series[,2]
  # Moves trend to greater than 1 due to log
  if(min(stlTrend) < 1) {
    stlTrend <- stlTrend + (1 - min(stlTrend))
  }

  times <- seq(frequency,(length(stlTrend)+frequency-1))/frequency

  # Fiting a linear model
  lin <- lm(stlTrend~times)
  len <- length(stlTrend)
  # y axis intercept
  cut <- coef(lin)[1]
  # slope
  m <- coef(lin)[2]

  # Shows some debugging information
  if(debug) {
    plot(lin)
    lines(c(1,(len+frequency-1/frequency)),c(m*frequency+cut,(len+frequency-1)*m+cut),col="red")
  }

  # Fitting an exponential model
  stlTrendLOG <- log(stlTrend)
  lmEXP <- lm(stlTrendLOG~times)
  # y axis intercept
  cutEXP <- coef(lmEXP)[1]
  # slope
  mEXP <- coef(lmEXP)[2]
  if(debug) {
    plot(stlTrendLOG)
    lines(c(1,(len+frequency-1)/frequency),c(mEXP*frequency+cutEXP,(len+frequency-1)*mEXP+cutEXP),col="red")
  }

  # Creates both models
  Counts.exponential2 <- exp(predict(lmEXP,list(Time=stl$time.series[,0])))
  Counts.linear <- predict(lin,list(Time=stl$time.series[,0]))

  # plot all types of models
  if(debug){
    plot(exp(stlTrendLOG))
    expModel <- matrix(nrow=len,ncol = 2)
    for(i in frequency:(frequency+len-1)) {
      expModel[(i-frequency+1),1] <- i/frequency
      expModel[(i-frequency+1),2] <- exp(i*(mEXP)+cutEXP)
    }
    lines(expModel[,1],expModel[,2],col="blue")
    lines(times, Counts.exponential2,lwd=2, col = "green", xlab = "Time (s)", ylab = "Counts")
    lines(times, Counts.linear,lwd=2, col = "purple", xlab = "Time (s)", ylab = "Counts")
  }

  # Calculates RSME between the models and the trend
  RMSElog <- sqrt(mean((stlTrend-Counts.exponential2)^2))
  RMSE <- sqrt(mean((stlTrend-Counts.linear)^2))

  if(debug){
    print(paste("RMSE exp fit:",RMSElog))
    print(paste("RMSE lin fit:",RMSE))
  }

  if (RMSE <= difFactor*RMSElog) {
    # if the difference is small, the trend model estimation is risky
    if(RMSE <= 0.8*difFactor*RMSElog) {
      risky_trend_model <- FALSE
    } else {
      risky_trend_model <- TRUE
    }
    ret <- "linear"
  } else {
    # if the difference is small, the trend model estimation is risky
    if(RMSE > 1.2*difFactor*RMSElog) {
      risky_trend_model <- FALSE
    } else {
      risky_trend_model <- TRUE
    }
    ret <- "exp"
  }
  return(
    list(
      "trendmodel" = ret,
      "risky_trend_model" = risky_trend_model
    )
  )
}

