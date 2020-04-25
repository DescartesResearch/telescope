#' @author Andre Bauer

# library(forecast)
# library(pracma)
# library(proxy)
# library(tseries)
#library(catboost)
# library(Cubist)
# library(evtree)
# library(randomForest)
# library(rpart)
# library(e1071)
# library(xgboost)

# catboostprediction <- function(train,label,test){ 
#   model <- catboost.train(catboost.load_pool(train, label = label), params = list(iterations=100 , logging_level="Silent"))
#   return(catboost.predict(model,catboost.load_pool(test)))
# }

#' @description Trains a cubist model.
#'
#' @title Trains a cubist model.
#' @param train A matrix containing the seasonal component and fourier terms for training.
#' @param label A vector containing the detrended time series.
#' @param test A matrix containing the seasonal component and fourier terms for the prediction. 
#' @return The future detrended time series. 
cubistprediction <- function(train,label,test){ 
  model <- cubist(x=train, y=label)
  return(predict(model,test))
}

#' @description Trains a evtree model.
#'
#' @title Trains a evtree model.
#' @param train A matrix containing the seasonal component and fourier terms for training.
#' @param label A vector containing the detrended time series.
#' @param test A matrix containing the seasonal component and fourier terms for the prediction. 
#' @return The future detrended time series.
evtreeprediction <- function(train,label,test){ 
  data <- as.data.frame(cbind(train,label))
  colnames(data) <- c(colnames(train), 'Target')
  model <- evtree(Target ~ ., data = data, control = evtree.control(minsplit = 2))
  return(predict(model,as.data.frame(test)))
}

#' @description Trains a nnetar model.
#'
#' @title Trains a nnetar model.
#' @param train A matrix containing the seasonal component and fourier terms for training.
#' @param label A vector containing the detrended time series.
#' @param test A matrix containing the seasonal component and fourier terms for the prediction. 
#' @return The future detrended time series.
nnetarprediction <- function(train,label,test){
  model <- nnetar(y=label,xreg=train, MaxNWts=2000)
  return(forecast(model, xreg=test,h=nrow(test))$mean)
}

#' @description Trains a random forest model.
#'
#' @title Trains a random forest model.
#' @param train A matrix containing the seasonal component and fourier terms for training.
#' @param label A vector containing the detrended time series.
#' @param test A matrix containing the seasonal component and fourier terms for the prediction. 
#' @return The future detrended time series.
randomforestprediction <- function(train,label,test){
  model <- randomForest(y=label,x=train)
  return(predict(model,test))
}

#' @description Trains a rpart model.
#'
#' @title Trains a rpart model.
#' @param train A matrix containing the seasonal component and fourier terms for training.
#' @param label A vector containing the detrended time series.
#' @param test A matrix containing the seasonal component and fourier terms for the prediction. 
#' @return The future detrended time series.
rpartprediction <- function(train,label,test){
  data <- as.data.frame(cbind(train,label))
  colnames(data) <- c(colnames(train), 'Target')
  model <- rpart(Target ~ ., data = data, method="anova",control = rpart.control(minsplit = 2, maxdepth = 30, cp = 0.000001))
  return(predict(model,as.data.frame(test)))
}

#' @description Trains a svr model.
#'
#' @title Trains a svr model.
#' @param train A matrix containing the seasonal component and fourier terms for training.
#' @param label A vector containing the detrended time series.
#' @param test A matrix containing the seasonal component and fourier terms for the prediction. 
#' @return The future detrended time series.
svrprediction <- function(train,label,test){
  model <- svm(y=label,x=train)
  return(predict(model,test))
}

#' @description Trains a xgboost model.
#'
#' @title Trains a xgboost model.
#' @param train A matrix containing the seasonal component and fourier terms for training.
#' @param label A vector containing the detrended time series.
#' @param test A matrix containing the seasonal component and fourier terms for the prediction. 
#' @return The future detrended time series.
xgboostprediction <- function(train,label,test){
  model <- doXGB.train(myts = label, cov = train, booster = "gblinear", verbose = 0)
  return(predict(model,test))
}

#' @description Calculates the error of the forecast regarding the acutal values.
#'
#' @title Calculates forecasting error.
#' @param forecast The forecast of the detrended time series.
#' @param test The acutual detrended time series. 
#' @return The mape measure.
mape <- function(forecast,test){
  
  test <- as.vector(test)
  forecast <- as.vector(forecast)
  
  nonZeros <- which(test != 0, arr.ind = T)
  forecast <- forecast[nonZeros]
  test <- test[nonZeros]
  
  et <- abs(test-forecast)
  pt <- 100*et/test
  meanMAPE <- mean(abs(pt))
  
  return(meanMAPE)
}

#' @description Fits a sinus on the data and returns the fitted errors.
#'
#' @title Fitting a sinus.
#' @param y A vector containg the time series.
#' @return The errors of the fit.
fitSin <- function(y){
  x <- seq(1:length(y))
  k <- mean(y)
  y.n <- y - k
  a <- max(abs(range(y.n)))
  z <- asin((y-k)/a)
  z.lm <- lm(x ~ z)
  f <- 1/z.lm$coefficients[2]
  phi <- z.lm$coefficients[1]
  e <- y - a*sin(f*(x-phi))+k
  return(e)
}

#' @description Quantifies how well the time series can be approximated by a sinus wave. 
#'
#' @title Checks the auto-correlation of the fitted errors.
#' @param e The error of a fit.
#' @return The auto-correlation of the errors of the fit
durbinwatson <- function(e){
  
  tryCatch({
    return(sum(diff(e)^2)/sum(e^2))
  }, error = function(e) {
    return(-1)
  })
}


#' @description Calculates the time series characteristics of the time series. 
#'
#' @title Calculates time series characteristics.
#' @param tvp.stl The STL decomposistion of the time series
#' @param freq The frequency of the time series
#' @return The time series characteritics
calculateCharacteristics <- function(tvp.stl, freq){
  
  hist.length <- length(tvp.stl[,1])
  
  remainder <- tvp.stl[,3]
  
  season <- tvp.stl[,1]
  
  detrend <-  remainder + season
  
  # Measures on original data
  xbar <- mean(detrend,na.rm=TRUE)
  s <- sd(detrend,na.rm=TRUE)
  
  # Skewness
  sk <- abs(mean((detrend-xbar)^3,na.rm=TRUE)/s^3)
  
  # Kurtosis
  k <- mean((detrend-xbar)^4,na.rm=TRUE)/s^4
  
  # Remainder skewness
  rsk <- abs(mean((remainder-xbar)^3,na.rm=TRUE)/s^3)
  
  # Remainder kurtosis
  rk <- mean((remainder-xbar)^4,na.rm=TRUE)/s^4
  
  
  # Entropy
  entropy <- c()
  
  periods <- split(detrend,ceiling(seq_along(detrend)/freq))
  for(o in 1:length(periods)){
    if(length(periods[[o]]) == freq)
      entropy <- c(entropy, approx_entropy(periods[[o]]))
  }
  
  entropy.mean <- mean(entropy)
  entropy.varcof <- sd(entropy)/mean(entropy)
  
  # Mean cosine similarity
  cos.sim <- mean(simil(array(unlist(periods), dim=c(length(periods), freq)),method="cosine"))
  
  # Durbinwatson
  dfit <- durbinwatson(fitSin(train$Season))
  
  # Season
  season <- ifelse(var(detrend,na.rm=TRUE) < 1e-10, 0,
                   max(0,min(1,1-var(remainder)/var(detrend,na.rm=TRUE))))
  
  # Serial correlation == autocorrelation
  Q <- Box.test(detrend,lag=10)$statistic/(hist.length*10)
  
  # Remainder Serial correlation
  rQ <- Box.test(remainder,lag=freq)$statistic/(hist.length*freq)
  
  # Nonlinearity
  N <- tseries::terasvirta.test(na.contiguous(ts(detrend)))$statistic
  
  # Remainder nonlinearity
  rN <- tseries::terasvirta.test(na.contiguous(ts(remainder)))$statistic
  
  # Hurst=d+0.5 where d is fractional difference.
  H <- fracdiff::fracdiff(na.contiguous(detrend),0,0)$d + 0.5
  
  PGram <- spec.pgram(x = ts(detrend),plot=FALSE)
  SpecSortPgram <- sort(PGram$spec)
  
  power2 <- 1/PGram$freq[which(PGram$spec==SpecSortPgram[length(SpecSortPgram)-1])]
  power3 <- 1/PGram$freq[which(PGram$spec==SpecSortPgram[length(SpecSortPgram)-2])]
  
  power.max <- max(SpecSortPgram)
  
  m <- cbind(hist.length, sk, k, rsk, rk, entropy.mean, entropy.varcof, cos.sim, 
             dfit, season, Q, rQ, N, rN, H,
             power2, power3, power.max)
             
             
             
             
             
  colnames(m) <- c('length', 'skewness', 'kurtosis', 'rem.skewness', 'rem.kurtosis', 'mean.entropy', 'varcoef.entropy', 'mean.cosine.similarity', 
                    'durbinwatson', 'seasonal', 'autocorrelation',  'rem.autocorrelation', 'non.linear', 'rem.non.linear', 'Hurst',  
                    'X2nd.freq', 'X3rd.freq', 'max.spec')                                                                                    
  
  return(m)
  
}


#' @description Quantifies how well each method has performed on the given time series. 
#'
#' @title Calculates the forecasting error of each method.
#' @param data The features of the time series.
#' @param label The detrended time series.
#' @return The auto-correlation of the errors of the fit
calculateAccuracies <- function(data, label){
  
  methods <- c(#"Catboost", 
    "Cubist", "Evtree", "Nnetar", "RF", "Rpart", "SVR", "XGBoost")
  
  # splits the time series into train and test data
  hist.length <- (ceiling(nrow(data) * 0.8))
  
  train <- data[1:hist.length,]
  test <- data[(hist.length+1):nrow(data),]
  
  train.label <- label[1:hist.length]
  test.label <- label[(hist.length+1):nrow(data)]
  
  # calulating the forecasting error of each method
  acc <- cbind(#mape(catboostprediction(train,train.label,test),test.label), 
               mape(cubistprediction(train,train.label,test),test.label), 
               mape(evtreeprediction(train,train.label,test),test.label), 
               mape(nnetarprediction(train,train.label,test),test.label), 
               mape(randomforestprediction(train,train.label,test),test.label), 
               mape(rpartprediction(train,train.label,test),test.label), 
               mape(svrprediction(train,train.label,test),test.label), 
               mape(xgboostprediction(train,train.label,test),test.label))
  
  colnames(acc) <- methods
  
  return(acc)
  
}

#' @description Builds the meta-level data set. That is, calculating the forecasting error and time series characteritics for each time series.
#'
#' @title Builds the meta-level data set
#' @param timeseries A list of time series
#' @param natural Optional parameter: A flag indicating wheter only natural frequencies (e.g., daily, hourly, ...) or all found frequencies shall be considered.
#' @param boxcox Optional parameter: A flag indicating if the Box-Cox transofrmation should be performed. It is not recommend to disable the transformation. TRUE by default.
#' @return The accuracies and characteristics of the time series.
buildMetaLevelDataset <- function(timeseries,natural=TRUE,boxcox=TRUE){
  if(!is.list(timeseries)){
    stop("The time series have to be stored within a list!")
  }
  
  check.seasonality <- FALSE
  
  timeseries.chars <- c()
  accuricies <- c()
  
  for(i in 1:length(timeseries)){
    tvp <- timeseries[[i]]
    if(!is.null(dim(tvp))){
      tvp <- tvp[,1]
    }
    print(paste("training progress:",i/length(timeseries)))
    # considers only seasonal time series
    if(frequency(tvp)>1){
      check.seasonality <- TRUE
      
      # stl requires at least 2 full periods
      if(length(tvp) <= 2*frequency(tvp)){
        warning(paste('The time series with index', i, 'ignored!'))
      } else {
        # adjusts the time series and retrieves dominant frequencies
        prep <- preprocessing(tvp,natural,boxcox)
        tvp <- prep$tvp
        
        tvp.mul <- msts(tvp, seasonal.periods = prep$freqs)
        
        # stl decomposition
        tvp.stl <- stl(tvp,s.window = "periodic",t.window=length(tvp)/2)
        # retrives the fourier terms based on the most dominant frequencies
        fourier.terms <- fourier(tvp.mul, K = rep(1,length(prep$freqs)))
        
        data <- as.matrix(cbind(tvp.stl$time.series[,1], fourier.terms))
        colnames(data) <- c('Season', colnames(fourier.terms))
        
        label <- as.vector(tvp - tvp.stl$time.series[,2])
        
        hist.length <- (ceiling(nrow(data) * 0.8))
        
        # calculates the time series characteristics
        timeseries.chars <- rbind(timeseries.chars,calculateCharacteristics(tvp.stl$time.series[1:hist.length,],frequency(tvp)))
        # quantifies the forecasting errof of each method
        accuricies <- rbind(accuricies, calculateAccuracies(data,label))
        
        
      }
      
      
    } else {
      warning(paste('The time series with index', i, 'ignored!'))
    }
  }
  
  if(!check.seasonality){
    stop("There are none seasonal time series!")
  }
  
  
  return(list("chars"=timeseries.chars, "acc"=accuricies))
}



