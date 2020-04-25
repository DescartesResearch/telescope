#' @author Andre Bauer

set.seed(200)

#' @description Trains a recommendations for chossing the best-suited machine learning approach for a given time series.
#'
#' @title Trains the recommendation system
#' @param timeseries A list of time series
#' @param natural Optional parameter: A flag indicating wheter only natural frequencies (e.g., daily, hourly, ...) or all found frequencies shall be considered.
#' @param boxcox Optional parameter: A flag indicating if the Box-Cox transofrmation should be performed. It is not recommend to disable the transformation. TRUE by default.
#' @return An object of the class telescope.recommendation 
#' @examples
#' telescope.trainrecommender(list(AirPassengers, forecast::taylor))
#' @export
telescope.trainrecommender <- function(timeseries,natural=TRUE,boxcox=TRUE){
  if(!is.list(timeseries)){
    stop("The time series have to be stored within a list!")
  }
  
  # gets the meta-level data set
  dataset <- buildMetaLevelDataset(timeseries,natural,boxcox)
  return(trainingrecommender(dataset$chars, dataset$acc))
  
}

#' @description Trains a recommendations for chossing the best-suited machine learning approach for a given time series.
#'
#' @title Trains the recommendation system.
#' @param timeseries.chars A mxn matrix with n time series characteristics of m time series. More precisely, each row reflects the characteristics of a time series.
#' @param accuracies A mxn matrix with accuries measures of n machine learning methods based on m time series.
#' @return An object of the class telescope.recommendation. 
trainingrecommender <- function(timeseries.chars, accuracies){
  
  if(nrow(timeseries.chars)!=nrow(accuracies)){
    stop('The rows of both inputs do not match')
  }
  
  # gets method names
  methods <- colnames(accuracies)
  
  # calculates the degredation compared to the best method
  degr <- t(apply(accuracies,1, function(x) x/min(x)))
  
  est.degr <- c()
  
  # trains each a random forest method for estimation the degredation
  models <- list()
  for(m in methods){
    train <- cbind(timeseries.chars, degr[,m])
    colnames(train) <- c(colnames(timeseries.chars), "degr")
    model <- randomForest(degr ~ ., data=train, mtry = floor(sqrt(ncol(train))), sampsize = floor(0.9*nrow(train)), replace=TRUE )
    models[length(models)+1] <- list(model)
    est.degr <- cbind(est.degr,predict(model,timeseries.chars))
  }
  
  # retrieves the best method for each time series
  best <- c()
  for(j in 1:nrow(degr)){
    
    ind <- which.min(degr[j,])
    best <- c(best, methods[ind])
  }
  
  # trains a random forest to map the estimated degredations of each method to the best method
  train <- as.data.frame(est.degr)
  train$best <- factor(best)
  colnames(train) <- c(methods, 'Best')
  model <- randomForest(Best ~ ., data=train, mtry = floor(sqrt(ncol(train))), sampsize = floor(0.9*nrow(train)), replace=TRUE )
  
  return(structure(list("reg.models"=models, "clas.model"=model, "names"=methods, "chars"=colnames(timeseries.chars)), class = 'telescope.recommendation'))
  
}


#' @description Consults the recommendation system to retrieve the best method for a given time series.
#'
#' @title Retrieves the best method for a given time series.
#' @param tvp A time series. Will be ignored if tvp.stl is set.
#' @param tvp.stl Optional parameter: A stl decomposition of the time series. 
#' @param model A object of the class telescope.recommendation containing the hybrind random forest model.
#' @param natural Optional parameter: A flag indicating wheter only natural frequencies (e.g., daily, hourly, ...) or all found frequencies shall be considered.
#' @param boxcox Optional parameter: A flag indicating if the Box-Cox transofrmation should be performed. It is not recommend to disable the transformation. TRUE by default.
#' @return The best method as string.
consultrecommender <- function(tvp,tvp.stl=NULL,model,natural=TRUE,boxcox=TRUE){
  if(!is(model, 'telescope.recommendation')){
    stop('The model must be of class: telescope.recommendation!')
  }
  if(is.null(tvp.stl) && !is.ts(tvp)){
    stop('tvp has to be a time series!')
  }
  
  
  # if tvp.stl is null decompses the time series 
  if(is.null(tvp.stl)){
    if(frequency(tvp)==1){
      stop('The time series must be seasonal. That is, frequency > 1')
    }
    prep <- preprocessing(tvp,natural,boxcox)
    tvp <- prep$tvp
    tvp.stl <- stl(tvp,s.window = "periodic",t.window=length(tvp)/2)
  } else {
    if(frequency(tvp.stl$time.series[,1])==1)
      stop('The time series must be seasonal. That is, frequency > 1')
  }
  
  
  est.degr <- c()
  
  # calculates the time series characteristics for the time series
  chars <- calculateCharacteristics(tvp.stl$time.series, frequency(tvp.stl$time.series[,1]))
  
  # estimates the degregations of each method
  for(i in 1:length(model$names)){
    est.degr <- cbind(est.degr,predict(model$reg.models[[i]],chars))
  }
  
  colnames(est.degr) <- model$names
  
  # guesses the best method for the given time series
  method <- as.character(predict(model$clas.model,est.degr))
  
  
  return(method)
  
  
}
