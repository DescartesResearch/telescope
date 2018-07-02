#' @author Marwin Zuefle, Andre Bauer

#' @description Estimates the frequency of the time series based on a periodogram
#'
#' @title Guess the Frequency of the Time Series
#' @param timeValuePair Timestamps with raw values (not as time series object)
#' @param asInteger Optional parameter: Boolean indicating whether frequerncy needs to be an integer
#' @param difFactor Optional parameter: A factor to determine whether a possible frequency is accepted
#' @param ithBest Optional parameter: i-th most likely frequency is tested. 1 by default.
#' @param spans Optional parameter: vector of odd integers giving the widths of modified Daniell smoothers to be used to smooth the periodogram. Null by default
#' @param debug Optional parameter: If TRUE, debugging information will be displayed. FALSE by default
#' @param PGramTvp Optional parameter: An already created periodogram. Null by default
#' @return The found frequency and the created periodogram
guessFrequencyPeriodogram <- function(timeValuePair, asInteger = TRUE, difFactor = 0.5, ithBest = 1, spans = NULL, debug = FALSE, PGramTvp = NULL) {

  num <- ithBest-1

  # Is there already a periodogram?
  if(num==0 || is.null(PGramTvp)) {
    PGramTvp <- spec.pgram(x = ts(timeValuePair),plot=FALSE,spans = spans)
  }
  # Sort the spectrum
  SpecSortPgramTvp <- sort(PGramTvp$spec)
  # Calculates the frequencies
  ProbFreq <- 1/PGramTvp$freq[which(PGramTvp$spec==SpecSortPgramTvp[length(SpecSortPgramTvp)-num])]

  if(length(ProbFreq) == 0 || ProbFreq >= length(timeValuePair)) {
    freqPeriodo <- -1
    if(debug) print("Guessed Frequency longer than time series length!")
  } else if(ProbFreq >= length(timeValuePair)/2) {
    freqPeriodo <- -1
    if(debug) print("Guessed Frequency too long: we need more than 2 periods worth of data!")
  } else {

    if(asInteger) {
      IntProbFreq <- round(ProbFreq)
      DiffTvp <- diff(timeValuePair,IntProbFreq)
    } else {
      DiffTvp <- diff(timeValuePair,ProbFreq)
    }

    # Create the periodogram for the difference of the original time series using the estimated frequency as lag
    PGramTvpDiff <- spec.pgram(ts(DiffTvp),plot=FALSE,spans = spans)

    # Find the closest frequency in the periodogram of diffs compared to the estimated frequency of the original observations
    if((min(abs(PGramTvpDiff$freq-(1/ProbFreq)))+1/ProbFreq) %in% PGramTvpDiff$freq) {
      ProbFreqDiff <- 1/(min(abs(PGramTvpDiff$freq-(1/ProbFreq)))+1/ProbFreq)
    } else if((-min(abs(PGramTvpDiff$freq-(1/ProbFreq)))+1/ProbFreq) %in% PGramTvpDiff$freq) {
      ProbFreqDiff <- 1/(-min(abs(PGramTvpDiff$freq-(1/ProbFreq)))+1/ProbFreq)
    } else {
      stop("Error determining frequency! Could not find equivalent frequency in diff.")
    }

    # Determine the spectrum of the estimated frequency for the original observations as well as the maximum spectrum
    lowerBound <- 1/(ProbFreq+0.005)
    upperBound <- 1/(ProbFreq-0.005)
    maxSpec <- PGramTvp$spec[which(PGramTvp$freq>lowerBound & PGramTvp$freq<upperBound)]
    maxSpecTest <- max(PGramTvp$spec)

    # Show debugging information
    if(debug){
      if (num==0) {
        if (maxSpec==maxSpecTest) {
          print(paste("Maxim Spectrum found:",maxSpec))
        } else {
          print("Problem with maximum spectrum!!!!!")
        }
      }
    }

    # Determine the spectrum of the estimated frequency for the diffed observations as well as the maximum spectrum
    if(0.05*ProbFreqDiff<1) {
      lowerBoundDiff <- 1/(ProbFreqDiff+1)
      upperBoundDiff <- 1/(ProbFreqDiff-1)
    } else {
      lowerBoundDiff <- 1/(ProbFreqDiff)-0.005
      upperBoundDiff <- 1/(ProbFreqDiff)+0.005
    }
    maxSpecDiff <- PGramTvpDiff$spec[which(PGramTvpDiff$freq>lowerBoundDiff & PGramTvpDiff$freq<upperBoundDiff)]
    maxSpecDiff <- max(maxSpecDiff)

    # Compare the spectrum of the original periodogram to the one of the diffed observations
    if(maxSpecDiff < difFactor*maxSpec) {
      if(asInteger) {
        freqPeriodo <- IntProbFreq
      } else {
        freqPeriodo <- ProbFreq
      }
    } else {
      if(debug) print("Some problems occured! The max spectrum is not a real season!")
      freqPeriodo <- -1
    }
  }

  return(list("frequency" = freqPeriodo, "pgram" = PGramTvp))
}

#' @description Calculates the frequency for the time series
#'
#' @title Calcuates the Frequency of the Time Series
#' @param timeValuePair Timestamps with raw values (not as time series object)
#' @param asInteger Optional parameter: Boolean indicating whether frequerncy needs to be an integer. TRUE by default
#' @param difFactor Optional parameter: A factor to determine whether a possible frequency is accepted. 0.5 by default
#' @param tolerance Optional parameter: A tolerance factor for matching estimated frequencies to listed frequencies. 0.05 by default
#' @param maxIters Optional parameter: Number of maximal rounds to estimate the frequency. 10 by default.
#' @param ithBest Optional parameter: i-th most likely frequency is tested. 1 by default.
#' @param PGramTvp Optional parameter: An already created periodogram. Null by default
#' @return The frequency that appears in the periodogram and the list of most frequent frequencys, the created perodogram and the last iteration
calcFrequencyPeriodogram <- function(timeValuePair, asInteger = TRUE, difFactor = 0.5, tolerance = 0.05, maxIters = 10, ithBest = 1, PGramTvp = NULL, debug = FALSE) {

  # typical measurements per hour
  mPerH <- c(1,2,4,6,12,60)

  # typical first seasons
  sea <- c(1,24)

  # typical second seasons
  sea2 <- c(1,7,30,365)

  # Create vector with all combinations
  allComb <- c()
  for(i in 1:length(mPerH)) {
    for(j in 1:length(sea)) {
      for(k in 1:length(sea2)) {
        allComb <- c(allComb,mPerH[i]*sea[j]*sea2[k])
      }
    }
  }

  allComb <- sort(unique(allComb))
  allComb <- allComb[which(allComb>=4)]


  frequency <- -1
  numIters <- 1
  while(frequency==-1) {
    # If there is no suitable frequency found during the maximum amount of iterations, the most dominant frequency, the one
    # determined using ithBest = 1, is returned
    if(numIters>maxIters) {
      if(debug) print("Periodogram could not find dominant frequency")
      freq <- guessFrequencyPeriodogram(timeValuePair,asInteger,difFactor,ithBest = 1,PGramTvp = PGramTvp)
      PGramTvp <- freq$pgram
      break()
    }

    # Add a span if there is no suitable frequency found during the first two iterations
    if(numIters<=2) {
      freq <- guessFrequencyPeriodogram(timeValuePair,asInteger,difFactor,ithBest = ithBest,PGramTvp = PGramTvp)
    } else {
      freq <- guessFrequencyPeriodogram(timeValuePair,asInteger,difFactor,ithBest = ithBest,spans = 5,PGramTvp = PGramTvp)
    }

    if(debug) print(paste("Iteration:",numIters, "testing frequency:",freq$frequency))
    PGramTvp <- freq$pgram
    a <- freq$frequency

    # Create a tolerance area for matching estimated frequencies to listed frequencies
    a.p <- c((1-tolerance)*a,(1+tolerance)*a)
    if(tolerance*a<1) {
      a.p <- c(a-1,a+1)
    }

    # Find all suitable frequencies which are between the bounds of a (a.p)
    a.p <- round(a.p)
    possible.freqs <- allComb[which(findInterval(allComb,a.p)==1)]
    # Calculate the differences of the frequencies found in the interval and the estimated frequency
    deltas <- abs(possible.freqs-a)

    # If there are any frequencies found in the interval, select the one with the smallest difference
    # Otherwise, there is no suitable frequency found in this iteration
    if(length(deltas)>0) {
      frequency <- possible.freqs[deltas==min(deltas)]
    } else {
      frequency <- c()
    }
    if(length(frequency)==0) {
      frequency <- -1
    }
    lastIterFreq1 <- numIters
    numIters <- numIters + 1
    ithBest <- ithBest + 1
  }

  # If there is no "good" frequency found, set frequency to 2 as STL requries at least this value
  if(frequency == -1) {
    frequency = 2
    if(debug) print("No frequency found. Set frequency to: 2")
  } else {
    if(debug) print(paste("Accepting frequency:",frequency))
  }

  return(list("frequency" = frequency, "pgram" = freq$pgram, "lastIterfreq" = lastIterFreq1))
}
