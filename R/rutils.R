# codeutils ------------------------------------


#' @title classDF - tabluate the class of each column in a da
#'
#' @description classDF - tabluate the class of each column in a dataframe.
#'
#' @param dataframe - the input dataframe for examination
#' @return generates paired column names with their classes
#' @export classDF
#' @examples
#' \dontrun{
#'  data(ChickWeight)
#'  classDF(ChickWeight)
#' }
classDF <- function(dataframe) {
  nvar <- dim(dataframe)[2]
  for (i in 1:nvar) cat(colnames(dataframe)[i],class(dataframe[,i]),"\n")
} # end of class_DF

#' @title countones used in apply to count the number of ones in a vector
#'
#' @description countones used in apply to count number of ones in a vector
#' @param invect vector of values
#' @return A single value of zero or the number of ones
#' @export countones
#' @examples
#' \dontrun{
#'   set.seed(12346)
#'   x <- trunc(runif(10)*10)
#'   x
#'   countones(x)  # should be 2
#' }
countones <- function(invect) {
  pick <- which(invect == 1)
  return(length(pick))
}

#' @title countzeros used in apply to count the number of zeros in a vector
#'
#' @description countzeros used in apply to count zeros in a vector
#' @param invect vector of values
#' @return A single value of zero or the number of zeros
#' @export countzeros
#' @examples
#' \dontrun{
#'   set.seed(12346)
#'   x <- trunc(runif(10)*10)
#'   x
#'   countzeros(x)  # should be 1
#' }
countzeros <- function(invect) {
  pick <- which(invect == 0.0)
  return(length(pick))
}

#' @title countgtzero used in apply to count the number >0 in a vector
#'
#' @description countgtzero used in apply to count number >0 in a vector
#' @param invect vector of values
#' @return A single value of number of values > 0
#' @export countgtzero
#' @examples
#' \dontrun{
#'   set.seed(12346)
#'   x <- trunc(runif(10)*10)
#'   x
#'   countgtzero(x)  # should be 9
#' }
countgtzero <- function(invect) {
  pick <- which(invect > 0)
  return(length(pick))
}

#' @title countNAs used in apply to count the number of NAs in a vector
#'
#' @description countNAs used in apply to count number of NAs in a vector
#' @param invect vector of values
#' @return A single value of zero or the number of NAs
#' @export countNAs
#' @examples
#' \dontrun{
#'   set.seed(12346)
#'   x <- trunc(runif(10)*10)
#'   x[c(3,7)] <- NA
#'   countNAs(x)  # should be 2
#' }
countNAs <- function(invect) {
  pick <- which(is.na(invect))
  return(length(pick))
}

#' @title countgtOne used in apply to count the number > 1 in a vector
#'
#' @description countgtOne used in apply to count the number > 1 in a vector
#' @param invect vector of values
#' @return A single value of zero or the number of NAs
#' @export countgtOne
#' @examples
#' \dontrun{
#'   set.seed(12346)
#'   x <- trunc(runif(10)*10)
#'   x
#'   countgtone(x)  # should be 7
#' }
countgtOne <- function(invect) {
  pick1 <- which(invect > 1.0)
  return(length(pick1)/length(invect))
}

#' @title facttonum converts a vector of numeric factors into numbers
#'
#' @description facttonum converts a vector of numeric factors into numbers.
#'     If the factors are not numeric then the outcome will be a series of 
#'     NA. It is up to you to apply this function only to numeric factors. 
#'     A warning will be thrown if the resulting output vector contains NAs
#'
#' @param invect vector of numeric factors to be converted back to numbers
#'
#' @return an output vector of numbers instead of the input factors
#' @export
#'
#' @examples
#' \dontrun{
#'  DepCat <- as.factor(rep(seq(100,600,100),2)); DepCat
#'  5 * DepCat[3]
#'  as.numeric(levels(DepCat))  # #only converts levels not the replicates
#'  DepCat <- facttonum(DepCat)
#'  5 * DepCat[3]
#'  x <- factor(letters)
#'  facttonum(x)
#' }
facttonum <- function(invect){
  if (class(invect) == "factor") {
    outvect <- suppressWarnings(as.numeric(levels(invect))[invect])
  }
  if (class(invect) == "numeric") outvect <- invect
  if (any(is.na(outvect)))
    warning("NAs produced, input vector may have non-numbers present \n")
  return(outvect)
} # end of facttonum

#' @title freqMean calculates the mean and stdev of count data
#'
#' @description freqMean calculates the mean and stdev of count data
#'     it requires both the values and their associated counts and
#'     return a vector of two numbers.
#'
#' @param values the values for which there are counts
#' @param infreqs the counts for each of the values empty cells can be
#'     either 0 or NA
#'
#' @return a vector containing the mean and st.dev.
#' @export
#'
#' @examples
#' \dontrun{
#' vals <- c(1,2,3,4,5)
#' counts <- c(3,NA,7,4,2)
#' freqMean(vals,counts)  # should give 3.125 and 1.258306
#' }
freqMean <- function(values,infreqs) {
  N <- length(values)
  if (N != length(infreqs)) {
    cat("vectors have different lengths \n")
    ans <- c(NA,NA)
    names(ans) <- c("mean","stdev")
  } else {
    nobs <- sum(infreqs,na.rm=T)
    sumX <- sum(values * infreqs,na.rm=T)
    av <- sumX/nobs
    if (length(infreqs[infreqs > 0.01]) > 1) {
      sumX2 <- sum(values * values * infreqs,na.rm=T)
      stdev <- sqrt((sumX2 - (sumX * sumX)/nobs)/(nobs-1))
    } else { stdev <- NA
    }
    ans <- c(av,stdev)
    names(ans) <- c("mean","stdev")
  }
  return(ans)
} # end of freqMean

#' @title geomean log-normal bias corrected geometric mean of a vector
#'
#' @description Calculates log-normal bias corrected geometric mean of a 
#'     vector. NAs and zeros are removed from consideration.
#' @param invect is a vector of numbers in linear space.
#' @return The bias-corrected geometric mean of the vector
#' @export geomean
#' @examples
#' \dontrun{
#'  x <- c(1,2,3,4,5,6,7,8,9)
#'  geomean(x)
#' }
geomean <- function(invect) {
  pick <- which((invect <= 0.0))
  if (length(pick) == 0) {
    avCE <- mean(log(invect),na.rm=TRUE)
    stdev <- sd(log(invect),na.rm=TRUE)
  } else {
    avCE <- mean(log(invect[-pick]),na.rm=TRUE)
    stdev <- sd(log(invect[-pick]),na.rm=TRUE)
  }
  gmean <- exp(avCE + (stdev^2)/2)
  return(gmean)
}  # end of geomean

#' @title getmin generates the lower bound for a plot
#'
#' @description getmin generates lower bound for a plot where it is unknown
#'     whether the minimum is less than zero of not. If less than 0 then
#'     multiplying by the default mult of 1.05 works well but if the outcome
#'     if > 0 then the multiplier needs to be adjusted appropriately so 
#'     the minimum is slightly lower than the minimum of the data
#'
#' @param x the vector of data to be tested for its minimum
#' @param mult the multiplier for both ends, defaults to 1.05 (=0.95 if >0)
#'
#' @return a suitable lower bound for a plot if required
#' @export
#'
#' @examples
#' \dontrun{
#' vect <- rnorm(10,mean=0,sd=2)
#' sort(vect)
#' getmin(vect,mult=1.0)
#' }
getmin <- function(x,mult=1.05) {
  ymin <- min(x,na.rm=TRUE)
  if (ymin < 0) {
    ymin <- ymin * mult
  } else {
    ymin <- ymin * (2 - mult)
  }
  return(ymin)
} # end of getmin

#' @title getmax generates the upper bound for a plot
#'
#' @description getmax generates upper bound for a plot where it is unknown
#'     whether the maximum is greater than zero of not. If > 0 then
#'     multiplying by the default mult of 1.05 works well but if the outcome
#'     if < 0 then the multiplier needs to be adjusted appropriately so the 
#'     maximum is slightly higher than the maximum of the data
#'
#' @param x the vector of data to be tested for its maximum
#' @param mult the multiplier for both ends, defaults to 1.05 (=0.95 if < 0)
#'
#' @return a suitable upper bound for a plot if required
#' @export
#'
#' @examples
#' \dontrun{
#'  vect <- rnorm(10,mean=0,sd=2)
#'  sort(vect,decreasing=TRUE)
#'  getmax(vect,mult=1.0)
#'  vect <- rnorm(10,mean = -5,sd = 1.5)
#'  sort(vect,decreasing=TRUE)
#'  getmax(vect,mult=1.0)
#' }
getmax <- function(x,mult=1.05) {
  ymax <- max(x,na.rm=TRUE)
  if (ymax > 0) {
    ymax <- ymax * mult
  } else {
    ymax <- ymax * (2 - mult)
  }
  return(ymax)
} # end of getmax

#' @title getseed generates a random number seed
#' 
#' @description getseed generates a seed for use within set.seed. 
#'     It produces up to a 6 digit integer from the Sys.time. This
#'     Initially, at the start of a session there is no seed; a new one 
#'     is created from the current time and the process ID when one is 
#'     first required. Here, in getseed, we do not use the process ID so 
#'     the process is not identical but this at least allows the 
#'     set.seed value to be stored should the need to repeat a set of 
#'     simulations arise. The process generates up to a six digit number
#'     it then randomly reorders those digits and that becomes the seed.
#'     That way, if you were to call getseed in quick succession the
#'     seeds generated should differ even when they are generated close
#'     together in time.
#'
#' @return  an integer up to 7 digits long 
#' @export
#'
#' @examples
#' useseed <- getseed()
#' set.seed(useseed)
#' rnorm(5)
#' set.seed(12345)
#' rnorm(5)
#' set.seed(useseed)
#' rnorm(5)
getseed <- function() {
  pickseed <- as.character(as.integer(Sys.time()))
  nc <- nchar(pickseed)
  if (nc > 7) pickseed <- substr(pickseed,(nc-6),nc)
  nc <- nchar(pickseed)  
  pseed <- unlist(strsplit(pickseed,split=character(0)))
  pseed <- sample(pseed,nc)
  newseed <- paste(pseed,collapse="")
  newseed <- as.numeric(newseed)
  return(newseed)
} # end of getseed

#' @title gettime calculates time in seconds passed each day
#' 
#' @description gettime is a function designed to facilitate the measurement
#'     of time between intervals within R software that are expected to
#'     take a maximum of hours. It calculates the time as seconds elapsed 
#'     from the start of each day. As long as the timing of events does not
#'     pass from one day to the next accurate results will be generated.
#'
#' @return the time in seconds from the start of a day
#' @export
#'
#' @examples
#' \dontrun{
#'   begin <- gettime()
#'   for (i in 1:1e6) sqrt(i)
#'   finish <- gettime()
#'   print(finish - begin)
#' }
gettime <- function() {
  tim <- unlist(as.POSIXlt(Sys.time()))
  hr <- as.numeric(tim["hour"])*3600
  min <- as.numeric(tim["min"])*60
  sec <- as.numeric(tim["sec"])
  return(hr+min+sec)
} # end of gettime

#' @title greplow - uses tolower in the search for the pattern
#'
#' @description greplow a grep implementation that ignores the case of 
#'     either the search pattern or the object to be search. Both are 
#'     converted to lower case before using grep.
#' @param pattern - the text to search for in x
#' @param x - the vector or object within which to search for 'pattern' once
#'    both have been converted to lowercase.
#'
#' @return the index location within x of 'pattern', if it is present, 
#'     an empty integer if not
#' @export greplow
#'
#' @examples
#' \dontrun{
#' txt <- c("Long","Lat","LongE","LatE","Depth","Zone","Effort","Method")
#' greplow("zone",txt)
#' greplow("Zone",txt)
#' greplow("long",txt)
#' }
greplow <- function(pattern,x) {
  return(grep(tolower(pattern),tolower(x)))
}

#' @title info gets the dimension or length of a matrix, array, data.frame or list
#' 
#' @description info gets the dimension or length of a matrix, array, 
#'    data.frame or list. It is safer than dim because if the object is a
#'    list dim fails.
#'
#' @param invar an object that is either a matrix, an array, a data.frame or
#'     a list
#' @param verbose should the head of the object be printed to console. 
#'     default=FALSE
#'
#' @return the dimensions of the object
#' @export
#'
#' @examples
#' x <- array(rnorm(125,mean=5,sd=1),dim=c(5,5,5))
#' info(x,FALSE)
#' x <- list(x=5,y=6,z=7)
#' info(x)
info <- function(invar,verbose=FALSE) {
  cat("Class: ",class(invar),"\n")
  str(invar,max.level=1)
  cat("\n")
  categories <-  c("matrix","array","data.frame")
  if (class(invar) %in% categories) {
    cat("Dimension: ",dim(invar),"\n")
    if (verbose) print(head(invar,2))
  } else {
    cat("Length: ",length(invar),"\n")
  }
} # end of info

#' @title makelabel generates a label from text and values
#'
#' @description makelabel It is common to want a label with text and a series 
#'     of values. But paste and paste0 cycles the text and the values. To
#'     avoid this makelabel first combines the values as text and then
#'     adds the input text to the front of the values
#'
#' @param txt the input text for the label, can be empty
#' @param vect the series of values to be included in the label
#' @param sep the separator for the components; defaults to  '_'
#' @param digits how many significant digits for the values; default = 3
#'
#' @return a character string made up of text and values
#' @export
#'
#' @examples
#' pars <- c(18.3319532,33.7935124,3.0378107,6.0194465,0.5815360,0.4270468)
#' makelabel("Cohort1",pars[c(1,3,5)],sep="__")
#' makelabel("",pars[c(1,3,5)],sep="__",digits=4)
makelabel <- function (txt, vect, sep = "_", digits = 3) {
  label <- round(vect[1], digits)
  if (length(vect) > 1) {
    nnum <- length(vect)
    for (i in 2:nnum) label <- paste(label, round(vect[i], digits), sep = sep)
  }
  if (nchar(txt) > 0) label <- paste0(txt,sep,tmp)
  return(label)
} # end of makelabel

#' @title makeUnit generates a unit matrix whose diagonal can be changed
#' 
#' @description makeUnit generates a unit matrix but includes the facility
#'     to alter the diagonal value away from 1.0 if desired.
#'
#' @param N the order of the matrix
#' @param diagvalue defaults to 1.0, but otherwise can be a different 
#'     constant or a vector of dimension N
#'
#' @return a square matrix defaulting to a unit matrix
#' @export
#'
#' @examples
#' \dontrun{
#'   makeUnit(4)
#'   surv <- exp(-0.2)
#'   makeUnit(4,surv)
#' }
makeUnit <- function(N,diagvalue=1.0) {
  N <-trunc(N)
  UnitM <- matrix(0,nrow=N,ncol=N,dimnames=list(1:N,1:N))
  diag(UnitM) <- diagvalue
  return(UnitM)
}  # end of makeUnit


#' @title magnitude returns the magnitude of numbers
#'
#' @description magnitude is useful when using an
#'     optimizer such as optim, which uses a parscale parameter.
#'     magnitude can determine the respective parscale value for each
#'     parameter value.
#'
#' @param x the vector of numbers (parameters) whose magnitudes are
#'     needed
#'
#' @return a vector of magnitudes
#' @export
#'
#' @examples
#' \dontrun{
#'   x <- c(0,0.03,0.3,3,30,300,3000)
#'   magnitude(x)
#' }
magnitude <- function(x) {
  return(10^(floor(log10(abs(x)))))
}


#' @title outfit tidy print of output from optim, nlminb, or nlm
#'
#' @description outfit takes in the output list from either optim,
#'     nlminb, or nlm and prints it more tidily to the console, In the
#'     case of nlm it also prints the conclusion regarding the
#'     solution. It might be more effective to implement an S3 method.
#'
#' @param inopt the list object output by nlm, nlminb, or optim
#' @param backtran a logical default = TRUE If TRUE it assumes
#'     that the parameters have been log-transformed for stability
#'     and need back-transforming
#' @param digits the number of digits to round the backtransformed 
#'     parameters. defaults to 5.
#' @param title character string used to label the output if desired,
#'     default = empty character string
#' @param parnames default="" which means the estimated parameters
#'     will merely be numbered. If a vector of names is given 
#'     then this will be used instead, at least, for nlm and optim.
#'
#' @return nothing but it does print the list to the console tidily
#' @export
#'
#' @examples
#'  x <- 1:10  # generate power function data from c(2,2) + random
#'  y <- c(2.07,8.2,19.28,40.4,37.8,64.68,100.2,129.11,151.77,218.94)
#'  alldat <- cbind(x=x,y=y)
#'  pow <- function(pars,x) return(pars[1] * x ^ pars[2])
#'  ssq <- function(pars,indat) {
#'     return(sum((indat[,"y"] - pow(pars,indat[,"x"]))^2))
#'  }  # fit a power curve using normal random errors
#'  pars <- c(2,2)
#'  best <- nlm(f=ssq,p=pars,typsize=magnitude(pars),indat=alldat)
#'  outfit(best,backtran=FALSE) #a=1.3134, b=2.2029 ssq=571.5804
outfit <- function(inopt,backtran=TRUE,digits=5,title="",
                   parnames=""){
  #  inopt=bestvB; backtran = FALSE; digits=5; title=""; parnames=""
  nlmcode <- c("gradient close to 0, probably solution",
               ">1 iterates in tolerance, probably solution",
               "Either ~local min or steptol too small",
               "iteration limit exceeded",
               "stepmax exceeded ,5 times")
  if (length(grep("value",names(inopt))) > 0) { # optim
    cat("optim solution: ", title,"\n")
    cat("minimum     : ",inopt$value,"\n")
    cat("iterations  : ",inopt$counts," iterations, gradient\n")
    cat("code        : ",inopt$convergence,"\n")
    if (backtran) {
      ans <- cbind(par=inopt$par,transpar=round(exp(inopt$par),digits))
    } else {
      ans <- t(inopt$par)
    }
    if ((length(parnames) > 1) & (length(parnames) == length(inopt$par))) {
      rownames(ans) <- parnames
    } else {
      rownames(ans) <- 1:length(inopt$par)
    }
    print(ans)
    cat("message     : ",inopt$message,"\n")
  } # end of optim
  if (length(grep("minimum",names(inopt))) > 0) {  # nlm - preferred
    cat("nlm solution: ", title,"\n")
    cat("minimum     : ",inopt$minimum,"\n")
    cat("iterations  : ",inopt$iterations,"\n")
    cat("code        : ",inopt$code,nlmcode[inopt$code],"\n")
    if (backtran) {
      ans <- cbind(par=inopt$estimate,gradient=inopt$gradient,
                   transpar=round(exp(inopt$estimate),digits))
    } else {
      ans <- cbind(par=inopt$estimate,gradient=inopt$gradient)
    }
    if ((length(parnames) > 1) & 
        (length(parnames) == length(inopt$estimate))) {
      rownames(ans) <- parnames
    } else {
      rownames(ans) <- 1:length(inopt$estimate)
    }
    print(ans)
  } # end of nlm
  if (length(grep("objective",names(inopt))) > 0) {
    cat("nlminb solution: ", title,"\n")   # nlminb seems to be deprecated
    cat("par        : ",inopt$par,"\n")
    cat("minimum    : ",inopt$objective,"\n")
    cat("iterations : ",inopt$iterations,"\n")
    cat("code       : ",inopt$evaluations," iterations, gradient\n")
    cat("message    : ",inopt$message,"\n")
  }
  if (length(grep("hessian",names(inopt))) > 0) {
    cat("hessian     : \n")
    print(inopt$hessian)
  }
} # end of outfit

#' @title printV returns a vector cbinded to 1:length(invect)
#'
#' @description printV takes an input vector and generates another vector of
#'     numbers 1:length(invect) which it cbinds to itself. This is primarily
#'     useful when trying to print out a vector which can be clumsy to read 
#'     when print across the screen. applying printV leads to a single 
#'     vector being printed down the screen
#'
#' @param invect the input vector to be more easily visualized, this can be
#'     numbers, characters, or logical. If logical the TRUE and FALSE are
#'     converted to 1's and 0's
#' @param label the column labels for vector, default is index and value
#'
#' @return a dataframe containing the vector 1:length(invect), and invect.
#' @export
#'
#' @examples
#' \dontrun{
#' vec <- rnorm(10,mean=20,sd=2)
#' printV(vec)
#' vec <- letters
#' printV(vec)
#' vec <- c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE)
#' printV(vec,label=c("index","logicstate"))
#' }
printV <- function(invect,label=c("value","index")) {
  n <- length(invect)
  outvect <- as.data.frame(cbind(invect,1:n))
  colnames(outvect) <- label
  return(outvect)
} # end of print_V

#' @title properties - used to check a data.frame before standardization
#'
#' @description properties - used to check a data.frame before
#'     standardization. The maximum and minimum are constrained to four
#'     decimal places. It allows for columns of NAs and for Posix 
#'     columns.
#' @param indat the data.frame containing the data fields to be used
#'     in the subsequent standardization. It tabulates the number of
#'     NAs and the number of unique values for each variable and finds
#'     the minimum and maximum of the numeric variables
#' @param dimout determines whether or noth the dimensions of the data.frame
#'     are printed to the screen or not; defaults to FALSE
#' @return a data.frame with the rows being each variable from the input
#'     input data.frame and the columns being the number of NAs, the
#'     number of unique values, and minimum and maximum (where possible).
#' @export properties
#' @examples
#' \dontrun{
#'  data(abdat)
#'  properties(abdat$fish)
#' }
properties <- function(indat,dimout=FALSE) {  # indat=sps1; dimout=FALSE
  dominmax <- function(x) {
    if (length(which(x > 0)) == 0) return(c(NA,NA))
    mini <- min(x,na.rm=TRUE)
    maxi <- max(x,na.rm=TRUE)
    return(c(mini,maxi))
  }
  if(dimout) print(dim(indat))
  isna <- sapply(indat,function(x) sum(is.na(x)))
  uniques <- sapply(indat,function(x) length(unique(x)))
  columns <- length(indat)
  clas <- character(columns)
  for (i in 1:columns) {
    clas[i] <- class(indat[,i])[1]
  }
  numbers <- c("integer","numeric")
  pick <- which(clas %in% numbers)
  minimum <- numeric(length(uniques))
  maximum <- minimum
  for (i in 1:length(pick)) {
    minmax <- dominmax(indat[,pick[i]])
    minimum[pick[i]] <- minmax[1]
    maximum[pick[i]] <- minmax[2]
  }
  pick <- which((clas == "character") & (isna == 0))
  if (length(pick) > 0) {
    for (i in 1:length(pick)) {
      pickna <- which(indat[,pick[i]] == "NA")
      if (length(pickna) > 0) isna[pick[i]] <- length(pickna)
    }
  }  
  index <- 1:length(isna)
  props <- as.data.frame(cbind(index,isna,uniques,clas,round(minimum,4),
                               round(maximum,4),t(indat[1,])))
  colnames(props) <- c("Index","isNA","Unique","Class","Min",
                       "Max","Example")
  return(props)
} # end of properties


#' @title quants used in apply to estimate quantiles across a vector
#'
#' @description quants used in 'apply' to estimate quantiles across a vector
#' @param invect vector of values
#' @param probs the quantiles wanted in the outputs; default = 
#'     c(0.025,0.05,0.5,0.95,0.975)
#' @return a vector of the c(0.025,0.05,0.5,0.95,0.975) quantiles or
#'     whatever is input to probs
#' @export quants
#' @examples
#' \dontrun{
#'  x <- runif(1000)
#'  quants(x)
#'  quants(x,probs=c(0.075,0.5,0.925))
#' }
quants <- function(invect,probs = c(0.025,0.05,0.5,0.95,0.975)) {
  ans <- quantile(invect,probs =probs,na.rm=T)
  return(ans)
}

#' @title removeEmpty removes empty strings from a vector of strings
#'
#' @description removeEmpty removes empty strings from a vector of strings.
#'     Such spaces often created by spurious commas at the end of lines. It
#'     also removes strings made up only of spaces and removes spaces from
#'     inside of inidivdual chunks of text.
#'
#' @param invect vector of input strings, possibly containing empty strings
#'
#' @return a possibly NULL vector of strings
#' @export
#'
#' @examples
#' \dontrun{
#' x <- c("1","","2","","   ","3"," ","4","","a string","end")
#' x
#' length(x)
#' length(removeEmpty(x))
#' removeEmpty(x)
#' }
removeEmpty <- function(invect) {
  tmp <- gsub(" ","",invect)
  tmp <- tmp[nchar(tmp) > 0]
  return(tmp)
}

#' @title revsum generates a vector of the cumulative sum from n to 1 
#'
#' @description revsum generates a vector of the cumulative sum of an input
#'     vector from n to 1 rather than from 1 - n, as in cumsum.
#' 
#' @param x an input vector
#'
#' @return a vector of cumulative values from n to 2
#' @export
#'
#' @examples
#' x <- c(1,2,3,4,5)/15
#' print(round(cbind(x,cumsum(x),revsum(x)),3))
revsum <- function(x) {
  n <- length(x)
  if (n < 2) warning("Input vector in revsum less than length 2.  \n")
  ans <- numeric(n)
  ans[n] <- x[n]
  for (i in (n-1):1) ans[i] <- ans[i+1] + x[i] 
  return(ans)
} # end of revsum


#' @title str1 a simple replacement for str(x,max.level=1)
#' 
#' @description str1 an abbreviated replacement for str(x,max.level=1), which I 
#'     put together because to often I make a typo when typing out the full
#'     str syntax. Hence I find str1 helpful
#'
#' @param x the object whose structure is to be listed
#'
#' @return str(x,max.level=1)
#' @export
#'
#' @examples
#' x <- matrix(rnorm(25,mean=5,sd=1),nrow=5,ncol=5)
#' str1(x)
str1 <- function(x){
  return(str(x,max.level=1))
}

#' @title str2 a simple replacement for str(x,max.level=2)
#' 
#' @description str2 an abbreviated replacement for str(x,max.level=2), which I 
#'     put together because to often I make a typo when typing out the full
#'     str syntax. For when str1 is not detailed enough.
#'
#' @param x the object whose structure is to be listed
#'
#' @return str(x,max.level=2)
#' @export
#'
#' @examples
#' x <- matrix(rnorm(25,mean=5,sd=1),nrow=5,ncol=5)
#' str2(x)
str2 <- function(x){
  return(str(x,max.level=2))
}

#' @title tidynames can replace awkward data.frame names with better ones
#'
#' @description tidynames can replace awkward or overly long data.frame
#'     column names with better ones that are easier to use. It also
#'     permits one to maintain the same set of column names within an
#'     analysis even when the source data.frame includes alterations.
#'
#' @param columns the vector of names that should include the ones to be
#'     altered
#' @param replace the names to be changed, as a vector of character
#'     strings
#' @param repwith the replacement names as a vector of character strings
#'
#' @return a vector of new columns names
#' @export
#'
#' @examples
#'  print("wait")
tidynames <- function(columns,replace,repwith) {
  nreplace <- length(replace)
  if (nreplace != length(repwith))
    stop("Different number of names in replace and repwith \n")
  for (i in 1:nreplace) {
    pick <- grep(replace[i],columns)
    #cat(i,pick,"\n")
    if (pick[1] > 0) {
      columns[pick[1]] <- repwith[i]
    } else {
      warning(paste0(replace[i]," not in the dataset"))
    }
  }
  return(columns)
} # end of tidynames

#' @title toXL copies a data.frame or matrix to the clipboard
#'
#' @description toXL copies a data.frame or matrix to the clipboard
#'    so one can then switch to Excel and just type ctrl + V to paste the
#'    data in place
#'
#' @param x a vector or matrix
#' @param output - a boolean determining whether to print the object to the
#'    screen as well as the clipboard; defaults to FALSE
#' @return Places the object 'x' into the clipboard ready for pasting
#' @export toXL
#' @examples
#' datamatrix <- matrix(data=rnorm(100),nrow=10,ncol=10)
#' colnames(datamatrix) <- paste0("A",1:10)
#' rownames(datamatrix) <- paste0("B",1:10)
#' toXL(datamatrix,output=TRUE)
toXL <- function(x,output=FALSE) {
  write.table(x,"clipboard",sep="\t")
  if(output) print(x)
}

#' @title which.closest find the number closest to a given value
#'
#' @description which.closest finds either the number in a vector which is
#'     closest to the input value or its index value
#'
#' @param x the value to lookup
#' @param invect the vector in which to lookup the value x
#' @param index should the closest value be returned or its index; 
#'     default=TRUE
#'
#' @return by default it returns the index in the vector of the value 
#'     closest to the input value
#' @export
#'
#' @examples
#' \dontrun{
#' vals <- rnorm(100,mean=5,sd=2)
#' pick <- which.closest(5.0,vals,index=TRUE)
#' pick
#' vals[pick]
#' which.closest(5.0,vals,index=FALSE)
#' }
which.closest <- function(x,invect,index=T) {
  pick <- which.min(abs(invect-x))
  if (index) {
    return(pick)
  } else {
    return(invect[pick])
  }
} # end of which_.closest

#' @title wtedmean calculates the weighted mean of a set of values and weights
#'
#' @description wtedmean solves the problem of calculating a weighted mean
#'     value from a set of values with different weights. Within the aMSE this
#'     is common when trying to summarize across populations within an SAU or
#'     summarize SAU within a zone by finding a mean value weighted by the
#'     respective catch from each related population or SAU.
#'
#' @param x the values whose weighted mean is wanted
#' @param wts the weights to use, often a set of catches
#'
#' @return a single real number
#' @export
#'
#' @examples
#' saucpue <- c(91.0,85.5,88.4,95.2)
#' saucatch <- c(42.0,102.3,75.0,112.0)
#' wtedmean(saucpue,saucatch)
#' saucatch/sum(saucatch)  # the relative weights
wtedmean <- function(x,wts) {
  pwts <- wts/sum(wts,na.rm=TRUE)
  ans <- sum((x * pwts),na.rm=TRUE)
  return(ans)
} # end of wtedmean



# fileutils -----------------------------------------------

#' @title describefunctions lists all R functions in a set of files
#' 
#' @description describefunctions lists all the R functions in a set of R files
#'     along with their syntax, the linenumber in each file, the filename, the
#'     function name, and the functions within the set of R files that each 
#'     function calls. In addition, there is now a crossreference column,
#'     which identifies which functions call each function. If just the indir
#'     is provided then all R files in that directory will be examined. .Rmd
#'     files will not be considered but any other file type starting with .R
#'     may cause trouble until I find a fix!
#'
#' @param indir the directory in which to find the R files
#' @param files a vector of filenames, as character, within which to search for 
#'     functions, default="", which means all R files in indir will be used
#' @param outfile the full path and name of the CSV file to which the results 
#'     should be saved. default="", which means the output will only be 
#'     returned invisibly. If outfile has a fullpath csv filename then it will
#'     also be written to that file as well as retunred invisibly
#'
#' @return It can produce a csv file but also returns the results invisibly 
#' @export
#'
#' @examples
#' filen <- tempfile("test",fileext=".R")
#' txt <- c("# this is a comment",
#'  "#' @title ...",
#'  "dummy <- function() {",
#'  "  out <- anotherdummy()",
#'  "  return(out)",
#'  "}",
#'  "# a possibly confusing use of function",
#'  "#' @title ...",
#'  "anotherdummy <- function() {",
#'  "  return(NULL)",
#'  "}",
#'  "  ")
#'  write(txt,file=filen)
#'  usedir <- paste0(tempdir(),"//")
#'  filename <- tail(unlist(strsplit(filen,"\\",fixed=TRUE)),1)
#'  x <- describefunctions(indir=usedir,files=filename,outfile="")
#'  x
describefunctions <- function(indir,files="",outfile="") {
  if (nchar(files[1]) == 0) {
    dirfiles <- dir(indir)
    pickfiles <- grep(".R",dirfiles,ignore.case=TRUE)
    files <- dirfiles[pickfiles]
    pickRmd <- grep(".Rmd",dirfiles,ignore.case=TRUE)
    if (length(pickRmd) > 0) files <- files[-pickRmd]
  }
  nfiles <- length(files)
  numfuns <- matrix(0,nrow=nfiles,ncol=1,dimnames=list(files,c("nfuns")))
  allfiles <- NULL
  for (i in 1:nfiles) { # i = 1
    outfuns <- listfuns(paste0(indir,files[i]))
    if (nrow(outfuns) > 1) {
      numfuns[i,1] <- nrow(outfuns)
    } else {
      if (nchar(outfuns[1,"functions"]) > 0) numfuns[i,11] <- 1
    }
    allfiles <- rbind(allfiles,outfuns)
  }
  allfilesort <- allfiles[order(allfiles[,"functions"]),]
  allrefs <- matrix(0,nrow=0,ncol=1)
  for (i in 1:nfiles) {# i = 1
    if (numfuns[i] > 0) {
       allrefs <- rbind(allrefs,findfuns(indir,files[i],
                                         allfilesort[,"functions"]))
    } else {
      allrefs <- rbind(allrefs,", , ")
    }
  }
  allfiles[,"references"] <- allrefs
  x <- allfiles[order(allfiles[,"functions"]),]
  x[,"crossreference"] <- ""
  nfun <- nrow(x)
  for (i in 1:nfun) {
    pickf <- grep(x[i,"functions"],x[,"references"])
    if (length(pickf) > 0) {
      if (length(pickf) == 1) x[i,"crossreference"] <- x[pickf,"functions"]
    } else {
      x[i,"crossreference"] <- paste0(x[pickf,"functions"],collapse=", ")
    }
  }
  if (nchar(outfile) > 5) write.csv(x,file = outfile)
  return(invisible(x))
} # end of describefunctions

#' @title extractpathway traces the sequence of functions calls within a function
#'
#' @description extractpathway is used when documenting the sequence of function
#'     calls within a set of functions within a package one is developing. It
#'     needs to know the location of the R directory (indir) for the package,
#'     the starting functions at the beginning of a particular algorithm, and
#'     a listing from the rutilsMH function describefunctions. Then it traces
#'     the sequential usage of all known package functions, ignoring base R
#'     functions. The final output is a vector of function names, starting with
#'     the top-level function.
#'
#' @param indir the R directory within an R package being documented
#' @param infun the name of the top-level function whose sequential function use
#'     is being explored
#' @param allfuns the output of applying readLines to a text file containing
#'     R code.
#'
#' @return a vector of function names in the sequence in which they are used
#'     within the first names function
#' @export
#'
#' @examples
#' print("Wait on complex use of tempDir; see describefunctions for an example")
extractpathway <- function(indir,infun,allfuns) {  #  indir=tempdir(),infun="dummy",allfuns=x
  functions <- allfuns[,"functions"]
  pickrow <- which(functions == infun)
  if (length(pickrow) == 0)
    stop("Input function to extractpathway, not in allfuns. /n")
  rfile <- allfuns[pickrow,"file"]
  infile <- paste0(indir,rfile,".R")
  solution <- infun
  content <- readLines(con=infile)
  funLines <- identifyfuns(content=content)
  begin <- grep(paste0(infun," <- function"),content)
  finish <- funLines[which(funLines == begin) + 1] - 1
  funcont <- content[(begin+1):finish]
  testhash <- substr(funcont,1,4)
  omit <- grep("#",testhash)
  if (length(omit > 0)) funcont <- funcont[-omit]
  funcont <- removeEmpty(funcont)
  nline <- length(funcont)
  for (i in 1:nline) { #   i =2
    txt <-  removeEmpty(unlist(strsplit(funcont[i],split=character(0))))
    bracket <- match("(",txt)
    if (!is.na(bracket)) {
      loc2 <- match("(",txt)[1] - 1 # get end of object name
      loc1 <- grep("<",txt) + 2        # get putative start of object name
      if (length(loc1) == 0) {
        fun <- paste0(txt[1:loc2],collapse="")
      } else {
        if (loc1 > loc2) {
          fun <- paste0(txt[1:loc2],collapse="")
        } else {
          fun <- paste0(txt[loc1:loc2],collapse="")
        }
      }
      if (fun %in% functions) solution <- c(solution,fun)
    }
  }
  return(unique(solution))
} # end of extractpathway


#' @title extractRcode pulls out the r-code blocks from Rmd files
#' 
#' @description extractRcode pulls out the r-code blocks from Rmd files and 
#'     saves them into a separate R file. 
#'
#' @param indir the directory in which the rmd file is to be found and into
#'     which the output file will be placed.
#' @param rmdfile the name of the Rmd file whose R code is to be extracted
#' @param filename the name of the R file into which the r-code is to go. 
#'
#' @return generates an R file in the working directory, otherwise returns nothing
#' @export
#'
#' @examples
#' print("wait on a real example")
extractRcode <- function(indir,rmdfile,filename="out.R") { # indir=indir; rmdfile=inrmd; filename="out.R"
  infile <- paste0(indir,"/",rmdfile)
  fileout <- paste0(indir,"/",filename)
  cat("# R-code from the file ",rmdfile,"\n\n",file=fileout,append=FALSE)
  txt <- readLines(infile)
  pick <- grep("```",txt)
  steps <- length(pick)
  for (i in seq(1,steps,2)) {
    begin <- pick[i]
    cat("#",txt[begin],"\n",file=fileout,append=TRUE)
    finish <- pick[i+1]
    for (j in (begin+1):(finish-1)) {
      cat(txt[j],"\n",file=fileout,append=TRUE)
    }
    cat("#",txt[finish],"\n\n\n\n",file=fileout,append=TRUE)
  }
} # end of extractRcode


#' @title findfuns finds references to other functions within other functions
#' 
#' @description findfuns is used when developing a complex project containing 
#'     many R files, each containing many R functions. Given a file that 
#'     contains a set of functions (infile) and a data.frame of all functions  
#'     from the project (allfuns), which is obtained using listfuns, then 
#'     findfuns searches each function for references to any of the projects
#'     functions. This allows them to be cross referenced
#'
#' @param indir the directory in which the file identified in 'infile' is
#'     located
#' @param infile the filename of the R file within which to search for the 
#'     functions listed in the allfuns data.frame derived from the listfuns
#'     function
#' @param allfuns a data.frame of functions and their properties listed in 
#'     the order of the sorted function names in the 'function' column 
#'
#' @return the same data.frame except that the references column will have been
#'     populated
#' @export
#'
#' @examples
#' print("wait on suitable data-set")
findfuns <- function(indir,infile,allfuns) { # indir=indir;infile=files[i]; allfuns=allfilesort[,"functions"]
  # indir=ddir;infile=files[1]; allfuns=allfilesort[,"functions"]
  infile <- file.path(indir,infile)
  numfun <- length(allfuns)
  content <- readLines(con=infile)
  rfun <- tail(unlist(strsplit(infile,"/")),1)
  rfile <- substr(rfun,1,nchar(rfun)-2)
  funLines <- grep("function",content)
  titles <- grep("@title",content)
  testhash <- substr(content[funLines],1,4)
  omit <- grep("#",testhash)
  if (length(omit) > 0) {
    funLines <- funLines[-omit]
    testhash <- testhash[-omit]
  }
  omit2 <- grep("  ",testhash) # remove functions internal to other functions
  if (length(omit2) > 0) funLines <- funLines[-omit2]
  nfun <- length(funLines)
  outf <- as.data.frame(matrix("",nrow=nfun,ncol=1))
  bounds <- matrix(0,nrow=nfun,ncol=2,
                   dimnames=list(paste0(rfile,1:nfun),c("start","end")))
  bounds[,1] <- funLines + 1
  if (nfun > 1) {
    bounds[,2] <- c((titles[2:nfun] - 2),length(content))
  } else {
    bounds[,2] <- length(content)
  }
  for (i in 1:nfun) { # i=1
    funname <- removeEmpty(unlist(strsplit(content[funLines[i]],"<-"))[1])
    funcont <- content[bounds[i,1]:bounds[i,2]]
    testhash <- substr(funcont,1,5)
    omit <- grep("#",testhash)
    if (length(omit) > 0) funcont <- funcont[-omit]
    whichfun <- ", "
    for (j in 1:numfun) {  #  j = 65
      if (allfuns[j] != funname)
        if (length(grep(allfuns[j],funcont)) > 0)
          whichfun <- paste0(whichfun,allfuns[j],", ")
    }
    outf[i,] <- whichfun
  }
  return(outf)
} # end of findfuns


#' @title getDBdir identifies the DropBox path
#'
#' @description getDBdir is needde where multiple computers have different 
#'    names.
#'
#' @return the path to the DroBox directory
#' @export
#'
#' @examples
#' getDBdir()
getDBdir <- function() {
  if (dir.exists("C:/Users/Malcolm/Dropbox")) {
    prefixdir <- "C:/Users/Malcolm/Dropbox/"
  } else { 
    if (dir.exists("C:/Users/had06a/DropBox")) {
      prefixdir <- "C:/Users/had06a/DropBox/" 
    } else {
      prefixdir <- "C:/Users/User/Dropbox/"
    }
  }
  return(prefixdir)
} # end of getDBdir


#' @title getname returns the name of a variable as character
#'
#' @description getname runs 'deparse(substitute(x))' to get the
#'     name of the input variable. Saves remembering the syntax
#'
#' @param x any variable whose name is wanted as a character string
#'
#' @return a character string with the name of input variable
#' @export
#'
#' @examples
#' \dontrun{
#' a_variable <- c(1,2,3,4,5,6,7,8)
#' getname(a_variable)
#' }
getname <- function(x) {
  return((deparse(substitute(x))))
}

#' @title getnamespace returns the namespace for a given function
#'
#' @description getnamespace searches the loaded NameSpaces and returns the
#'     name of the NameSpace or package for the input function. This is used
#'     in by 'network'. If the namespace is not loaded this will not be able
#'     to be found.
#'
#' @param fun the name of the function of interest. It must be of class
#'     character, which can be obtained using 'getname'
#'
#' @return the name of the loaded NameSpace or package within which a 
#'     function can be found.
#' @export
#'
#' @examples
#' \dontrun{
#'    getnamespace(getname(lm))
#'    getnamespace(getname(anova))
#' }
getnamespace <- function(fun) {
  if (nchar(fun) == 0) return(NA)
  nss <- loadedNamespaces()
  envs <- c(lapply(nss,.getNamespace))
  return(nss[vapply(envs, function(env) exists(fun, env, inherits = FALSE),logical(1))])
} # end of get_namespace

#' @title identifyfuns uses text from readLines to identify function beginnings
#' 
#' @description identifyfuns is used when tracing the interactions between 
#'     functions within R packages. It uses the vector of character vectors
#'     that is produced by readLines and identifies the starting lines of all
#'     functions. It ignores all functions defined within comments, as well as
#'     ignoring all functions defined internally to other functions. It does the
#'     latter by testing for a couple of spaces at the start of a line 
#'     containing a function definition, which functions defined within another
#'     function should have.
#'
#' @param content the output of applying readLines to a text file containing 
#'     R code.
#'
#' @return a vector of line numbers identifying the start of all functions 
#'     within the content. This may be a vector of zero length if there are no
#'     functions.
#' @export
#'
#' @examples
#' txt <- c("# this is a comment",
#' "dummy <- function() { return(NULL) }",
#' "# a possibly confusing use of function",
#' "anotherdummy <- function() { return(NULL) }")
#' identifyfuns(txt)
identifyfuns <- function(content) {
  funLines <- grep("function",content)
  testhash <- substr(content[funLines],1,4)
  omit <- grep("#",testhash)
  if (length(omit) > 0) {
    funLines <- funLines[-omit]
    testhash <- testhash[-omit]
  }
  omit2 <- grep("  ",testhash) # remove functions internal to other functions
  if (length(omit2) > 0) funLines <- funLines[-omit2]
  return(funLines)
} # end of identifyfuns

#' @title pkgfuns names all functions within a package
#'
#' @description pgkfuns when given the name of a loaded library gives the 
#'     names of all functions within that library sorted in alphebetical 
#'     order.
#'
#' @param packname the name of the package as character
#'
#' @return a character vector containing the names of all functions in the 
#'     named package
#' @export
#'
#' @examples
#' \dontrun{
#'   pkgfuns("graphics")
#'   pkgfuns("rutilsMH")
#' }
pkgfuns <- function(packname) { # packname=pkgname
  funcs <- names(.getNamespace(packname))
  pick <- grep("__",funcs)
  funcs <- funcs[-pick]
  pick <- which(funcs == ".packageName")
  funcs <- funcs[-pick]
  return(sort(funcs))
} # end of pgkfuns

#' @title splitDate - Generates a vector of date and time components
#'
#' @description splitDate - Generates a vector of date and time components,
#'     perhaps for inclusion in filenames or other labels; helpful for
#'     keeping different run outputs seperate and identifiable.
#' @param dat - a system time from Sys.time() to be broken in components;
#'     defaults to NA, whereupon the current time is used.
#' @return a vector od characters relating to 'Year', 'Month', 'Day','Time',
#'     and a DateTime, which is a combination of all of these suitable for
#'     inclusion in a filename.
#' @export
#' @examples
#' \dontrun{
#' tmp <- splitDate()
#' print(tmp)
#' print(names(tmp))
#' print(as.numeric(tmp[1:3]))
#' print(tmp["DateTime"])
#' }
splitDate <- function(dat=NA) {
  if(is.na(dat)) dat <- as.POSIXlt(Sys.time())
  out <- unlist(dat)
  tim <- paste(trunc(as.numeric(out[3])),trunc(as.numeric(out[2])),
               "_",trunc(as.numeric(out[1]),1),sep="")
  day <- as.character(trunc(as.numeric(out[4])))
  month <- as.character(trunc(as.numeric(out[5])) + 1)
  if (nchar(month) == 1) month <- paste(0,month,sep="")
  year <- as.character(trunc(as.numeric(out[6])) - 100)
  combined <- paste(year,month,day,"_",tim,sep="")
  ans <- c(year,month,day,tim,combined)
  names(ans) <- c("Year","Month","Day","Time","DateTime")
  return(ans)
} # end of split_Date

# rmdutils------------------------------------


#' @title digitsbyrow a helper function for knitr, to specify formats by row
#'
#' @description digitsbyrow is a solution obtained from StackOverFlow, suggested
#'     by Tim Bainbridge in 11/12/19. knitr formats table columns as a whole,
#'     which can be a problem if one wants to mix integers with real numbers in
#'     the same columns. This first transposes the data.frame/matrix being
#'     printed, fixes the formats, and then transposes it back. In knitr one
#'     then needs to use the align argument to fix the alignment. In may version
#'     I have conserved both rownames and colnames for both data.frames and
#'     matrices (the original only did so for data.frames but I often print
#'     matrices). digitsbyrow converts all entries to character so knitr becomes
#'     necessary for printing.
#'
#' @param df the data.frame or matrix to be printed by knitr
#' @param digits a vector of the digits wanted for each row of the df or matrix
#'
#' @return a formatted data.frame or matrix depending on input
#' @export
#'
#' @examples
#' x <- matrix(c(rnorm(5,mean=5,sd=1),seq(1,10,1)),nrow=3,ncol=5,byrow=TRUE,
#'             dimnames=list(1:3,1:5))
#' digitsbyrow(x, c(3,0,0))
#' # needs knitr to use kable
#' # kable(digitsbyrow(x, c(3,0,0)),align='r',row.names=TRUE)
digitsbyrow <- function(df, digits) {
  tmp0 <- data.frame(t(df))
  tmp1 <- mapply(
    function(df0, digits0) {
      formatC(df0, format="f", digits=digits0)
    },
    df0=tmp0, digits0=digits
  )
  tmp1 <- data.frame(t(tmp1))
  rownames(tmp1) <- rownames(df)
  colnames(tmp1) <- colnames(df)
  if (class(df)[1] == "matrix") tmp1 <- as.matrix(tmp1)
  return(tmp1)
} # end of digitsbyrow

#' @title halftable halves the height of a tall narrow data.frame
#'
#' @description halftable would be used when printing a table using kable
#'     from knitr where one of the columns was Year. The objective would be 
#'     to split the table in half taking the bottom half and attaching it on
#'     the right hand side of the top half. The year column would act as the
#'     index.
#'
#' @param inmat the data.frame to be subdivided
#' @param yearcol the column name of the year field
#' @param subdiv the number of times the data.frame should be subdivided;
#'     the default is 3 but the numbers can only be 2 or 3.
#'
#' @return a data.frame half the height and double the width of the original
#' @export
#'
#' @examples
#' \dontrun{
#' x <- as.data.frame(matrix(runif(80),nrow=20,ncol=4))
#' x[,1] <- 1986:2005
#' x[,4] <- paste0("text",1:20)
#' halftable(x,yearcol="V1",subdiv=2)
#' halftable(x[,c(1,2,4)],yearcol="V1")
#' x1 <- rbind(x,x[1,])
#' x1[21,"V1"] <- 2006
#' halftable(x1,yearcol="V1",subdiv=3)
#' }
halftable <- function(inmat,yearcol="Year",subdiv=3) {
  if (!(subdiv %in% c(2,3))) stop("\n subdiv must be 2 or 3 \n")
  numrow <- dim(inmat)[1]
  numcol <- dim(inmat)[2]
  extra <- rep(NA,numcol)
  if ((numrow %% subdiv) == 0) {
    newnr <- numrow/subdiv
    incomplete <- FALSE
  } else {
    newnr <- trunc(numrow/subdiv) + 1
    incomplete <- TRUE
  }
  # years <- inmat[,yearcol]
  first <- inmat[1:newnr,]
  if (subdiv == 2) {
    second <- inmat[-c(1:newnr),]
    diff <- (nrow(first) - nrow(second))
    if (diff > 0) {
      numcol <- ncol(inmat)
      third <- rbind(second,extra)
    } else {
      third <- second
    }
  } else {
    second <- inmat[c(newnr+1):c(2*newnr),]
    first <- cbind(first,second)
    third <- inmat[-c(1:(2*newnr)),]
    diff <- nrow(first) - nrow(third)
    if (diff > 0) third <- rbind(third,extra)
    if (diff > 1) third <- rbind(third,extra)
  }
  outmat <- cbind(first,third)
  rownames(outmat) <- 1:newnr
  return(outmat)
} # end of halftable





#' @title kablerow a replacement for knitr::kable which enables row formatting
#'
#' @description knitr::kable enables one to round the number of digits for each
#'     column of a table. However, sometimes one wants to format the rows and
#'     not the columns. kablerow enables that while using the kable function.
#'     It rounds the rows to the desired number of digits and then converts
#'     those rounded values to characters, which kable can then print more
#'     appropriately.
#'
#' @param x an input matrix or data.frame
#' @param rowdigits the number of digits desired for each row
#' @param namerows should row.names be printed; default=NA. change to TRUE for
#'     row.names printing
#' @param namecols should col.names be printed; default=NA (which prints V1, V2
#'     ,V3, ...)
#'
#' @return Nothing but it does use knitr::kable to print a formatted matrix
#' @export
#'
#' @examples
#' x <- matrix(rnorm(25,mean=5,sd=1),nrow=5,ncol=5)
#' colnames(x) <- 1:5
#' numdig <- c(2,3,4,3,2)
#' rownames(x) <- c("a","b","c","d","e")
#' kablerow(x,rowdigits=c(2,3,4,3,2),namerows=TRUE)
kablerow <- function(x,rowdigits,namerows=NA,namecols=NA) { # x=x; rowdigits=c(2,3,4,3,2); namerows=FALSE
  xr <- as.data.frame(x)
  num <- nrow(x)
  for (i in 1:num) {
    x[i,] <- round(x[i,],rowdigits[i])
    xr[i,] <- as.character(x[i,])
  }
  kable(xr,align="r",row.names=namerows,col.names=namecols)
} # end of kablerow


#' @title listExamples lists all the examples in a package R file
#'
#' @description listExamples lists all the examples in a package R file. It
#'     comments out the first line number and any dontrun statements along 
#'     with their following curly bracket.
#'
#' @param infile - a character variable containing the path and filename
#' @return Creates an R file in the working directory and prints its name to
#'     the console
#'
#' @export
#' @examples
#' \dontrun{
#' txt <- vector("character",4)
#' txt[1] <- "#' @examples "
#' txt[2] <- "#' /dontrun{"
#' txt[3] <- "#' print("This is an example of using listExamples")"
#' txt[4] <- "#' }"
#' infile <- textConnection(txt)
#' listExamples(infile)
#' }
listExamples <- function(infile) {  
  outfile <- paste0("examples_",tail(unlist(strsplit(infile,"/")),1))
  cat("All the example code from  \n",file=outfile,append=FALSE)
  cat(infile,"\n\n",file=outfile,append=TRUE)
  content <- readLines(con=infile)
  egLines <- grep("@examples",content)
  funlines <- grep("<- function",content)
  nline <- length(egLines)
  for (i in 1:nline) {
    # find extent of example  i = 1
    count <- 1
    cat("#Linenumber: ",egLines[i] + count,"\n",file=outfile,append=TRUE)
    repeat {  # i=1; count=1
      tmpline <- content[egLines[i] + count]
      if (substr(tmpline,1,1) == "#") {
        lenc <- nchar(tmpline)
        if ((length(grep("dontrun",tmpline)) > 0) |
            (length(grep("#' }",tmpline)) > 0)) {
          cat("# ",tmpline,file=outfile,append=TRUE)
        }
        cat(substr(tmpline,3,lenc),"\n",file=outfile,append=TRUE)
        count <- count + 1
      } else {
        cat("# ",tmpline,"\n",file=outfile,append=TRUE)
        cat("\n\n\n\n",file=outfile,append=TRUE)
        break()
      }
    }
  }
  print(outfile)
} # end of list_Examples

#' @title lininterpol - linearly interpolate values in a vector with NAs
#'
#' @description lininterpol - linearly interpolate values in a vector with 
#'     NAs. A common problem when plotting up time series is where there are
#'     missing values or NAs the plotted line will have gaps, one can always
#'     plot points on top of a line to identify where there are missing 
#'     values but an alternative would be to interpolate the missing values 
#'     linearly and plot that line as a dashed line. This function generates
#'     those linear interpolations. The input vector cannot have missing 
#'     values at the beginning or the end. If there are no missing values 
#'     the original vector is returned
#'
#' @param invect - the vector of values including missing values
#'
#' @return invect but with NAs replaced with linearly interpolated values.
#' @export
#'
#' @examples
#' \dontrun{
#'  Expt <- c(20102,18465,16826,15333,14355,NA,13843.7,NA,NA,NA,15180)
#'  lininterpol(Expt)
#' }
lininterpol <- function(invect) { 
  npt <- length(invect)
  answer <- invect
  pickNA <- which(is.na(invect))
  nna <- length(pickNA)
  if (nna == 0) return(invect)  # no NAs
  if ((pickNA[1] == 1) | (pickNA[nna] == nna))
    #   picknNA <- which(invect > 0)
    stop("input vector in lin-interpol cannot start or end with an NA")
  # identify groups of NAs
  group <- c(pickNA[1])
  count <- 1
  ans <- vector("list",nna) # possible each NA is an individual
  for (i in 2:nna) {
    if ((pickNA[i] - pickNA[(i-1)]) > 1) {
      ans[[count]] <- group
      group <- pickNA[i]
      count <- count + 1
    } else {
      group <- c(group,pickNA[i])
    }
  }
  ans[[count]] <- group
  for (i in 1:count) {  # i <- 2
    pickNA <- ans[[i]]
    begin <- (pickNA[1] - 1)
    finish <- (tail(pickNA,1) + 1)
    first <- invect[begin]
    second <- invect[finish]
    answer[begin:finish] <- seq(first,second,length=(length(pickNA) + 2))
  }
  return(answer)
}  # end of lin_interpol

#' @title listfuns produces a listing of all functions in an input R file
#'
#' @description listfuns reads in a given R file and then identifies each
#'     function header within it and pulls out the function name, its syntax,
#'     the line-number in the file, and associates that with the filename.
#'
#' @param infile the R file to be examined
#'
#' @return a data.frame of syntax, function name, line number, and file name
#' @export
#'
#' @examples
#' print("wait for an example")
listfuns <- function(infile) { # infile=paste0(ddir,filen[1]); 
  content <- readLines(con=infile)
  if (length(grep("/",infile) > 0)) {
    rfun <- tail(unlist(strsplit(infile,"/")),1)
    rfile <- substr(rfun,1,nchar(rfun)-2)
  } else {
    rfile <- infile
  }
  funLines <- grep("function",content)
  testhash <- substr(content[funLines],1,4)
  omit <- grep("#",testhash)
  if (length(omit) > 0) {
    funLines <- funLines[-omit]
    testhash <- testhash[-omit]
  }
  omit2 <- grep("  ",testhash) # remove functions internal to other functions
  if (length(omit2) > 0) funLines <- funLines[-omit2]
  nLine <- length(funLines)
  delF <- NULL
  if (nLine > 0) {
    for (i in 1:nLine) {
      tmpLine <- gsub(" ","",content[funLines[i]])
      if ((length(grep("function\\(",tmpLine)) == 0) |
          (substr(tmpLine,1,2) == "#'") |
          (length(grep("<-function",tmpLine)) == 0) |
          (length(grep("} #",tmpLine)) > 0)) delF <- c(delF,i)
    }
  }  
  ndelF <- length(delF)
  if (ndelF > 0) {
    funLines <- funLines[-delF]
  }
  if (ndelF == nLine) {
    txt <- paste0(infile,"  contained no recognizable functions")
    warning(cat(txt,"\n"))
    out <- "NA"
    funnames <- ""
    funLines <- 1
    n <- 1
  } else {
    outlines <- sort(c(funLines))
    out <- content[outlines]
    funnames <- out
    n <- length(out)
    for (i in 1:n) {  # i=1
      out[i] <- gsub(" ","",(unlist(strsplit(out[i],"\\{")))[1])
      funnames[i] <- removeEmpty(unlist(strsplit(out[i],"<-"))[1])
      out[i] <- gsub("<-function","",out[i])
    }
  }
  columns <- c("syntax","linenumber","file","functions","references")
  rows <- paste0(rfile,1:n)
  outfuns <- as.data.frame(matrix(NA,nrow=n,ncol=length(columns),
                                  dimnames=list(rows,columns)))
  outfuns[,"syntax"] <- out
  outfuns[,"functions"] <- funnames
  outfuns[,"linenumber"] <- funLines
  outfuns[,"file"] <- rfile
  return(outfuns)
} # end of listfuns


#' @title rmdcss generates some initial css style code for HTML Rmd files
#' 
#' @description rmdcss generates some initial css style code for HTML Rmd files
#'     as well as a mathjax script that will generate equation numbers for any
#'     display equations in the document. This prints the css style code and
#'     the mathjax script to the console from where it should be pasted into the
#'     Rmd file immediately following the YAML header. It now contains font 
#'     sizes for the h1 heaqder and the .inline and .display math classes
#'
#' @return nothing but it prints css style code and a mathjax script to the 
#'     console
#' @export
#'
#' @examples
#' rmdcss()
rmdcss <- function() {
  cat('<style type="text/css"> \n',
      '  body, td { \n',
      '  font-size: 16px; \n',
      '  font-family: "Times New Roman", Times, serif; \n',
      '} \n')
  cat('code.r{ \n',
      '  font-size: 15px; \n',
      '} \n')
  cat('pre {  \n',
      'font-size: 8px  \n',
      '}  \n')
  cat('h1 {  \n',
      '  font-size: 32px  \n',
      '}  \n')
  cat('.inline{font-size: 15px; } \n',
      '.display{font-size: 18px;} \n',
      '<','/style>  \n')
  cat('\n\n')
  cat('<script type="text/x-mathjax-config">  \n',
      '  MathJax.Hub.Config({  \n',
      '    TeX: {   \n',
      '      equationNumbers: {   \n',
      '        autoNumber: "all",  \n',
      '        formatNumber: function (n) {return ',3.,'+n}  \n',
      '      }  \n', 
      '    }  \n',
      '  });  \n',
      '<','/script>  \n')
} # end of rmdcss

#' @title setuprmd sets up and Rmd file ready to generate an HTML file
#' 
#' @description setuprmd sets up a custom Rmd file for generating an HTML file,
#'     which better suits my own preferences
#'
#' @param filen the full path filename for the final Rmd file. Ensure its 
#'     filetype = .Rmd. The default = "", which write the custom text to the 
#'     console.
#'
#' @return nothing but it does write a file to one's hard drive in the location
#'    listed in filen
#' @export
#'
#' @examples
#' setuprmd(filen="")
setuprmd <- function(filen="") {
  cat('--- \n',
      'title: "Title" \n',
      'author: Malcolm Haddon \n',
      'date: "`r Sys.time()`"  \n',
      'output: \n',
      '  html_document:   \n',
      '    df_print: paged    \n',
      '    fig_caption: yes \n',
      '    fig_height: 5.5 \n',
      '    fig_width: 6.5\n',
      '    number_section: yes \n',
      # '    toc: yes \n',
      # '    toc_depth: 2 \n',
      '---  \n',
      sep = "", file=filen, append=FALSE)
  cat('  \n',
      '```{r setup, include=FALSE}  \n',
      'knitr::opts_chunk$set(  \n',
      '  echo = FALSE,  \n',
      '  message = FALSE,  \n',
      '  warning = FALSE)  \n\n',
      'options(knitr.kable.NA = "", \n',
      '        knitr.table.format = "pandoc")  \n',
      '  \n\n',
      sep = "", file=filen, append=TRUE)
  cat('options("show.signif.stars"=FALSE,  \n',
      '        "stringsAsFactors"=FALSE,   \n',
      '        "max.print"=50000,          \n',
      '        "width"=240)                \n',
      '```  \n\n',
      sep = "", file=filen, append=TRUE)
  # cat('   \n',
  #     '<style type="text/css">  \n',
  #     '   body, td, h1, h2, h3, h4 {  \n',
  #     '   font-size: 16px;   \n',
  #     '   font-family: "Times New Roman" Times, serif; \n',
  #     '}   \n',
  #     '< /style>   \n\n\n',
  #     sep = "", file=filen, append=TRUE)
} # end of setuprmd

# diagrams ------------------------------------------------

#' @title circle draws a circle with a given origin and radius
#' 
#' @description circle provides the means of drawing a circle of a given
#'     radius and origin within a diagram ready for the addition of text.
#'
#' @param origx the final x origin
#' @param origy the final y origin
#' @param radius the radius of the circle
#' @param col the col of the circle
#' @param lwd the line width of the circle
#'
#' @return the matrix of x and y values invisibly  
#' @export
#'
#' @examples
#'   makecanvas()
#'   circle(origx=35,origy=70,radius=30,lwd=2,col=1)
#'   circle(origx=65,origy=60,radius=30,lwd=2,col=2)
#'   circle(origx=45,origy=40,radius=30,lwd=2,col=4)
circle <- function(origx=50,origy=50,radius=10,col=1,lwd=1) {
  ans <- pol2cart(angle=seq(0,360,0.1),dist=radius,xorig=origx,yorig=origy)
  lines(ans[,"x"],ans[,"y"],lwd=lwd,col=col)
  return(invisible(ans))
} # end of circle

#' @title cart2pol converts cartesian coordinates into the polar angle
#' 
#' @description cart2pol as a step in converting cartesian coordinates into
#'     polar coordinates this calculates the angle, in degrees, from x y
#'     values
#'
#' @param x either a vector of two values of a matrix of pairs of values
#'
#' @return a single angle of vector of angles
#'
#' @examples
#' \dontrun{
#'   cart2pol(c(3,3))  # should be 45
#'   dat <- matrix(c(3,4,5,7),nrow=2,ncol=2,byrow=TRUE)
#'   print(dat)
#'   cart2pol(dat)     # should be 36.8699 twice.
#' }
cart2pol <- function(x){
  if (is.vector(x)) angle <- 180 * (atan2(x[1],x[2])) / pi
  if (is.matrix(x)) angle <- 180 * (atan2(x[,1],x[,2])) / pi
  return(angle=angle)
} # end of cart2pol

#' @title diagrams provides the syntax of functions for making diagrams
#' 
#' @description diagrams provides the syntax of functions for making diagrams
#'
#' @return nothing but it write syntax for diagram functions to the console
#' @export
#'
#' @examples
#' diagrams()
diagrams <- function() {
  cat('circle(origx = 50, origy = 50, radius = 10, col = 1, lwd = 1) \n')
  cat('makecanvas(xstart = 0, xfinish = 100, ystart = 0, yfinish = 100) \n')
  cat('makerect(left, xinc, top, yinc, linecol = "grey", lwd = 1) \n')
  cat('makevx(init, inc) \n')
  cat('makevy(init, inc) \n')
  cat('plotoblong(x0, x1, y0, y1, border = 1, col = 0, lwd = 1)  \n')
} # end of diagrams

#' @title makecanvas sets up a plotting area ready for the flowchart
#'
#' @description makecanvas sets up a plotting areas ready for a flowchart
#'     made up of shapes, circles, polygons, rectangles, text, and arrows
#'
#' @param xstart x-origin value defaults = 0
#' @param xfinish maximum of x axis defaults = 100
#' @param ystart y-origin value default = 0
#' @param yfinish y-axis maximum default = 100
#'
#' @return nothing but plots an empty graph ready for polygons and text
#' @export
#'
#' @examples
#' \dontrun{
#'   makecanvas(ystart=50,yfinish=93.5)
#'   polygon(makevx(2,27),makevy(90,6),col=0,lwd=1,border=1)
#' }
makecanvas <- function(xstart=0,xfinish=100,ystart=0,yfinish=100) {
  par(mfrow=c(1,1),mai=c(0.1,0.1,0.1,0.1),oma=c(0.0,0,0.0,0.0))
  par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
  plot(seq(xstart,xfinish,length=101),seq(ystart,yfinish,length=101),
       type="n",xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
} # end of makecanvas

#' @title makevx make an x values vector
#'
#' @description makevx takes the left x value of a rectangle and the
#'     increment rightwards that defines a vector describing the four
#'     vertices of the rectangle topleft, topright, bottomright,
#'     bottomleft, topleft. when matched with makevy generates the
#'     descriptor for a complete rectangle.
#'
#' @param init x-value for the left-hand edge of a rectangle
#' @param inc the x-increment added to init to define the right-hand edge
#'
#' @return a vector of y-values
#' @export
#'
#' @examples
#' \dontrun{
#'  plot(0:100,seq(58,93.5,length=101),type="n",xaxt="n",yaxt="n",
#'  xlab="",ylab="", bty="n")
#'  polygon(makevx(2,27),makevy(90,6),col=0,lwd=1,border=1)
#' }
makevx <- function(init,inc) {
  return(c(init,init+inc,init+inc,init,init))
}


#' @title makevy make a y values vector
#'
#' @description makevy takes the top y value of a rectangle and the
#'     vertical increment downwards and defines a vector describing the four
#'     vertices of the rectangle topleft, topright, bottomright,
#'     bottomleft. topleft, when matched with makevx generates the
#'     descriptor for a complete rectangle.
#'
#' @param init y-value for the top edge of a rectangle
#' @param inc the y-increment subtracted from init to define the lower edge
#'
#' @return a vector of y-values
#' @export
#'
#' @examples
#' \dontrun{
#'  canvas(ystart=50,yfinish=93.5)
#'  polygon(makevx(2,27),makevy(90,6),col=0,lwd=1,border=1)
#' }
makevy <- function(init,inc) {
  return(c(init,init,init-inc,init-inc,init))
}


#' @title makerect draws a rectangle once a plot is available
#'
#' @description makerect draws a rectangle after canvas has been called
#'
#' @param left defines lefthand edge of rectangle
#' @param xinc left + xinc defines right-hand edge or rectangle
#' @param top defines top edge of rectangle
#' @param yinc top - yincdefines bottom edge of rectangle
#' @param linecol colour of line. default="grey"
#' @param lwd the width of the line, default=1
#' @param col the fill colour of the polygon drawn. default=NULL so not filled
#'
#' @return a vector denoting the center (x,y) of the rectangle
#' @export
#'
#' @examples
#' \dontrun{
#'    canvas(ystart=50,yfinish=93.5)
#'    makerect(left=2,xinc=27,top=90,yinc=6)
#' }
makerect <- function(left,xinc,top,yinc,linecol="grey",lwd=1,col=NULL) {
  polygon(makevx(left,xinc),makevy(top,yinc),col=col,
          lwd=lwd,border=linecol)
  centerx <- (left * 2 + xinc)/2
  centery <- (top * 2 - yinc)/2
  return(invisible(c(centerx,centery)))
}

#' @title plotoblong generates an oblong from x0,x1,y0,y1
#' 
#' @description plotoblong generates an oblong from x0,x1,y0,y1
#'
#' @param x0 x-axis left
#' @param x1 x-axis right
#' @param y0 yaxis bottom
#' @param y1 yaxis top
#' @param border colour of the border, default=black=1
#' @param col colour of fill, default = 0 =  empty
#' @param lwd width of the line,default=1
#'
#' @return nothing but it plots a polygon
#' @export
#'
#' @examples
#' \dontrun{
#'   canvas()
#'   plotoblong(1,50,1,50,lwd=3,linecol=4)
#' }
plotoblong <- function(x0,x1,y0,y1,border=1,col=0,lwd=1) {
  x <- c(x0,x0,x1,x1,x0); y <- c(y0,y1,y1,y0,y0)
  polygon(x,y,lwd=lwd,border=border,col=col)
}

#' @title pol2cart polar to cartesian coordinates
#' 
#' @description pol2cart translate polar coordinates of angles (as degrees)
#'     and a distance = radius, into cartesian coordinates of x and y. The
#'     option of using arbitrary origin coordinates is included
#'
#' @param angle the angle in degrees, either a single number of a vector
#' @param dist the length of the line or radius, a single number
#' @param xorig the final xorigin
#' @param yorig the final yorigin
#'
#' @return a matrix of 1 or more rows depending on length of angle
#' @export
#'
#' @examples
#' \dontrun{
#'   ans <- pol2cart(angle=seq(0,360,15),dist=20,xorig=30,yorig=30)
#'   print(ans)
#' }
pol2cart <- function(angle,dist,xorig=0,yorig=0){
  #  angle=45:50; dist=10; xorig=0; yorig=0
  numang <- length(angle)
  coord <- matrix(0,nrow=numang,ncol=2,dimnames=list(1:numang,c("x","y")))
  angler <- angle*pi/180
  for (i in 1:numang) {
    coord[i,] <- c(xorig + dist * sin(angler[i]),
                   yorig + dist * cos(angler[i]))  
  }
  return(coord) #output the new x and y coordinates
} # end of pol2cart

#' @title pythag calculates Pythagorus' theorum on a vector of two values
#' 
#' @description pythag Pythagorus' theorum states that the length of the
#'     hypotheneuse between two lines at right angels to each other (that
#'     is in cartesian coordinates) is the sqrt of the sum of their squares.
#'
#' @param x a vector of two numbers or a matrix of pairs of numbers
#'
#' @return a single number or a vector depending on input
#' @export
#'
#' @examples
#' \dontrun{
#'  pythag(c(3,4))  # should be 5
#'  dat <- matrix(c(3,4,5,7),nrow=2,ncol=2,byrow=TRUE)
#'  print(dat)
#'  pythag(dat)     # should be 5 and 10
#' }
pythag <- function(x) {  # x = ans
  if (is.vector(x)) ans <- sqrt((x[1]^2 + x[2]^2))
  if (is.matrix(x)) ans <- sqrt((x[,1]^2 + x[,2]^2))
  return(ans) 
}

#' @title '%ni%' identifies which element in x is NOT in y
#'
#' @param x a vector of elements which can be numeric or character
#' @param y a vector of elements which can be numeric or character
#'
#' @export
#' 
#' @examples
#'   x <- 1:10
#'   y <- 6:18
#'   x %ni% y
#'   pick <- (x %ni% y)
#'   x[pick]
`%ni%` <- function(x,y) {
  !(x %in% y)
}
