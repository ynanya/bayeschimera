library(HDInterval)
maphdi <- function(fit,credMass=0.99,include=TRUE,pars=NA){
  if(credMass<0|credMass>1) print("warning : credMass values are not valid")
  Class <- class(fit)[1]
  if(Class=="stanfit"){
    e <- as.data.frame(rstan::extract(fit))
  }else if(Class=="rjags"){
    e <- as.data.frame(fit$BUGSoutput$sims.list)
  }else if(Class=="data.frame"){
    e <- fit
  }else{
    print("error : this object cannot be processed")
    stop()
  }

  if(is.na(pars[1])) pars <- colnames(e)
  if (include==TRUE){
    e <- e[,pars]
  }else{
    e <- e[,!colnames(e) %in% pars]
  }
  if(length(pars)==1){
    MAP <- density(e)$x[which.max(density(e)$y)]
    HDI95 <- hdi(e)
    HDIAny <- hdi(e,credMass = credMass)
    Result <- data.frame(MAP = MAP, Lower95 = HDI95[1], Upper95 = HDI95[2],
                         LowerAny = HDIAny[1],UpperAny = HDIAny[2])
    rownames(Result) <- pars
  }else{
    MAP <- apply(e,2,function(z) density(z)$x[which.max(density(z)$y)])
    HDI95 <- apply(e,2,HDInterval::hdi)
    HDIAny <- apply(e,2,function(x) HDInterval::hdi(x,credMass = credMass))
    Result <- data.frame(MAP = MAP, Lower95 = HDI95[1,], Upper95 = HDI95[2,],
                         LowerAny = HDIAny[1,],UpperAny = HDIAny[2,])
  }

  colnames(Result)[4:5] <- c(
    paste("Lower",round(credMass*100,2),sep=""),
    paste("Upper",round(credMass*100,2),sep="")
  )
  return(Result)
  }