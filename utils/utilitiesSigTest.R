
binaryByAbsThreshold <- function( eids, css.timecourse, data.matrix=dm,cutoff=300){
  csst.label <- css.timecourse$name
  cols <- css.timecourse[["DM Column"]]
  nt <- length(cols)
  if ( length(eids) == 1 ){
    eid <- eids
    tc <- data.matrix[eid,cols]
    retval <- 0
    if ( max(tc) > cutoff ){ retval <- 1}
    names(retval) <- csst.label
  } else if ( length(eids) > 1 ){
    tcs <- data.matrix[eids,cols]
    maxabs <- apply(tcs,1,max)
    retval <- ( maxabs > cutoff ) + 0
  } else {
    stop("Error:Problem with Entrez IDs")
  }
  return (retval)
}

binaryByRatioThreshold <- function( eids, css.timecourse, data.matrix=dm,cutoff=3 ){
  csst.label <- css.timecourse$name
  cols <- css.timecourse[["DM Column"]]
  nt <- length(cols)
  if ( length(eids) == 1 ){
    eid <- eids
    tc <- data.matrix[eid,cols]
    maxrat <- max((tc/tc[1])[2:nt])
    retval <- 0
    if ( maxrat > cutoff ){ retval <- 1}
    names(retval) <- csst.label
  } else if ( length(eids) > 1 ){
    tcs <- data.matrix[eids,cols]
    ratmat <- (tcs/tcs[,1])[,2:nt]
    if ( nt > 2 ){
      maxrats <- apply(ratmat,1,max)
    } else { ## if only one ratio for t >0 
      maxrats <- ratmat
    }
    retval <- ( maxrats > cutoff ) + 0
  } else {
    stop("Error:Problem with Entrez IDs")
  }
  return (retval)
}

 
binarizeCSSTs <- function ( eids, cssts,data.matrix=dm ,ratio.cutoff=3, abs.cutoff=300 ) {
  if ( length(eids) == 1){
    rvec <- numeric()
    for ( csst in cssts ){    
      rval <- binaryByRatioThreshold(eid,csst,data.matrix=data.matrix,cutoff=ratio.cutoff)
      aval <- binaryByAbsThreshold(eid,csst,data.matrix=data.matrix,cutoff=abs.cutoff)
      rvec <- c(rvec,rval&aval )
    }
    return(rvec)
  }  else if ( length(eids) > 1 ){
    rmat <- numeric()
    for ( csst in cssts ){    
      rvals <- binaryByRatioThreshold(eids,csst,data.matrix=data.matrix,cutoff=ratio.cutoff)
      avals <- binaryByAbsThreshold(eids,csst,data.matrix=data.matrix,cutoff=abs.cutoff)
      combvals <- rvals & avals
      rmat <- cbind(rmat, combvals)
    }
    colnames(rmat) <- names(cssts)
    return(rmat)
  } else {
    stop("Error:Problem with Entrez IDs")
  }

}


maxRatioCSST <- function( eids, css.timecourse, data.matrix=dm ){
  csst.label <- css.timecourse$name
  cols <- css.timecourse[["DM Column"]]
  nt <- length(cols)
  if ( length(eids) == 1 ){
    eid <- eids
    tc <- data.matrix[eid,cols]
    maxrat <- max((tc/tc[1])[2:nt])
    retval <- maxval
    names(retval) <- csst.label
  } else if ( length(eids) > 1 ){
    tcs <- data.matrix[eids,cols]
    ratmat <- (tcs/tcs[,1])[,2:nt]
    if ( nt > 2 ){
      maxrats <- apply(ratmat,1,max)
    } else { ## if only one ratio for t >0 
      maxrats <- ratmat
    }
    retval <- maxrats
  } else {
    stop("Error:Problem with Entrez IDs")
  }
  return (retval)
}


maxRatioMultipleCSSTs <- function ( eids, cssts, data.matrix=dm ) {
  if ( length(eids) == 1){
    rvec <- numeric()
    for ( csst in cssts ){    
      rval <- maxRatioCSST(eid,csst,data.matrix=data.matrix)
    }
    return(rvec)
  }  else if ( length(eids) > 1 ){
    rmat <- numeric()
    for ( csst in cssts ){
      rvals <- maxRatioCSST(eids,csst,data.matrix=data.matrix)
      rmat <- cbind(rmat, rvals)
    }
    colnames(rmat) <- names(cssts)
    return(rmat)
  } else {
    stop("Error:Problem with Entrez IDs")
  }

}


meanAbsCSST <- function( eids, css.timecourse, data.matrix=dm ){
  csst.label <- css.timecourse$name
  cols <- css.timecourse[["DM Column"]]
  nt <- length(cols)
  if ( length(eids) == 1 ){
    eid <- eids
    tc <- data.matrix[eid,cols]
    meanabs <- mean(tc)
    retval <- meanabs
    names(retval) <- csst.label
  } else if ( length(eids) > 1 ){
    tcs <- data.matrix[eids,cols]
    meanabs <- apply(tcs,1,mean)
    retval <- meanabs
  } else {
    stop("Error:Problem with Entrez IDs")
  }
  return (retval)
}

meanAbsMultipleCSSTs <- function ( eids, cssts, data.matrix=dm ) {
  if ( length(eids) == 1){
    avec <- numeric()
    for ( csst in cssts ){    
      aval <- meanAbsCSST(eid,csst,data.matrix=data.matrix)
    }
    return(avec)
  }  else if ( length(eids) > 1 ){
    amat <- numeric()
    for ( csst in cssts ){
      avals <- meanAbsCSST(eids,csst,data.matrix=data.matrix)
      amat <- cbind(amat, avals)
    }
    colnames(amat) <- names(cssts)
    return(amat)
  } else {
    stop("Error:Problem with Entrez IDs")
  }

}

