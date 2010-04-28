
searchByNameValue <- function( termlistHttpForm ){
  datatosend <- paste("PATH=/sampleData/microarray/samples",
                      "MATCHING_ALL_TERMS=true",
                      "PROJECTION=raw_data_path",
                      "PROJECTION=organism",
                      "PROJECTION=Cell Type",
                      "PROJECTION=Strain",
                      "PROJECTION=Stimulus 1",
                      "PROJECTION=Time 1",
                      "PROJECTION=Stimulus 2",
                      "PROJECTION=Time 2",
                      "PROJECTION=chip",
                      "PROJECTION=Sex",
                      "PROJECTION=Investigator",
                      termlistHttpForm,
                      sep="&");
  
  sstring <- whiteSpaceURLform(paste("/addama-rest/primary-repo/search?",datatosend,sep=""))
  response <-  http.get("innate-immunity.systemsbiology.net",sstring,port=9080)
  
  ## 
  ## Last working, but slow version
  ##response = simplePostToHost("innate-immunity.systemsbiology.net", "/addama-rest/primary-repo/search", datatosend, referer="", port=9080);
  ## Port 80:Believe this version is preferable and general
  ##response = simplePostToHost("innate-immunity.systemsbiology.net", "/addama-rest/primary-repo/search", datatosend, referer="", port=80);
  ## Improvement pending
  ##response = http.post("innate-immunity.systemsbiology.net", "/addama-rest/primary-repo/search", datatosend, port=9080)
  
  ## chop HTTP header stuff in response
  resp = substr(response,regexpr("\\{",response)[1],nchar(response))

  resp = sub("}\r.*","}",resp)

  ##resp <- gsub("\r\n1d6f\r\n","",resp)
  resp <- gsub("\r\n[[:alnum:]]{1,50}\r\n","",resp)
  
  object = fromJSON(resp);

  return(object$result)
}

termHttpForm <- function(namedValue){
  return(paste(names(namedValue),namedValue,sep="="))
}
 
termListHttpForm <- function(termlist){
  return(paste(unlist(lapply(termlist,termHttpForm)),collapse="&"))
}

## Needs "bare" uuid
getNodeObject <- function ( uuid ){
  ##response <- getToHost("innate-immunity.systemsbiology.net",paste("/addama-rest/primary-repo/",uuid,sep=""),"",port=9080)
  response <- http.get("innate-immunity.systemsbiology.net",paste("/addama-rest/primary-repo/",uuid,sep=""),port=9080)
  ## chop HTTP header stuff in response
  resp <- substr(response,regexpr("\\{",response)[1],nchar(response))
  resp <- sub("}\r.*","}",resp)
  resp <- gsub("\r\n[[:alnum:]]{1,50}\r\n","",resp)
  object <- fromJSON(resp);
  return (object)
}  

## Needs full path from root
## whitespace OK
getNodeObjectByPath <- function ( path ){
  path <- whiteSpaceURLform(path)
  ##response <- getToHost("innate-immunity.systemsbiology.net",paste("/addama-rest/primary-repo?PATH=",path,sep=""),"",port=9080)
  response <-  http.get("innate-immunity.systemsbiology.net",paste("/addama-rest/primary-repo?PATH=",path,sep=""),port=9080)
  resp <- substr(response,regexpr("\\{",response)[1],nchar(response))
  resp <- sub("}\r.*","}",resp)
  resp <- gsub("\r\n[[:alnum:]]{1,50}\r\n","",resp)
  object <- fromJSON(resp);
  return(object)
}
  
meta <- function(stimulus1,time1,strain="Bl6",cell.type="Dendritic Cell",sex="Female",stimulus2=NULL,time2=NULL){
  outlist <- list()
  outlist[["Stimulus 1"]] <- stimulus1
  outlist[["Time 1"]] <- time1
  outlist[["Strain"]] <- strain
  outlist[["Cell Type"]] <- cell.type
  outlist[["Sex"]] <- sex
  outlist[["Stimulus 2"]] <- stimulus2
  outlist[["Time 2"]] <- time2
  return(outlist)
}


## Convert list metadata to string
## Strings match column headers on
## exon.form = TRUE: string form consistent with exon array output Data_Matrix.tsv
## exon.form = FALSE: string form consistent with Three Prime array output AllGenes.txt
## AllGenes.txt displays "Stimulus 1 Type", which is incorporated with a kludge below
metaToString <- function ( ob, exon.form=TRUE ){
  b <- unlist(ob)
  if (!("Stimulus 2" %in% names(b)) ){
    b["Stimulus 2"] <- ""
  }
  
  if (!("Time 2" %in% names(b)) ){
    b["Time 2"] <- ""
  }
  b2 <- sapply(b,dashify)
  if ( exon.form ){
    return(paste(b2[c("Cell Type","Strain","Stimulus 1","Time 1","Stimulus 2","Time 2","Sex")],collapse="_"))
  }
  else {
    if ( b["Stimulus 1"] %in% c("salmonella.flgB","salmonella.prgH","salmonella.wt","salmonella.flgBpEm62") ){
      insert <- ""
    } else { insert <- "in-vitro"}

    if ( b["Time 2"] != "" ) {
      return(paste(c(b2["Strain"],paste(b2["Stimulus 1"],"+",b2["Stimulus 2"],sep=""),insert,b2[c("Time 1")],insert,b2[c("Time 2","Cell Type")],"Mouse"),collapse="_"))
    } else {
      return(paste(c(b2[c("Strain","Stimulus 1")],insert,b2[c("Time 1","Stimulus 2","Time 2","Cell Type")],"Mouse"),collapse="_"))
    }
  }
}
  
dashify <- function( inString ){
  vecform <- strsplit(inString,split = " ")[[1]]
  return(paste(vecform,collapse="-"))
}

## convert whitespaces to %20
whiteSpaceURLform <- function ( instring ){
 return(gsub(" ","%20",instring))
}
 
arrayStats <- function ( arrayList, startDate, endDate ){

  ##startDate <- "1/1/06"
  ##endDate <- "12/31/06"
  n.arrays <- length(arrayList )
  
  chips <- character(length=n.arrays)
  dates <- vector(length=n.arrays)
  
  for ( n in 1:n.arrays ){
    
    uri.long <- listi[[n]]$uri

    uri.long <- gsub("\r\n98\r\n","",uri.long)
    uri.long <- gsub("\r\n171b\r\n","",uri.long)
    uri.long <- gsub("\r\n[[:alnum:]]{1,50}\r\n","",uri.long)
    
    toks <- strsplit(uri.long,"/")[[1]]
    ntoks <- length(toks)
    uri <- toks[ntoks]

    aobj <- getNodeObject(uri)

    chip.long <- aobj$chip
    toks <- strsplit(chip.long,"/")[[1]]
    ntoks <- length(toks)
    chip <- toks[ntoks]
    chips[n] <- chip

    date <- aobj[["submission_date"]]
    date <- dateConvert1(date)
    wir <- whichInRange(date,startDate,endDate)
    dates[n] <- length(wir)

    
  }    

  return(chips[which(dates==1)])
  ##return(cbind(dates,chips))
}
 
getURIsByDate <- function ( arrayList, startDate, endDate ){

  ##startDate <- "1/1/06"
  ##endDate <- "12/31/06"
  n.arrays <- length(arrayList )
  
  chips <- character(length=n.arrays)
  dates <- vector(length=n.arrays)
  uris <- character(length=n.arrays)
  
  for ( n in 1:n.arrays ){
    
    uri.long <- arrayList[[n]]$uri

    uri.long <- gsub("\r\n98\r\n","",uri.long)
    uri.long <- gsub("\r\n171b\r\n","",uri.long)
    uri.long <- gsub("\r\n93\r\n","",uri.long)

    toks <- strsplit(uri.long,"/")[[1]]
    ntoks <- length(toks)
    uri <- toks[ntoks]
    uris[n] <- uri

    aobj <- getNodeObject(uri)

    chip.long <- aobj$chip
    toks <- strsplit(chip.long,"/")[[1]]
    ntoks <- length(toks)
    chip <- toks[ntoks]
    chips[n] <- chip

    date <- aobj[["submission_date"]]
    date <- dateConvert1(date)
    wir <- whichInRange(date,startDate,endDate)
    dates[n] <- length(wir)

    
  }    

  return(uris[which(dates==1)])
}


filterArrayListByDate <- function ( arrayList, startDate, endDate ){
  ##startDate <- "1/1/06"
  ##endDate <- "12/31/06"
  n.arrays <- length(arrayList)
  outList <- list()
  for ( n in 1:n.arrays ){
    arrayObj <- arrayList[[n]]
    uri.long <- arrayList[[n]]$uri
    uri.long <- gsub("\r\n98\r\n","",uri.long)
    uri.long <- gsub("\r\n171b\r\n","",uri.long)
    uri.long <- gsub("\r\n93\r\n","",uri.long)
    toks <- strsplit(uri.long,"/")[[1]]
    ntoks <- length(toks)
    uri <- toks[ntoks]
    aobj <- getNodeObject(uri)
    date <- aobj[["submission_date"]]
    date <- dateConvert1(date)
    wir <- whichInRange(date,startDate,endDate)
    if ( length(wir) > 1 ){ cat("What's up with that?\n") }
    if ( length(wir) == 1 ){ outList[[length(outList)+1]] <- arrayObj }
  }    
  return(outList)
}


##
## Finding dates in an interval
## library(chron)
## which(!is.na(cut(dates("07/01/92") + 0:14, breaks <- dates(c( "07/03/92", "07/07/92" )))))

## Convert "Wed Dec 20 08:00:00 GMT 2006"
## To "Dec 20 2006"
dateConvert1 <- function ( dstring ) {
  toks <- strsplit(dstring," ")[[1]]
  m <- toks[2]
  d <- toks[3]
  y <- toks[6]
  return(paste(c(m,d,y),collapse= " "))
}

## convert "Dec 20 2006" to corresponding chron object
## also works where instring is a vector
dateConvert2 <- function ( instring ){
  return( dates(instring ,format="month day year"))
}

## input: date vector in form of "Dec 20 2006"
## startDate: earliest date, of form "1/1/06"
## endDate: end date, of form "12/1/06"
## output: indices of dates within range 

whichInRange  <- function( dates, startDate, endDate ){
  dates <- dateConvert2(dates)
  inds <- which(!is.na(cut(dates, breaks <- dates(c(startDate,endDate )))))
  return(inds)
}



## logical TRUE if query (list) is a subset of target (list)
inMeta <- function(query,target){
  naymes <- names(query)
  return( identical( query, target[naymes] ))
}


## logical TRUE if query (any object) is found in targetlist (list of lists)
inList <- function(query,targetList){
  logvec <- unlist(lapply(targetList,identical,query))
  return( TRUE %in% logvec )
}


## INDICES of query (any object), if found in targetlist (list of lists)
## Allows for a partial match, meaning that the
## matched object in the targetlist can have additional name-value pairs

inListSoft <- function(query,targetList){
  naymes <- names(query)
  indvec <- numeric()
  counter <- 0 
  for ( targ in targetList ){
    counter <- counter + 1
    if ( inMeta(query,targ) ){
      indvec <- c(indvec,counter)
    }
  }
  return(indvec)
}

   
## Retain only single stim
removeStim2 <- function( inObj ){
  keep.these <- which(unlist(lapply(lapply(inObj,"[[","Stimulus 2"),is.null)))
  return( inObj[keep.these] )
}

timesAsString <- function ( obj ){
  tv <- table(unlist(lapply(obj,"[[","Time 1")))
  return(paste(paste(names(tv),"(",as.vector(tv),")",sep=""),collapse=";"))
}

# integer 20 to string 0020, e.g.
timeAsPaddedString <- function (t){
  pre <- paste(rep("0",4-nchar(as.character(t))),collapse="")
  term.t <- paste(c(pre,as.character(t)),collapse="")
  term.t
}

timesAsVector <- function ( obj ){
  tv <- table(unlist(lapply(obj,"[[","Time 1")))
  return(tv)
}

timesByPlatformAsString <- function ( obj ){
  chips <- gsub(Estring,"E",gsub(Tstring,"T",unlist(lapply(obj,"[[","chip"))))
  ecourses <- obj[which(chips=="E")]
  tcourses <- obj[which(chips=="T")]
  tstring <- ""
  if ( length(ecourses)>0 ){
    te <- timesAsString(ecourses)
    tstring <- paste(tstring,"E:",te," ",sep="")
  }
  if ( length(tcourses)>0 ){
    tt <- timesAsString(tcourses)
    tstring <- paste(tstring,"T:",tt," ",sep="")
  }
  return(tstring)
}


timesByPlatformAsVectors <- function ( obj ){
  chips <- gsub(Estring,"E",gsub(Tstring,"T",unlist(lapply(obj,"[[","chip"))))
  ecourses <- obj[which(chips=="E")]
  tcourses <- obj[which(chips=="T")]
  returnList <- list()
  if ( length(ecourses)>0 ){
    tv <- timesAsVector(ecourses)
    returnList[["E"]] <- tv
  }
  if ( length(tcourses)>0 ){
    tv <- timesAsVector(tcourses)
    returnList[["T"]] <- tv
  }
  return(returnList)
}


platformCountAsString <- function ( obj ){
  tv <- table(unlist(lapply(obj,"[[","chip")))
  outString <- paste(paste(names(tv),"(",as.vector(tv),")",sep=""),collapse=";")
  outString <- gsub(Estring,"E",outString)
  outString <- gsub(Estring.broken,"E",outString)
  outString <- gsub(Tstring,"T",outString)
  return(outString)
}
