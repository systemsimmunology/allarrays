
##Requires
##library(rjson)
##source("httpget.R")

## Parse Sample group object retrieved from JCR
## Input: Sample Group object from JCR
## Output: The actual sample groups, which
## are child nodes of the Input
parseSampleGroup <- function ( sampleGroup ){
  sgkids <- sampleGroup[["children"]]
  group.names <- unlist(lapply(sgkids,"[[","name"))
  sgkids.uri <- unlist(lapply(sgkids,"[[","uri"))
  sglist <- list()
  for ( i in 1:length(group.names)){
    group <- group.names[i]
    uri.long <- sgkids.uri[i]
    toks <- strsplit(uri.long,"/")[[1]]
    ntoks <- length(toks)
    uri <- toks[ntoks]
    sglist[[group]] <- getNodeObject(uri)
  }
  return(sglist)
}

## Remove Stim2 Sample Groups
removeStim2SG <- function ( sglist.in ){
  ## These include the string "NULL"
  sg.stim2 <- as.character(lapply(sglist.in,"[[","Stimulus 2"))
  have.stim2 <- which(sg.stim2 != "NULL")
  ## Let's remove these
  sglist <- sglist.in[-have.stim2 ]
  return(sglist)
}  


## Remove time1 != time2  Sample Groups
removeUnequalTime <- function ( sglist.in ){
  ## These include the string "NULL"
  sg.time1 <- as.character(lapply(sglist.in,"[[","Time 1"))
  sg.time2 <- as.character(lapply(sglist.in,"[[","Time 2"))
  have.time1 <- (sg.time1 != "NULL")
  have.time2 <- (sg.time2 != "NULL")
  unequal <- sg.time1 != sg.time2
  chuckthese <- which(have.time1 & have.time2 & unequal)
  if ( length(chuckthese) > 0 ){
    sglist <- sglist.in[-chuckthese ]
  } else {
    sglist <- sglist.in
  }
  return(sglist)
}  


## This is temporary code
## Creates Male and Female sample group
## replacing Sex=NULL cases
repairSGgender <- function ( sglistIn ){

  sglistOut <- list()
  
  for ( i in 1:length(sglistIn) ) {
    sg <- sglistIn[[i]]
    sgn <- names(sglistIn[i])
    genderFine <- TRUE
    ## If Sex field is missing that's bad
    if ( !("Sex" %in% names(sg) )){
      genderFine <- FALSE
    }
    ## If it's there but is NULL, that is bad, too
    if ( "Sex" %in% names(sg) ){
      if ( is.null(sg[["Sex"]]) ){
        genderFine <- FALSE
      }
    }
    if ( genderFine ){
      sglistOut[[sg$name]] <- sg
    }
    if ( !genderFine ){
      sgM <- sg
      sgM[["Sex"]] <- "Male"
      sgnM <- paste(sgn,"_Male",sep="")
      sgM[["name"]] <- sgnM
      sgF <- sg
      sgF[["Sex"]] <- "Female"
      sgnF <- paste(sgn,"_Female",sep="")
      sgF[["name"]] <- sgnF
      ## out with the old, in with the (2) new
      ##sglist <- sglist[-i]
      ##sglist <- sglist[setdiff(names(sglist),sg$name)]
      sglistOut[[sgnM]] <- sgM
      sglistOut[[sgnF]] <- sgF
    }
  }
  return(sglistOut)
}
      

## For a Sample Group List Object
## Determine unique combinations of
## *C*ell Type, *S*train, *S*timulus 1 (and Sex)
## Returns only CSSs for which Stimulus 1 is not NULL, Unstimulated, etc.

getUniqueCSSs <- function ( sglist ){

  ## Build up unique css list
  labels <- c("Cell Type","Strain","Stimulus 1","Stimulus 2","Sex")
  
  uniqueCSSs <- list()
  for ( sg in sglist ) {
    query <- list()
    query[labels] <- sg[labels]

    if ( !is.null(query[["Stimulus 1"]]) ){## if not null, keep going
      if ( (query[["Stimulus 1"]] != "Unstimulated") & (query[["Stimulus 1"]] != "Mock infected") ){ ## if true stim, keep going
        ## For Stim 2, we don't require that it be 'non-zero'
        ## But we do want to use "", for consistency, for NULL, Unstim
        if ( is.null(query[["Stimulus 2"]]) ) {
          query[["Stimulus 2"]] = ""
        } else { 
          if ( (query[["Stimulus 2"]] != "Unstimulated") & (query[["Stimulus 2"]] != "Mock infected") ){ ## if true stim, keep going
            query[["Stimulus 2"]] = query[["Stimulus 2"]]
          } else {
            query[["Stimulus 2"]] =""
          }
        }
        stringform <- paste(query,collapse="_")
        if ( !inList(query,uniqueCSSs) ){
          uniqueCSSs[[stringform]] <- query
          uniqueCSSs[[stringform]][["name"]] <- stringform
        }
      }
    }
  }
  return(uniqueCSSs)
}

## For a Sample Group List
## Determine unique combinations of
## *C*ell Type, *S*train, (and Sex)
## for cases in which Stimulus 1 *is* not NULL, Unstimulated, etc.
## Stimulus 1 is set to "" in the output list

getUnstimmedCSSs <- function ( sglist ){
  nonstim.for.output <- ""
  ## Build up unique css list
  labels <- c("Cell Type","Strain","Stimulus 1","Sex")
  uniqueCSSs <- list()
  for ( sg in sglist ) {
    query <- list()
    query[labels] <- sg[labels]
    if ( is.null(query[["Stimulus 1"]]) ){ query[["Stimulus 1"]] <- nonstim.for.output }
    if ( query[["Stimulus 1"]]=="Unstimulated" ){ query[["Stimulus 1"]] <- nonstim.for.output }
    if ( query[["Stimulus 1"]]=="Mock infected" ){ query[["Stimulus 1"]] <- nonstim.for.output }
    if ( query[["Stimulus 1"]]=="" ){
      stringform <- paste(query,collapse="_")
        if ( !inList(query,uniqueCSSs) ){
          uniqueCSSs[[stringform]] <- query
          uniqueCSSs[[stringform]][["name"]] <- stringform
        }
    }
  }
  return(uniqueCSSs)
}

##
## Arrange Sample Groups into Timecourses for fixed
## Cell Type, Strain Stimululus
##
##
findTimeCourse <- function( CSSs, sglist ){

  CSSout <- list()
 
  for ( umname in names(CSSs) ){
    um <- CSSs[[umname]][c("Cell Type","Strain","Stimulus 1","Stimulus 2","Sex")]
    CSSout[[umname]] <- CSSs[[umname]] ## copy original metadata
    
    for ( sgname in names(sglist) ){

      sg <- sglist[[sgname]]
      ## Have to transform unstim Stim 2 
      if ( is.null(sg[["Stimulus 2"]]) ) {
        sg[["Stimulus 2"]] = ""
      } else {
        if ( sg[["Stimulus 2"]] == "Unstimulated" | sg[["Stimulus 2"]] == "Mock infected") {
          sg[["Stimulus 2"]] = ""
        }
      }      
      if ( inMeta(um,sg) ){
        CSSout[[umname]][["Sample Group"]] <- c(CSSout[[umname]][["Sample Group"]],sgname)
        CSSout[[umname]][["Time 1"]] <- c(CSSout[[umname]][["Time 1"]],as.numeric(sg[["Time 1"]]) )
        CSSout[[umname]][["Replicate Count"]] <- c(CSSout[[umname]][["Replicate Count"]],length(sg[["sampleUUIDs"]]))
      }
    }
  }
  ## Time sort
  for ( umname in names(CSSs) ){
    if ( length(CSSout[[umname]][["Time 1"]]) > 1 ){
      time.order <- sort(CSSout[[umname]][["Time 1"]],index.return=TRUE)$ix
      CSSout[[umname]][["Sample Group"]] <- CSSout[[umname]][["Sample Group"]][time.order]
      CSSout[[umname]][["Time 1"]] <- CSSout[[umname]][["Time 1"]][time.order]
      CSSout[[umname]][["Replicate Count"]]  <- CSSout[[umname]][["Replicate Count"]][time.order]
    }
  }

  return(CSSout)
}

##
## Find time zero sample groups
##
## Code divides into two according to whether there
## are special cases of hand-curated choices
findTimeZero <- function( CSSs, sglist, zfile = NULL ) {
  
  ## Attempt to find 0 time point for each time course
  unstim <- names(which(unlist(lapply(CSSs,"[[","Stimulus 1"))=="Unstimulated"))
  mock <- names(which(unlist(lapply(CSSs,"[[","Stimulus 1"))=="Mock infected"))
  nulls <- names(which(unlist(lapply(lapply(CSSs,"[[","Stimulus 1"),is.null))))
  stimmed <- setdiff(names(CSSs),union(unstim,union(nulls,mock)))
  ## These are stimmed CSSs

  sg.unstim <- names(which(unlist(lapply(sglist,"[[","Stimulus 1"))=="Unstimulated"))
  sg.nulls <- names(which(unlist(lapply(lapply(sglist,"[[","Stimulus 1"),is.null))))
  sg.mock <- names(which(unlist(lapply(sglist,"[[","Stimulus 1"))=="Mock infected"))
  possiblezeros <- union(sg.unstim,union(sg.nulls,sg.mock))
  ## The are candidate zero time cases among the sample groups
  
  if ( !is.null(zfile) ) {
    zd <- read.table(zfile,header=1,as.is=T,sep="\t")
    zerotime <- list()
    for ( umname in stimmed ){
      um <- CSSs[[umname]]
      query <- list()
      query[["Strain"]]  <- um[["Strain"]]
      query[["Cell Type"]] <- um[["Cell Type"]]
      ind <- which((zd[["Strain"]]==query[["Strain"]])&(zd[["Cell.Type"]]==query[["Cell Type"]]))
      if ( length(ind) >= 2 ){cat("Problem:More than one zero time match/n")}
      if ( length(ind) == 1 ){
        zerotime[[umname]] <- zd[ind,"TimeZeroChoice"]
      }
      if ( length(ind) == 0 ){
        for ( sg in sglist[possiblezeros] ){
          if ( inMeta ( query, sg ) ){
            zerotime[[umname]] <- sg$name
          }
        }
      }
    }
  } else { ## there is no zerofile
    zerotime <- list()
    for ( umname in stimmed ){
      um <- CSSs[[umname]]
      query <- list()
      query[["Strain"]]  <- um[["Strain"]]
      query[["Cell Type"]] <- um[["Cell Type"]]
      for ( sg in sglist[union(sg.unstim,union(sg.nulls,sg.mock))] ){
        if ( inMeta ( query, sg ) ){
          zerotime[[umname]] <- sg$name
        }
      }
    }
  } ## end ops for no zero file case

  return( zerotime )
  
}
  
##
## Three possible groups for Bl6 unstim
## "Bl6_Unstimulated_in-vitro_0000___BMDM_Mouse" 6 samples
## "Bl6______BMDM_Mouse" 1 sample
## "Bl6_Unstimulated_in-vitro_0240___BMDM_Mouse" 3 samps

## Let's use the 3 sample group for now "Bl6_Unstimulated_in-vitro_0240___BMDM_Mouse"

## Like wise for Strain "Dicer fl/fl LyszCre+/-"
## there are two "choices
## "Dicer-fl-fl-LyszCre+--______BMDM_Mouse"
## "Dicer-fl-fl-LyszCre+--_Unstimulated_in-vitro_0000___BMDM_Mouse"

## Hard to tell which is better !

## m["submission_date" ,] "2009-02-23T00:00:00.000-08:00" "2008-11-12T00:00:00.000-08:00"
## m["sbeams_user" ,]  "rsuen"    "jsissons"

## Let's take the newer one:
## "Dicer-fl-fl-LyszCre+--______BMDM_Mouse"

includeZeroTime <- function ( CSS.timecourse ) {
  
  for ( csst.name in names(CSS.timecourse) ){ ## these are just stimmed at this point
    csst <- CSS.timecourse[[csst.name]]
    zt <- zeroTimes[[csst.name]] ## sample group name
    CSS.timecourse[[csst.name]][["Time 1"]] <- c(0, CSS.timecourse[[csst.name]][["Time 1"]] )
    CSS.timecourse[[csst.name]][["Sample Group"]] <- c(zt,CSS.timecourse[[csst.name]][["Sample Group"]])
    ##  samps[[umname]] <- c(zt,samps[[umname]])
    ##  times[[umname]] <- c(0,times[[umname]])
    ##  nreps[[umname]] <- c(length(sglist[[zt]][["sampleUUIDs"]]),nreps[[umname]])
  }

  return(CSS.timecourse)
}

##
## Queries list for the existence of a defined element
## return empty string if not found
## 
nullToBlank <- function (inlabel,inList){
  if (!(inlabel %in% names(inList))  ){
    outVal <- ""
  } else {
    if ( is.null(inlabel)) {
      outVal <- ""
    } else {
      outVal <- inList[[inlabel]]
    }
  }
  return(outVal)
}


##
## Create mapping of sample group names to expression matrix columns
## 
DMcolsFromSG <- function ( sglist, columns, exon.form=TRUE  ) {

  dmcol.from.sg <- vector()
  for ( sg in sglist ){ 
    sg[["Stimulus 1"]] <- nullToBlank("Stimulus 1",sg)
    sg[["Time 1"]] <- nullToBlank("Time 1",sg)
    sg[["Stimulus 2"]] <- nullToBlank("Stimulus 2",sg)
    sg[["Time 2"]] <- nullToBlank("Time 2",sg)
    dm.stringform <- metaToString(meta(strain=sg[["Strain"]],stimulus1=sg[["Stimulus 1"]],stimulus2=sg[["Stimulus 2"]],time1=sg[["Time 1"]],time2=sg[["Time 2"]],cell.type=sg[["Cell Type"]],sex=sg[["Sex"]]),exon.form)
    is.there <- dm.stringform %in% columns  
    ## Temporary try to replace. Should work with both male and female
    if ( !is.there  & is.null(sg[["Sex"]]) ){
      dm.stringform <- metaToString(meta(strain=sg[["Strain"]],stimulus1=sg[["Stimulus 1"]],stimulus2=sg[["Stimulus 2"]],time1=sg[["Time 1"]],time2=sg[["Time 2"]],cell.type=sg[["Cell Type"]],sex="Female"),exon.form)
      is.there <- dm.stringform %in% columns
    }
    if ( !is.there & is.null(sg[["Stimulus 1"]]) ){
      dm.stringform <- metaToString(meta(strain=sg[["Strain"]],stimulus1="",stimulus2=sg[["Stimulus 2"]],time1="",time2=sg[["Time 2"]],cell.type=sg[["Cell Type"]],sex=sg[["Sex"]]),exon.form)
      is.there <- dm.stringform %in% columns
    }
    ## Hmmm...
    if ( dm.stringform == "Myd88-null_D4_in-vitro_0120___BMDM_Mouse" ){
      dm.stringform <- "Myd88-null_D4_in-vitro_0120_in-vitro__BMDM_Mouse"
      is.there <- dm.stringform %in% columns
    }
    if ( !is.there ){
      cat(sg$name, " ", dm.stringform,"\n")
    }
    if ( is.there ) {
      dmcol.from.sg[sg$name] <- dm.stringform
    }
  }
  ## Exceptions
##  dmcol.from.sg["IRF1-null______BMDM_Mouse"] <-  "BMDM_IRF1-null_____Female"
##  dmcol.from.sg["MyD88-Trif-null_Unstimulated_in-vitro_0000___Dendritic-Cell_Mouse"] <- "Dendritic-Cell_MyD88-Trif-null_Unstimulated_0240___Male$"
  
  return(dmcol.from.sg)

}
  
