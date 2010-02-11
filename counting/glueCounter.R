
library(rjson)
library(httpRequest)
library(chron)

util.dir <- file.path(Sys.getenv("AA"),"utils")
source(file.path(Sys.getenv("HOME"),"bin/R/MatrixPrintFormat.R"))
source(paste(util.dir,"httpget.R",sep="/"))
source(paste(util.dir,"utilitiesMetaData.R",sep="/"))

path = "/sampleData/microarray/samples/Aderem%20Lab"

##This report covers
##the period from August 30, 2008 - February 28, 2009.
## format m/d/y
startDate <- "8/30/08"
endDate <- "2/28/09"


##the period from March 1, 2008 - August 29, 2009 
## format m/d/y
startDate <- "3/1/09"
endDate <- "8/29/09"


##the period from August 29, 2009 - February 28, 2010 
## format m/d/y
startDate <- "8/29/09"
endDate <- "2/28/10"

chips <- c("Mouse Exon","Mouse 430 2.0","Mouse Promoter 1.0R","Isbimm","Human U133 Plus 2.0","Human U133A","Human Exon","Rhesus Macaque")

glueprojects <- c(
##"Aderem Custom ChIP",
"TLR Specificity",
"Litvak Samples",
"Exon Arrays",
"LPS Time Course",
"PIC Time Course",
##"Aderem Promoter Arrays",
"Influenza",
"Liz Gold Arrays",
"Dendritic Cell Arrays",
"Leishmania Infected Cells",
"Viral Infections",
"Glue Genomics",
"GLUE grant"  )

collect <- numeric()
keepers <- character()

for ( proj in glueprojects ){
##for ( i in 1:2 ){

  cat(proj,"\n")
  proj.p20 <- whiteSpaceURLform(proj)
  proj.path <- paste(path,proj.p20,sep="/")
  response = getToHost("innate-immunity.systemsbiology.net",paste("/addama-rest/primary-repo?PATH=",proj.path,sep=""),"",port=80)
  response = sub("1ff8","",response)
  response = sub("\r\n1f29\r\n","",response)
  response = gsub("\r\n1ff8\r\n","",response)
  response = gsub("\r\n93\r\n","",response)
  response = gsub("\r\n[[:alnum:]]{1,50}\r\n","",response)
  ## chop HTTP header stuff in response
  resp = substr(response,regexpr("\\{",response)[1],nchar(response))
  resp = sub("}\r.*","}",resp)
  object = fromJSON(resp)

  listi <- object$children

  if ( length(listi ) >= 1){

    oo <- table(arrayStats(listi,startDate,endDate))

    if ( length(oo)>=1 & sum(oo)>0 ) {
    
      blank <- numeric(length=length(chips))
      names(blank) <- chips
      blank[names(oo)] <- as.numeric(oo)   
      collect <- rbind(collect,blank)
      keepers <- c(keepers,proj)
    }
  }
  
  ##samplecounts[proj] <- length(object[["children"]])
  
}

rownames(collect) <- keepers

write.table(matrixPrintFormat(collect),file="1MarchTo26Aug2009.tsv",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)



