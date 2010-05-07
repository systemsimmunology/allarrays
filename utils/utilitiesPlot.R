
### Plot functions
### ( Try not to write too many of these  ! )

## tmax is the maximum time to plot, in hours

plotCSS <- function( eid, css.timecourse, data.matrix=dm, ymax=NULL,main=NULL,logscale=FALSE, tmax=350,inset.text=NULL ){

  csst.label <- css.timecourse$name
  
  x <- css.timecourse[["Time 1"]]/60. ## Time in hours
  nt <- length(x) ## number of time points
  maxind <- length(which(x <= tmax ))
  x <- x[1:maxind]
  
  ## columns of datamatrix
  cols <- css.timecourse[["DM Column"]][1:maxind]

  ## expression values 
  if (!logscale ){ 
    tc <- data.matrix[eid,cols]
  } else {
    tc <- log2(data.matrix[eid,cols])
  }

  xlab <- "Time [hr]"
  ylab <- "Intensity"
  if ( logscale ){ ylab <- "Log2(Intensity)" }
  if  ( is.null(main) ) {
    main <- paste(gene.symbol[eid],":",csst.label)
  }
  if ( is.null(ymax) ){
    ylim <- c(0.,max(tc))
  } else {
    ylim <- c(0.,ymax)
  }
  ##cat(length(x),length(tc),"\n")
  plot(tc,x=x,xlab=xlab,ylim=ylim,ylab=ylab,type='l',col='blue',main=main)
  points(tc,x=x,type='p',col='blue',pch=19)

  if ( !is.null(inset.text) ){
    text(0.1*max(x),0.9*ymax,inset.text,cex=2)
  }
  
}


## guess at an evenly spaced grid needed
## to accomodate inval subplots
evengrid <- function ( inval ){
  fsi <- floor(sqrt(inval))
  is.square <- fsi == sqrt(inval)
  splitter <- fsi*(fsi+1)
  if ( is.square ){
    return( c(fsi,fsi))
  }
  if ( inval <= splitter ){
    return( c(fsi,fsi+1))
  }
  else{
    return( c(fsi+1,fsi+1))
  }
}
  
## for a single ID and multiple css, do plots
## ( see gridPredPlotWrapperODE etc. in utilitiesODEsolve.R ) 
gridPlotCSS <- function ( eid, css.tcs, data.matrix=dm, nx=NULL, ny=NULL , ymax=NULL, labvec=NULL,main='', tmax=350 ){
  
  ## could do an input parameter check
  
  ## PNG - seem to work only one page at a time, and gives crummy resolution
  ## png(file,height=800,width=1100,quality=100)
  
  ## EPS - may work! EPS not rendered in Word ( mac )
  ## postscript(file,onefile = FALSE,paper='letter')

  ## PS
  ## file=inptfile 
  ##postscript(file,paper='letter')

  ##Can't get pdf to do 'landscape'
  ##pdf(file,paper='letter',height=8.5,width=11)

  if ( is.null(nx) & is.null(ny) ){
    mfrow <- evengrid(length(css.tcs))
  } else {
    mfrow=c(nx,ny)
  }
  
  par(mfrow=mfrow,mar=(c(3, 2, 2, 1)+0.1),mgp=c(1.75, 0.75, 0),oma=c(0,0,3,0))

  if ( is.null(ymax) ){
    ymax <- max(data.matrix[eid,unlist(lapply(css.tcs,"[[","DM Column"))])
    ##ymax <- max(data.matrix[eid,])
  } else {
    ymax <- ymax
  }

  for ( css in css.tcs ){
    plotCSS(eid,css,data.matrix=data.matrix,ymax=ymax, tmax=tmax, inset.text=labvec[css$name] )
  }
  mtext(main, adj=0.5, side=3, cex=1.5, outer=TRUE)
}


## colorPlot of metadata
metaColorMat <- function ( inMat, colormap ){
  par(mar=c(0,0,0,0), cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
  plot.new()
  nr <- nrow(inMat)
  nc <- ncol(inMat)
  for (j in 1:nr ) {
    for (i in 1:nc ) {
      label = inMat[j,i]
      kolor = colormap[[label]]
      rect(i/nc,j/nr,(i-1)/nc,(j-1)/nr,col=kolor)
      text((i-0.5)/nc,(j-0.5)/nr,label,col=ifelse(mean(col2rgb(kolor))<120,"white","black"))
    }
  }
}
