
### Plot functions
### ( Try not to write too many of these  ! )

## tmax is the maximum time to plot, in hours

plotCSS <- function( eid, css.timecourse, data.matrix=dm, ymax=NULL,main=NULL,logscale=FALSE, tmax=100 ){

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
}
 
## for a single ID and multiple css, do plots
## ( see gridPredPlotWrapperODE etc. in utilitiesODEsolve.R ) 
gridPlotCSS <- function ( eid, css.tcs, data.matrix=dm, nx=1, ny=1 , ymax=NULL, labvec=NULL,main='', tmax=100 ){

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
 
  par(mfrow=c(nx,ny),mar=(c(3, 2, 2, 1)+0.1),mgp=c(1.75, 0.75, 0),oma=c(0,0,3,0))

  if ( is.null(ymax) ){
    ymax <- max(data.matrix[eid,unlist(lapply(css.tcs,"[[","DM Column"))])
    ##ymax <- max(data.matrix[eid,])
  } else {
    ymax <- ymax
  }

  for ( css in css.tcs ){
    plotCSS(eid,css,data.matrix=data.matrix,main=labvec[css$name],ymax=ymax, tmax=tmax)
    ##plotCSS(eid,css,data.matrix=data.matrix,main=css$name,ylim=ylim)
  }

  mtext(main, adj=0.5, side=3, cex=1.5, outer=TRUE)
  
}

              
