

rawQCplot <- function( data, bottom_margin=27,o=NULL )
{
  if(is.null(o)) o <- 1:nrow(data$endog)
  n <- length(o)
  par(oma=c(bottom_margin,5,3,3))
  par(mfrow=c(3,1))

  p <- ceiling(log10(max(data$endog[o,])))
  q <- ceiling(log10(quantile(data$endog[o,],.05)))
  tt <- c(0,c(sapply(q:(p+1),function(r) c(.3,1)*10^r)))
  dr <- range(data$endog)
  
  par(cex.main=2)

  plot.new()
  plot.window(xlim=c(0,n+1),ylim=dr+.5,xaxs="i",log="y")
  boxplot.matrix( data$endog[o,]+.5, use.cols=F,add=T,
    xaxt="n", yaxt="n",frame.plot=F, main="endogenous" )
  axis(side=2,at=0.5+tt,labels=tt,las=2,cex.axis=1.5)
  mtext(side=2,at=sum(log10(dr))/2,text="raw counts",line=5)
  title(main="endogenous")

  plot.new()
  plot.window(xlim=c(0,n+1),ylim=dr+.5, 
    xaxs="i",
    ylab="raw counts",xlab="",main="+ spike-ins",log="y")
  boxplot.matrix( data$endog[o,]+.5, use.cols=F,add=T,xaxt="n", yaxt="n",frame.plot=F )
  axis(side=2,at=0.5+tt,labels=tt,las=2,cex.axis=1.5)
  mtext(side=2,at=sum(log10(dr))/2,text="raw counts",line=5)
  title(main="+ spike-ins")
  for(i in 1:ncol(data$pos))
    lines(1:n, .5 + data$pos[o,i],col=2,lwd=2)
  for(i in 1:ncol(data$neg))
    lines(1:n, .5 + data$neg[o,i],col=4,lwd=2)

  plot.new()
  plot.window(xlim=c(0,n+1),ylim=dr+.5,xaxs="i",
    ylab="raw counts",xlab="",main="+ housekeeping",log="y")
  boxplot.matrix( data$endog[o,]+.5, use.cols=F,add=T,xaxt="n", yaxt="n",frame.plot=F )
  axis(side=2,at=0.5+tt,labels=tt,las=2,cex.axis=1.5)
  title(main="+ housekeeping")
  mtext(side=2,at=sum(log10(dr))/2,text="raw counts",line=5)
  for(i in 1:ncol(data$hk))
    lines(1:n, .5 + data$hk[o,i],col="green4",lwd=2)

  mtext(side=1,at=1:n,text=rownames(data$endog)[o],line=3,las=2)
}
