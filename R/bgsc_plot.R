bgsc_plot <- function(data,i)
{
if(is.character(i))
  i <- grep(i,colnames(data$z))[1]

gpch <- c(endog=1,hk=20,neg=20,pos=20)[data$gtag]
gcol <- c(endog="darkgray",hk="green4",neg=4,pos=2)[data$gtag]
gcex <- c(endog=1,hk=2,neg=2,pos=2)[data$gtag]
tt <- c(0,c(sapply(0:5,function(r) c(1,3)*10^r)))

plot( data$mz, data$z[,i],pch=gpch,col=gcol,cex=gcex,axes=F,frame=T,
  main=colnames(data$z)[i],
  xlab="reference counts",ylab="sample counts")
abline(0,1,lty=2)
axis(side=1,labels=tt,at=log2(.5+tt))
axis(side=2,labels=tt,at=log2(.5+tt))
curve( log2( 0.5 + data$normpar[1,i] + 2^(x + data$normpar[2,i])),xlim=range(data$mz),col="orange3",add=T,lwd=2)
abline(h=log2(.5+data$normpar[1,i]),col="orange3",lty=2)
}
