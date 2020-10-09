# background and scale normalization
#
bgsc_norm <- function(data)
{
  z <- log2( 0.5 + t( cbind(data$endog,data$hk,data$neg,data$pos)))
  gtag <- c( rep("endog",ncol(data$endog)),
             rep("hk",ncol(data$hk)),
             rep("neg",ncol(data$neg)),
             rep("pos",ncol(data$pos)))
  names(gtag) <- rownames(z)

  w <- ifelse( grepl("^(POS)|(NEG)_",rownames(z)),0,1)
  
  mz <- apply(z,1,mean)
  mz <- mz-min(mz)+0.5

  lsf <- function(p,i) {
    ifelse(p[1] < 0, Inf, sum( w*(z[,i]-log2(0.5+p[1] + 2^(mz+p[2])))^2))
  }

  normpar <- sapply( 1:ncol(z),
    function(i) { optim( c(3,0),lsf,NULL,i)$par })

  zn <- sapply(1:ncol(z),function(i)
    z[,i] - log2(0.5+normpar[1,i]+2^(normpar[2,i]+mz)) +mz)

  structure(list(gtag=gtag,z=z,mz=mz,normpar=normpar,zn=zn),class="nsx_norm")
}
