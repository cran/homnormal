#' Multiple Normal Distribution Data Generation
#'
#'  This function generates data from multiple normal distribution.
#' @export
#'@param n Sample sizes of each group. n=c(n1,n2,...nk); for example: n=c(3, 4, 5).
#'
#'@param mu Mean  of each group.mu=c(mu1,mu2,...muk); for example: mu=c(3, 4, 5).
#'
#'@param sigma Standard deviation of each group.sigma=c(sigma1,sigma2,...sigmak);
#'for example: sigma=c(1, 2, 3).
#'@param tn Trial number for all groups. Default of the parameter is 1.
#' This parameter for use more than 1, is especially useful for resampling
#' such as Monte Carlo, Parametric Bootstrap.
#' 
#' @return a data matrix with size (n1,n2,...nk) with group number 1,2,...k at first row and random number with
#' mnean mu=(mu1,mu2,...muk) and standard deviation  sigma=(sigma1,sigma2,...sigmak)
#' 
#' @examples
#' 
#' n=c(3, 4, 5)
#' mu=c(3, 4, 5)
#' sigma=c(3, 4, 5)
#' F=datagen(n,mu,sigma);muh=F[1];S2h=F[2];x=F[3]
#' muh
#' S2h
#' x
#'
#' # Following example especially useful for simulation based tecnhiques
#' # such as Monte Carlo, Parametric Bootstrap and comparison studies
#' # by using simulation.
#'
#'  Fm=datagen(c(3, 4, 5),c(3, 4, 5),c(3, 4, 5),10);muhm=Fm[1];S2hm=Fm[2];xm=Fm[3]
#' muhm
#' S2hm
#' xm

datagen=function(n,mu,sigma,tn=1) {
  e=c(0, cumsum(n))
  k=length(n)
  a=matrix(,nrow=k, ncol=2)
  for (i in 1:k){
    a[i,1]= e[i]+1
    a[i,2]=e[i+1]
  }
  x=matrix(,nrow=tn+1, ncol=sum(n))
  mut=matrix(,nrow=tn, ncol=k)
  mu2t=matrix(,nrow=tn, ncol=k)
  s2t=matrix(,nrow=tn, ncol=k)
  for (i in 1:k){
    x[1,a[i,1]:a[i,2]]=i
    x[2:(tn+1),a[i,1]:a[i,2]]=matrix( rnorm(n[i]*tn,mu[i],sigma[i]), tn,n[i])
    if (tn>1) {
      mut[1:tn,i]=rowMeans(x[2:(tn+1),a[i,1]:a[i,2]])
      mu2t[1:tn,i]=rowMeans(x[2:(tn+1),a[i,1]:a[i,2]]^2)
      s2t[1:tn,i]=(mu2t[1:tn,i]-mut[1:tn,i]^2)*n[i]/(n[i]-1)
    } else {
      mut[1:tn,i]=mean(x[2:(tn+1),a[i,1]:a[i,2]])
      mu2t[1:tn,i]=mean(x[2:(tn+1),a[i,1]:a[i,2]]^2)
      s2t[1:tn,i]=(mu2t[1:tn,i]-mut[1:tn,i]^2)*n[i]/(n[i]-1)
    }
  }
  return(list(mut,s2t,x))
}