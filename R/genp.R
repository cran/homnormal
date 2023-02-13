#' Generalized p value Test for Homogeniety
#'
#' Tests the homogeniety of variances for more than two normal groups
#' using generalized p value test.
#' @export
#' @param x1 a numeric matrix containing the values of groups.
#'
#' @param x2  numeric matrix containing the values of group numbers.
#'
#' @param m  number of resampling.
#'
#' @param alfa significance level of the test. Default number is 0.05.
#'
#' @param table a logical variable that indicates table will appear or not. Default is TRUE.
#'
#' @param graph box plot of groups of raw or centered data.
#'
#'@return if table is TRUE, then it gives a detailed table,
#'else it gives a vector of r value(r=1 when null hypothesis was rejected and r=0 when null hypothesis was accepted)
#' p-value and test statistic value.
#'
#' @examples
#'   data(FH_data)
#'    x1=FH_data$SurvivalTime
#'    x2=FH_data$HospitalNo
#'    genp(x1,x2)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'    genp(x1,x2,alfa=0.10)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'    genp(x1,x2,alfa=0.10,m=5000)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'     genp(x1,x2,alfa=0.10,table=FALSE)
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     genp(x1,x2,alfa=0.10,table=FALSE,graph="raw")
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     genp(x1,x2,alfa=0.10,table=FALSE,graph="none")
#' # ---THIS VERSION IS ESPECIALLY USEFUL FOR COMPARISON STUDIES BY SIMULATION---
#' #    #first value of the vector is r value(r=1 when rejected and r=0 when accepted null hypothesis)
#' # second value of the vector is the p-value and third value is the tests statistic value
#' @seealso \code{\link[homnormal]{Brown_Forsythe}}, \code{\link[homnormal]{Cat_GG}}, \code{\link[homnormal]{Cat_LR}}, \code{\link[homnormal]{bdai}}, \code{\link[homnormal]{slrt}}, \code{\link[homnormal]{levene}}
#' @references Liu, X., & Xu, X. (2010). A new generalized p-value approach for testing the homogeneity of variances. Statistics & probability letters, 80(19-20), 1486-1491.
genp=function(x1,x2,alfa=0.05,m=2000,table=TRUE,graph="none") {
  mut=c()
  s2=c()
  n=c()
  for (i in factor(x2)) {
    n[as.numeric(i)]=(sum(x2==i))
  }
  k=length(n)
  for (i in 1:k) {
    mut[i]=mean(x1[x2==i])
    s2[i]=var(x1[x2==i])
  }
  E=c()
  VR=c()
  R=matrix(, nrow=k, ncol=m)
  H=matrix(, nrow=k-1, ncol=k)
  for (i in 1:k) {
    f<- function(x) log(x)*x^((n[i]-3)/2)*exp(-x/2)/2^((n[i]-1)/2)/gamma((n[i]-1)/2)
    g<-  function(x) log(x)^2*x^((n[i]-3)/2)*exp(-x/2)/2^((n[i]-1)/2)/gamma((n[i]-1)/2)
    E[i]=integrate(f,0,Inf)$value
    VR[i]=integrate(g,0,Inf)$value-E[i]^2
    }
  H[1:(k-1),1:(k-1)]=diag(k-1)
  H[1:(k-1),k]=-1
  for (i in 1:k) {
    R[i,]=log(n[i]*s2[i])-log(rchisq(m,n[i]-1))
  }
  ER=log(n*s2)-E
  muHR=H %*% matrix(ER,nrow=k)
  HR=H %*% R
  SHR=H %*% diag(VR) %*% t(H)
  TGP=diag(t(HR- matrix(rep(muHR,m),ncol=m))%*%solve(SHR)%*%(HR- matrix(rep(muHR,m),ncol=m)))
  T=t(muHR)%*%solve(SHR)%*%muHR
  p=mean(TGP>drop(T))
  r=p<alfa
  if(graph=="raw") {
    boxplot(x1~x2,xlab = "", ylab = "")
  }
  else if(graph=="centered") {
    boxplot(x1-rep(mut,ne)~x2,xlab = "", ylab = "")
  }
  else if (graph=="none")
  {}
  if (table==FALSE) {
    return(c(r,p,T))
  }
  else if (table==TRUE)
  {
    C1=1:length(n)
    C2=c(n)
    C3=c(mut)
    C4=c(s2)
    C5=c(rep(NA,ifelse(k%%2==0,ceiling(k/2),floor(k/2))),T,rep(NA,floor(k/2-1)))
    C6=c(rep(NA,ifelse(k%%2==0,ceiling(k/2),floor(k/2))),p,rep(NA,floor(k/2-1)))
    HT<-hux(Group_No=C1,Sample_Size=C2,Sample_Mean=C3,Sample_Var=C4,Test_Stat=C5,p_value=C6)
    align(HT)[,1:5]          <- 'center'
    number_format(HT)[,3:6]      <- 2
    number_format(HT)[,2]      <- 0
    
    
    print_md(HT, header = TRUE)
    return(HT)
  }
}