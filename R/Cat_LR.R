#' Computational Approach Test for Homogeniety
#'
#' Tests the homogeniety of variances for more than two normal groups
#' using standartized likelihood ratio test.
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
#' 
#'@return if table is TRUE, then it gives a detailed table,
#'else it gives a vector of r value(r=1 when null hypothesis was rejected and r=0 when null hypothesis was accepted)
#' p-value and test statistic value.
#'
#'
#'
#' @examples
#'
#'     data(FH_data)
#'    x1=FH_data$SurvivalTime
#'    x2=FH_data$HospitalNo
#'    Cat_LR(x1,x2)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'    Cat_LR(x1,x2,alfa=0.10)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'     Cat_LR(x1,x2,alfa=0.10,table=FALSE)
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     Cat_LR(x1,x2,alfa=0.10,table=FALSE,graph="raw")
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     Cat_LR(x1,x2,alfa=0.10,m=5000,table=FALSE,graph="none")
#' # ---THIS VERSION IS ESPECIALLY USEFUL FOR COMPARISON STUDIES BY SIMULATION---
#' #    #first value of the vector is r value(r=1 when rejected and r=0 when accepted null hypothesis)
#' # second value of the vector is the p-value and third value is the tests statistic value
#'
#'@seealso \code{\link[homnormal]{Brown_Forsythe}}, \code{\link[homnormal]{Cat_GG}}, \code{\link[homnormal]{bdai}}, \code{\link[homnormal]{genp}}, \code{\link[homnormal]{slrt}}, \code{\link[homnormal]{levene}}
#' @references Chang, C. H., Pal, N., & Lin, J. J. (2017). A revisit to test the equality of variances of several populations. Communications in Statistics-Simulation and Computation, 46(8), 6360-6384.

Cat_LR=function(x1,x2,alfa=0.05,m=2000,table=TRUE,graph="none"){
  mutt=c()
  s2=c()
  n=c()
  for (i in factor(x2)) {
    n[as.numeric(i)]=(sum(x2==i))
  }
  k=length(n)
  for (i in 1:k) {
    mutt[i]=mean(x1[x2==i])
    s2[i]=var(x1[x2==i])
  }
s2pooled=sum(s2*(n-1))/sum(n-1)
T=-sum(n*(log(s2)-log(s2pooled)))
av2=datagen(n,rep(0,k),rep(s2pooled,k),m)
a=av2[1];b=av2[2]
mut=a[[1]];s2t=b[[1]]
nn=kronecker(matrix(data=1,nrow=m,ncol=1),matrix(data=n,ncol=k,nrow=1))
s2p=rowSums((nn-1)*s2t)/sum(n-1)
s2pp=kronecker(matrix(data=s2p,nrow=m,ncol=1),matrix(data=1,ncol=k,nrow=1))
Tc=-rowSums(nn*(log(s2t)-log(s2pp)))
p=mean(Tc>T)
p
r=(p<alfa)
i
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
  C3=c(mutt)
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