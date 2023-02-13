#' Bahandary-Dai Test for Homogeniety
#'
#' Tests the homogeniety of variances for more than two normal groups
#' using Bahandary-Dai test.
#' @export
#' @param x1 a numeric matrix containing the values of groups.
#'
#' @param x2  numeric matrix containing the values of group numbers.
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
#'
#'    data(FH_data)
#'    x1=FH_data$SurvivalTime
#'    x2=FH_data$HospitalNo
#'    bdai(x1,x2)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'    bdai(x1,x2,alfa=0.10)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'     bdai(x1,x2,alfa=0.10,table=FALSE)
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     bdai(x1,x2,alfa=0.10,table=FALSE,graph="raw")
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     bdai(x1,x2,alfa=0.10,table=FALSE,graph="none")
#' # ---THIS VERSION IS ESPECIALLY USEFUL FOR COMPARISON STUDIES BY SIMULATION---
#' #    #first value of the vector is r value(r=1 when rejected and r=0 when accepted null hypothesis)
#' # second value of the vector is the p-value and third value is the tests statistic value
#'
#' @seealso \code{\link[homnormal]{Brown_Forsythe}}, \code{\link[homnormal]{Cat_GG}}, \code{\link[homnormal]{Cat_LR}}, \code{\link[homnormal]{genp}}, \code{\link[homnormal]{slrt}}, \code{\link[homnormal]{levene}}
#' @references Bhandary, M., & Dai, H. (2008). An alternative test for the equality of variances for several populations when the underlying distributions are normal. Communications in Statistics-Simulation and Computation, 38(1), 109-117.
bdai=function(x1,x2,alfa=0.05,table=TRUE,graph="none"){
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
sp=(sum(s2*(n-1))-s2*(n-1))/(sum(n-1)-(n-1))
F=s2/sp;
Fh=c(F,1/F)
ph=c(pf(F,n-1,(sum(n-1)-(n-1)),lower.tail = FALSE), pf(1/F,(sum(n-1)-(n-1)),n-1,lower.tail = FALSE))
sp=sort(ph)
spy=sort(ph)*(2*k)/(1:(2*k))
p=min(spy)
Th=Fh[ph==sp[spy==p]]
r=(p<alfa)

if(graph=="raw") {
  boxplot(x1~x2,xlab = "", ylab = "")
}
else if(graph=="centered") {
  boxplot(x1-rep(mut,ne)~x2,xlab = "", ylab = "")
}
else if (graph=="none")
{}


if (table==FALSE) {
  return(c(r,p,Th))
}
else if (table==TRUE)
{
  C1=1:length(n)
  C2=c(n)
  C3=c(mutt)
  C4=c(s2)
  C5=c(rep(NA,ifelse(k%%2==0,ceiling(k/2),floor(k/2))),Th,rep(NA,floor(k/2-1)))
  C6=c(rep(NA,ifelse(k%%2==0,ceiling(k/2),floor(k/2))),p,rep(NA,floor(k/2-1)))
  HT<-hux(Group_No=C1,Sample_Size=C2,Sample_Mean=C3,Sample_Var=C4,Test_Stat=C5,p_value=C6)
  align(HT)[,1:5]          <- 'center'
  number_format(HT)[,3:6]      <- 2
  number_format(HT)[,2]      <- 0
 
  
  print_md(HT, header = TRUE)
  return(HT)
}
}