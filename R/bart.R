#' Bartlett Test for Homogeniety
#'
#' Tests the homogeniety of variances for more than two normal groups.
#'
#' @export
#'
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
#' 
#'@return if table is TRUE, then it gives a detailed table,
#'else it gives a vector of r value(r=1 when null hypothesis was rejected and r=0 when null hypothesis was accepted)
#' p-value and test statistic value.
#'
#' @examples
#'     data(FH_data)
#'    x1=FH_data$SurvivalTime
#'    x2=FH_data$HospitalNo
#'    bart(x1,x2)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'    bart(x1,x2,alfa=0.10)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'     bart(x1,x2,alfa=0.10,table=FALSE)
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     bart(x1,x2,alfa=0.10,table=FALSE,graph="centerized")
#'      readline(prompt = "Pause. Press <Enter> to continue...")
#'     bart(x1,x2,alfa=0.10,table=FALSE,graph="none")
#' # ---THIS VERSION IS ESPECIALLY USEFUL FOR COMPARISON STUDIES BY SIMULATION---
#' #    #first value of the vector is r value(r=1 when rejected and r=0 when accepted null hypothesis)
#' # second value of the vector is the p-value and third value is the tests statistic value
#' 
#' @seealso \code{\link[homnormal]{levene}} \code{\link[homnormal]{Brown_Forsythe}}, \code{\link[homnormal]{Cat_GG}}, \code{\link[homnormal]{Cat_LR}}, \code{\link[homnormal]{genp}}, \code{\link[homnormal]{slrt}}, \code{\link[homnormal]{bdai}}
#' @references Bartlett, M. S. (1937). "Properties of sufficiency and statistical tests". Proceedings of the Royal Statistical Society, Series A 160, 268–282 JSTOR.

bart<-function(x1,x2,alfa=0.05,table=TRUE,graph="none"){
  ne=c()
  for (i in factor(x2)) {
    ne[as.numeric(i)]=(sum(x2==i))
  }
  xm=c()
  s2t=c()
  mut=c()
  s2tt=c()
  k=length(ne)
  for (i in 1:k) {
    mut[i]=mean(x1[x2==i])
    s2tt[i]=var(x1[x2==i])
  }

  
  d=ne-1
  xv=sum(d*s2tt)/(sum(ne)-k)
  t=sum(1/d)
  c=1+1/(3*(k-1))*(t-1/sum(d))
  B=log(xv)*sum(d)-sum(log(s2tt)*d)
  T=B/c
  p=1-pchisq(T,k-1)
  r=(p<alfa)
  c(r,p,T)
  
    
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
    C1=1:length(ne)
    C2=c(ne)
    C3=c(mut)
    C4=c(s2tt)
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