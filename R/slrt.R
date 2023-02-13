#' Standartized Likelihood Ratio Test for Homogeniety
#'
#' Tests the homogeniety of variances for more than two normal groups
#' using standartized likelihood ratio test.
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
#'@return if table is TRUE, then it gives a detailed table,
#'else it gives a vector of r value(r=1 when null hypothesis was rejected and r=0 when null hypothesis was accepted)
#' p-value and test statistic value.
#'
#' @examples
#'    data(FH_data)
#'    x1=FH_data$SurvivalTime
#'    x2=FH_data$HospitalNo
#'    slrt(x1,x2)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'    slrt(x1,x2,alfa=0.10)
#'    readline(prompt = "Pause. Press <Enter> to continue...")
#'     slrt(x1,x2,alfa=0.10,table=FALSE)
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     slrt(x1,x2,alfa=0.10,table=FALSE,graph=FALSE)
#'     readline(prompt = "Pause. Press <Enter> to continue...")
#'     slrt(x1,x2,alfa=0.10,table=FALSE,graph="none")
#' # ---THIS VERSION IS ESPECIALLY USEFUL FOR COMPARISON STUDIES BY SIMULATION---
#' #    #first value of the vector is r value(r=1 when rejected and r=0 when accepted null hypothesis)
#' # second value of the vector is the p-value and third value is the tests statistic value
#' @seealso \code{\link[homnormal]{Brown_Forsythe}}, \code{\link[homnormal]{datagen}},\code{\link[homnormal]{levene}},
#' \code{\link[homnormal]{Cat_LR}}, \code{\link[homnormal]{genp}}
#' @references Gökpınar, E. (2020). Standardized likelihood ratio test for homogeneity of variance of several normal populations. Communications in Statistics-Simulation and Computation, 1-11.
slrt<-function(x1,x2,alfa=0.05,table=TRUE,graph="none"){
  mut=c()
  s2tt=c()
  n=c()
  for (i in factor(x2)) {
    n[as.numeric(i)]=(sum(x2==i))
  }
  k=length(n)
  for (i in 1:k) {
    mut[i]=mean(x1[x2==i])
    s2tt[i]=var(x1[x2==i])
  }
  s2pooled=sum(s2tt*(n-1))/sum(n-1)
  k=length(n)
  E=-sum(n*(digamma((n-1)/2)-log(n)))+sum(n)*(digamma(sum(n-1)/2)-log(sum(n)))
  V=sum(n^2*trigamma((n-1)/2))-sum(n)^2*trigamma(sum(n-1)/2)
  Tl=-sum(n*(log(s2tt)-log(s2pooled)))
  T=(Tl-E)/sqrt(V)*sqrt(2*k-2)+k-1;
  p=1-pchisq(T,k-1)
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
    return(c(r,p,T))
  }
  else if (table==TRUE)
  {
    
    C1=1:length(n)
    C2=c(n)
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