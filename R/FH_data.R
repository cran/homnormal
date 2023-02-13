#' Fleming and Harrington Data
#' 
#' The data related to survival times of patients was collected from 4 hospitals, 
#' which was a part of the data by given Fleming and Harrington(1991). The data contain failure time of the
#' patients.
#' 

#' @docType data 
#'
#' @usage data(FH_data)
#' 
#' @format A dataframe with 21 rows 2 variables
#' \describe{
#' \item{HospitalNo}{Hospital No}
#' \item{SurvivalTime}{ Survival Time of Patients}
#' }
#' @source {T.R. Fleming and D.P. Harrington,  Counting processes and survival analysis.  Wiley Online Library, Vol. 8., 1991.}
#' @examples 
#'  data("FH_data")
#'    x1=FH_data$SurvivalTime
#'    x2=FH_data$HospitalNo
"FH_data"