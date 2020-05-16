\name{data.group}
\alias{data.group}
\docType{data}
\title{
SEER*Stat group survival data
}
\description{
An example SEER*Stat data set including the cause-specific group survival. 
}
\usage{data("data.group")}
\format{
  A data frame with 803 observations on the following 12 variables.
  \describe{
    \item{\code{Page_type}}{Page type with format 0=Life Page, 1=Summary Page, etc..}
    \item{\code{SEER_historic_stage_LRD}}{Stage varible with format 0=Localized, 1=Regional, 2=Distant.}
    \item{\code{Year_of_diagnosis_BC}}{Year variable with format 0=1985-89, 1=1990-94, 2=1995-99, 3=2000-13.}
    \item{\code{Age_breast}}{Age variable with format 0=15-59, 1=60-69, 2=70-84.}
    \item{\code{Interval}}{Follow-up years.}
    \item{\code{Alive_at_Start}}{The number of people who are alive at the beginning of interval.}
    \item{\code{Died}}{The number of cancer deaths in the interval.}
    \item{\code{Lost_to_Followup}}{The number of people lost to follow-up or dying from other causes in the interval.}
    \item{\code{CauseSpecific_Survival_Interval}}{The cause-specific survival for the interval.}
    \item{\code{CauseSpecific_Survival_Cum}}{The cause-specific cumulative survival.}
    \item{\code{CauseSpecific_SE_Interval}}{The standard error for cause-specific interval survival.}
    \item{\code{CauseSpecific_SE_Cum}}{The standard error for cause-specific cumulative survival.}
  }
}
\examples{
data("data.group")
}
\keyword{datasets}