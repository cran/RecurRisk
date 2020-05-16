\name{recurrisk.group}
\alias{recurrisk.group}
\title{
Estimate recurrence risk using grouped survival data
}
\description{
A function to estimate the risk of recurrence using cancer registry disease-specific grouped survival data.
}
\usage{
recurrisk.group(data, data.cansurv, stagevar, stage.dist.value, adj.r = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
   The SEER*Stat cause-specific or relative group survival data frame returned by function \code{read.seerstat}.
}
  \item{data.cansurv}{
   The data frame including parameter estimation and covariance matrix for the mixture cure survival model from Cansurv CSV format output. 
   It also contains the information on strata and covariates. 
}
  \item{stagevar}{
   The stage variable containing the distant stage from SEER*Stat data. 
}
  \item{stage.dist.value}{
   The numeric value of distant stage.
}
  \item{adj.r}{
   The adjustment factor used to adjust the registry-based survival curves for sensitivity analysis. The default value is 1.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  A data frame containing the following items.
  \item{link}{The parametric survival distribution among those not cured specified in CanSurv.}
  \item{cure}{The cure fraction estimated from the mixture cure survival model.}
  \item{lambda}{The estimated scale parameter of the survival distribution for those not cured.}
  \item{k}{The estimated shape parameter of the survival distribution for those not cured.}
  \item{theta}{The exponential hazard of the time from recurrence to cancer death.}
  \item{surv_curemodel}{The survival estimated from the mixture cure survival model.}  
  \item{surv_notcure}{The estimated survival for the non-cured fraction.}
  \item{median_surv_notcured}{The median survival time for the non-cured fraction.}
  \item{s1_numerical}{The numerical estimated survival to recurrence (recurrence-free survival) for the non-cured fraction.}
  \item{G_numerical}{The numerical estimated survival to recurrence.}
  \item{CI_numerical}{1-G_numerical, the numerical estimated cumulative incidence of recurrence which is the probability of progressing to cancer recurrence.}
  \item{s1_analytical}{The analytical estimated survival to recurrence (recurrence-free survival) for the non-cured fraction.}
  \item{G_analytical}{The analytical estimated survival to recurrence.}
  \item{CI_analytical}{1-G_analytical, the analytical estimated cumulative incidence of recurrence.}
  \item{se_CI_analytical}{The standard error of CI_analytical.}
  \item{obs_surv}{The observed survival from SEER*Stat.}
  \item{obs_dist_surv}{The observed survival for distant stage from SEER*Stat.}
}
%%\note{
%%  ~~further notes~~
%%}
\seealso{
\code{\link{data.group}}, \code{\link{read.seerstat}}.
}
\examples{
data("data.group")
data("data.cansurv")
stagevar<-"SEER_historic_stage_LRD"
stage.dist.value<-2
adj.r<-1.5
out<-recurrisk.group(data.group, data.cansurv, stagevar, stage.dist.value, adj.r)
}
