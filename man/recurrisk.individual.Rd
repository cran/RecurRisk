\name{recurrisk.individual}
\alias{recurrisk.individual}
\title{
Estimate recurrence risk using individual survival data
}
\description{
A function to estimate the risk of recurrence using cancer registry disease-specific individual survival data.
}
\usage{
recurrisk.individual(data, stratum, covar, timevar, eventvar,

stagevar, stage.dist.value, link = "Log-logistic", adj.r = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
   The individual survival data set created using the case-listing unformatted CSV data exported from SEER*Stat software 
   or user-generated individual data by other statistical software.
}
  \item{stratum}{
   The user-defined strata variables which should be categorical coded as integers. 
   One or more strata can be defined. If there is no stratum, user should define it as NULL.
}
  \item{covar}{
   The covariates defined in the mixture cure survival model which should be categorical coded as integers. 
   One or more covariates can be defined. If there is no covariate, user should define it as NULL.
}
  \item{timevar}{
   The variable specified as follow-up time in the survival model.
}
  \item{eventvar}{
   The status indicator, 0=alive, 1=dead (due to cancer). Note that, dead events due to other causes 
   should be defined as censoring events.
}
  \item{stagevar}{
   The stage variable defined in SEER*Stat data. If there are more than 1 stage variable, 
   the user will need to define the one which contains the distant stage.
}
  \item{stage.dist.value}{
   The numeric value of distant stage from the values of stage variable.
}
  \item{link}{
   The latency distribution for the cure model (non-cured survival). The current version can handle Weibull and log-logistic distributions. 
   The default link is log-logistic.
}
  \item{adj.r}{
   The adjustment factor used to adjust the registry-based survival curves for sensitivity analysis. The default value is 1.
}
}
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

\note{
Note that, the cure fraction and parametric survival distribution are estimated using \code{flexsurvcure} function in R package flexsurvcure for individual data.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{data.individual}}
}
\examples{
#data("data.individual")
#timevar<-"time"
#eventvar<-"status"
#stagevar<-"stage"
#stage.dist.value<-3
#stratum<-c("stage","agegroup")
#covar<-"yeargroup"
#link<-"Log-logistic"
#adj.r<-1.2
#out<-recurrisk.individual(data.individual, stratum, covar, timevar, eventvar,
#                          stagevar, stage.dist.value, link, adj.r)

}

