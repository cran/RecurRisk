\name{data.cansurv}
\alias{data.cansurv}
\docType{data}
\title{
CanSurv CSV format output 
}
\description{
An example CSV format output from CanSurv software including parameters for the mixture cure survival model.
}
\usage{data("data.cansurv")}
\format{
  A data frame with 99 observations on the following 15 variables. Two strata "SEER_historic_stage_LRD"
  and "Age_breast", and one covariate "Year_of_diagnosis_BC" are defined in the model.
  \describe{
    \item{\code{Link}}{The parametric survival distribution among those not cured specified in CanSurv.}
    \item{\code{Stratum_1__SEER_historic_stage_LRD}}{The stratum SEER_historic_stage_LRD value.}
    \item{\code{Stratum_2__Age_breast}}{The stratum Age_breast value.}
    \item{\code{estimate}}{The estimated parameters for the mixture cure survival model.}
    \item{\code{c_int}}{A numeric vector in covariance matrix wrt cure intercept.}
    \item{\code{c_1__Year_of_diagnosis_BC__0}}{A numeric vector in covariance matrix wrt cure parameter for Year_of_diagnosis_BC baseline.}
    \item{\code{c_1__Year_of_diagnosis_BC__1}}{A numeric vector in covariance matrix wrt cure parameter for Year_of_diagnosis_BC value 1.}
    \item{\code{c_1__Year_of_diagnosis_BC__2}}{A numeric vector in covariance matrix wrt cure parameter for Year_of_diagnosis_BC value 2.}
    \item{\code{c_1__Year_of_diagnosis_BC__3}}{A numeric vector in covariance matrix wrt cure parameter for Year_of_diagnosis_BC value 3.}
    \item{\code{u_int}}{A numeric vector in covariance matrix wrt mu intercept.}
    \item{\code{u_1__Year_of_diagnosis_BC__0}}{A numeric vector in covariance matrix wrt mu parameter for Year_of_diagnosis_BC baseline.}
    \item{\code{u_1__Year_of_diagnosis_BC__1}}{A numeric vector in covariance matrix wrt mu parameter for Year_of_diagnosis_BC value 1.}
    \item{\code{u_1__Year_of_diagnosis_BC__2}}{A numeric vector in covariance matrix wrt mu parameter for Year_of_diagnosis_BC value 2.}
    \item{\code{u_1__Year_of_diagnosis_BC__3}}{A numeric vector in covariance matrix wrt mu parameter for Year_of_diagnosis_BC value 3.}
    \item{\code{s_int}}{A numeric vector in covariance matrix wrt sigma intercept.}
  }
}
\details{
For each stratum, the variable order in the estimate column should be cure intercept, cure coefficients for each level of the covariate, 
mu intercept, mu coefficients for each level of the covariate, sigma intercept.
}

\examples{
data("data.cansurv")
}
\keyword{datasets}