\name{data.individual}
\alias{data.individual}
\docType{data}
\title{
Individual survival data
}
\description{
An example individual survival data including the least required variables for \code{recurrisk.individual} function, which is created using the case-listing unformatted SEER*Stat CSV data file.
}
\usage{data("data.individual")}
\format{
  A data frame with 421198 observations on the following 5 variables. 
  \describe{
    \item{\code{stage}}{Stage varible with format 0=Localized, 1=Regional, 2=Distant.}
    \item{\code{yeargroup}}{Year variable with format 0=1985-89, 1=1990-94, 2=1995-99, 3=2000-13.}
    \item{\code{agegroup}}{Age variable with format 0=15-59, 1=60-69, 2=70-84.}
    \item{\code{time}}{Follow-up years.}
    \item{\code{status}}{The status indicator, 0=alive, 1=dead due to cancer. }
  }
}
\details{
Note that, dead events due to other causes should be defined as censoring events.
}

\examples{
data("data.individual")
}
\keyword{datasets}
