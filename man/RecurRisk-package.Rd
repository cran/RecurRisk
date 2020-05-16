\name{RecurRisk-package}
\alias{RecurRisk-package}
\alias{RecurRisk}
\docType{package}
\title{
\packageTitle{RecurRisk}
}
\description{
This package is to project the risk of recurrence using disease-specific survival data.
}
\details{
This package has been developed to estimate the risk of progressing to distant recurrence using disease-specific survival from cancer registry data. 
The disease-specific survival is assessed via cause-specific survival or relative survival using SEER*Stat software. 
The cause-specific survival or relative survival is assumed to follow a mixture-cure model and the risk of recurrence is inferred from the survival among the non-cured fraction. 
The cure fraction and parametric survival distribution among those not cured are estimated using CanSurv software for group data and R flexsurvcure package for individual data. 
The current version can handle Weibull and log-logistic distributions for the non-cured survival.
}
\author{
Fanni Zhang <zhangf@imsweb.com>
}
\references{
Mariotto AB, Zou Z, Zhang F, Howlader N, Kurian AW, Etzioni R. 
Can We Use Survival Data from Cancer Registries to Learn about Disease Recurrence? The Case of Breast Cancer. 
Cancer Epidemiol Biomarkers Prev. 2018 Nov; 27(11):1332-1341.

Yu B, Tiwari RC, Cronin KA, McDonald C, Feuer EJ. 
CANSURV: A Windows program for population-based cancer survival analysis. 
Comput.Methods Programs Biomed. 2005;80(3):195-203.

Yu B, Tiwari RC, Cronin KA, Feuer EJ. 
Cure fraction estimation from the mixture cure models for grouped survival data. 
Stat.Med. 2004;23(11):1733-47.

De Angelis R, Capocaccia R, Hakulinen T, Soderman B, Verdecchia A. 
Mixture models for cancer survival analysis: application to population-based data with covariates. 
Stat. Med. 1999;18(4):441-54.

Capocaccia R. Relationships between incidence and mortality in non-reversible diseases. 
Stat.Med. 1993;12(24):2395-415.

Li CS, Taylor JMG, Sy JP. Identifiability of cure models. 
Statistics & Probability Letters. 2001;54(4):389-95.

Howlader N, Ries LA, Mariotto AB, Reichman ME, Ruhl J, Cronin KA. 
Improved estimates of cancer-specific survival rates from population-based data. 
J Natl Cancer Inst. 2010;102(20):1584-98.
}
