% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{proportionate}
\alias{proportionate}
\title{A Number Multiplies Proportions Under Limitations}
\usage{
proportionate(
  number,
  result_sum_strict = NULL,
  prop,
  prop_sum_strict = NULL,
  digits = 0
)
}
\arguments{
\item{number}{A number}

\item{result_sum_strict}{The limitation of the sum of results}

\item{prop}{The proportions to be multiplied by the number}

\item{prop_sum_strict}{The limitation of the sum of proportions}

\item{digits}{The number of decimal places}
}
\description{
A Number Multiplies Proportions Under Limitations
}
\examples{
## Save 1 decimal place
a <- proportionate(number = 355,
                   prop = c(0.2, 0.6, 0.15, 0.36),
                   digits = 1)
## The sum of the proportions is 1
b <- proportionate(number = 355,
                   prop = c(0.2, 0.6, 0.15, 0.05),
                   prop_sum_strict = 1,
                   digits = 1)
## Save 0 decimal place
c <- proportionate(number = 355,
                   prop = c(0.2, 0.6, 0.15, 0.05),
                   prop_sum_strict = 1,
                   digits = 0)
## The sum of the results is 355
d <- proportionate(number = 355,
                   result_sum_strict = 355,
                   prop = c(0.2, 0.6, 0.15, 0.05),
                   prop_sum_strict = 1,
                   digits = 0)
}
