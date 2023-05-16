\name{temp_graz}
\alias{temp_graz}
\docType{data}
\title{
Graz Temperature
}
\description{
Half-hourly temperature data collected in Graz (Austria) between July 1 and
August 31, 2022.
}
\usage{data("temp_graz")}
\format{
  A data frame with 2976 observations on the following 3 variables.
  \describe{
   \item{time}{Time}
   \item{temp_east}{Temperature (C°) measured in the east of Graz.}
   \item{temp_south}{Temperature (C°) measured in the south of Graz.}
  }
}

\source{
Land Steiermark (\href{https://data.steiermark.at/}{data.steiermark.gv.at}); licensed under \href{https://creativecommons.org/licenses/by/4.0/deed.en}{CC-BY-4.0}.
}

\examples{
data <- temp_graz
plot(as.POSIXct(data$time), data$temp_east, type = "l")
}
\keyword{datasets}
