\name{ovaHRDscar}
\alias{ovaHRDscar}
\title{get.ovaHRDscars}
\usage{
get.ovaHRDscars(a, chrominfo ="grch38")
}
\description{
Get the number of allelic imbalances associated with HRD in high grade serus ovarian cancer with the function get.ovaHRDscars
}
\examples{
a <-read.table("/examples/segments.txt", header=T)
get.ovaHRDscars(a, chrominfo ="grch38")
}
