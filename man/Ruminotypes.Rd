\docType{data}
\name{Ruminotypes}
\alias{Ruminotypes}
\title{Whole Ruminotypes dataset communities}
\description{Ruminotypes dataset contains communities (16_S, 18_S and Archaea) measuares from 65 loose-housed lactating Holstein cows. The study aims evaluating the relationships between communities and methane emmision yield }
\format{A list with seven components:
\describe{
\item{16_S}{matrix  with  61   rows  and  1198  columns. Each  row  represents  a  sample  and  each column represent one normalized OTU.}
\item{Archea}{a  matrix with 61 rows (samples) and 453 normalized columns (Archea IDs).}
\item{18_S}{a matrix with 61 rows and 107 normalized columns (protozoa level).}
\item{phenotype}{ data frame with 61 rows and 5 columns representing methane emission levels with a set of corrections.}
\item{Taxa_16S}{ a matrix with 1198 rows and 9 columns indicating the Taxa information for 16_S OTU.}
\item{Taxa_18S}{a matrix with 112 rows and 18 columns indicating the Taxa information for 18_S samples.}
\item{Taxa_Archea}{a matrix with 453 rows and 7 columns indicating the Taxa information for archea samples.}
}
}
\keyword{datasets}
\references{
\enumerate{
Ramayo-Caldas Y, Zingaretti LM, Bernard A, Estellé J,
Popova M, Pons N, Bellot P, Mach N, Rau A, Roume H,
Pérez-Enciso M, Faverdin N, Edouard N, Dusko S,
Morgavi DP, Renand G. Identification of rumen microbial
biomarkers linked to methane emission in Holstein dairy cows
In press.}}

\examples{
data(Ruminotypes)
}
\keyword{datasets}




