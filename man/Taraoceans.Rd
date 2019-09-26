\name{Taraoceans}
\alias{Taraoceans}
\docType{data}
\title{TARA Oceans datasets}
\description{
 TARA Oceans was an expedition allowing to the
 study of plankton communities and their interactions with enviromental variables.
 This dataset was taken from mixkernel package (https://cran.r-project.org/web/packages/mixKernel/index.html).
 Data constists on  139 prokaryotic-enriched samples collected from 68 stations and spread across threedepth layers:
 the surface (SRF), the deep chlorophyll maximum (DCM) layer and the mesopelagic(MES) zones.
 Samples were located in height different oceans or seas: Indian Ocean (IO),
 Mediter-ranean Sea (MS), North Atlantic Ocean (NAO), North Pacific Ocean (NPO),
 Red Sea (RS), SouthAtlantic Ocean (SAO), South Pacific Ocean (SPO) and South Ocean (SO).
 }
\usage{data("Taraoceans")}
\format{A list with seven components:
\describe{
\item{phychemdata}{matrix  with  139  rows  and  22  columns.   Each  row  represents  a  sample  and  each column
   an environmental variable.}
\item{pro.phylo}{a  matrix with 139 rows (samples) and 356 columns (prokaryotic OTUs).}
\item{taxonomy}{a matrix with 356 rows (prokaryotic OTUs) and 6 columns indicating the taxonomy of each OTU.}
\item{phylogenetic.tree}{a phylo object (see package ’ape’) representing the prokaryotic OTUs}
\item{pro.NOGs}{a matrix with 139 rows (samples) and 638 columns (NOGs).}
\item{GO}{a list with the names of Gene Ontologies.}
\item{sample}{a  list  containing  three  following  entries  (all  three  are  character  vectors):name(samplename),ocean(oceanic region of the sample) and depth(sample depth)}
}
}
\references{
 Sunagawa S., Coelho L.P., Chaffron S., Kultima J.R.,
 Labadie K., Salazar F., Djahanschiri B., ZellerG.,
 Mende D.R., Alberti A., Cornejo-Castillo F.,
 Costea P.I., Cruaud C., d’Oviedo F., Engelen S.,Ferrera I.,
 Gasol J., Guidi L., Hildebrand F., Kokoszka F.,
 Lepoivre C., Lima-Mendez G., PoulainJ., Poulos B.,
 Royo-Llonch M., Sarmento H., Vieira-Silva S., Dimier C.,
 Picheral M., Searson S.,Kandels-Lewis S.,TaraOceans coordinators,
 Bowler C., de Vargas C., Gorsky G., Grimsley N.,Hingamp P., Iudicone D.,
 Jaillon O., Not F., Ogata H., Pesant S., Speich S., Stemmann L., SullivanM.,
 Weissenbach J., Wincker P., Karsenti E., Raes J., Acinas S. and Bork P. (2015).
 Structure andfunction of the global ocean microbiome.Science,348, 6237 }


\examples{
data(Taraoceans)
}
\keyword{datasets}




