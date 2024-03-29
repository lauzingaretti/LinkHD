% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LinkData.R
\name{LinkData}
\alias{LinkData}
\title{LinkData: multiple heterogeneous dataset integration}
\usage{
LinkData(
  Data,
  Distance = c(),
  Center = FALSE,
  Scale = FALSE,
  CorrelVector = TRUE,
  nCluster = 0,
  cl_method = "pam"
)
}
\arguments{
\item{Data}{should be a list of dataframes or ExpressionSet data
with the same length of the number of tables to be integrate.
In each dataframe, the Observations (common elements on Statis)
should be in rows and the variables should be in columns.
Data also might be a MultiAssayExperiment object
from MultiAssayExperiment package, a software for
multi-omics experiments integration in Bioconductor.}

\item{Distance}{Vector indicating which distance (including scalar product)
should be applied to each study. If is missing,
the scalar product is used. The vector lenght must be equal
to the length of Data. Distance options: ScalarProduct, euclidean,
manhattan, canberra, pearson, pearsonabs, spearman,
spearmanabs, mahalanobis, BrayCurtis distance (please, use option Bray).
For binary data, the distance can be jaccard,
simple_matching, sokal_Sneath, Roger_Tanimoto, Dice,
Hamman, Ochiai, Phi_Pearson, 'Gower&Legendre.
Note that, use pre-processing option as compositional and
Euclidean is the same than use Aitchison distance for compositional data.}

\item{Center}{Logical. If TRUE, the data frame
is centered  by the mean. By default is FALSE.
If you have tables with different characteristics (continous phenotypes, frecuencies,
compositional data), we strongly recomendate normalize
datasets as a previous step through DataProcessing option.}

\item{Scale}{A logical value indicating whether the column vectors should be
standardized by the rows weight, by default is FALSE.
Note that all data into the list will be scaled.
If you don't need normalizing all data, you
could set this parameter as False and perform the normalization step
externally by using DataProcessing function.
If you have tables with different characteristics (continous phenotypes, frecuencies,
compositional data), we strongly recomendate normalize datasets
as a previous step through DataProcessing option.}

\item{CorrelVector}{Logical. If TRUE (default), the RV matrix is
computed using vectorial correlation, else
the Hilbert-Smith distance is used.}

\item{nCluster}{this variable indicates if common
elements on the dataset should be grouped (by default is zero, i.e. no-cluster).}

\item{cl_method}{categorical (pam or kmeans). pam is a robust
version of classical kmeans algorithm.}
}
\value{
\item{LinkData}{DistStatis class object with the
corresponding completed slots according to the given model}
}
\description{
Integrating multiple Heterogeneous Datasets
stored into a list. This function makes Statis using Distances options.
Statis is part of the PCA family and is based on singular value decomposition
(SVD) and the generalized singular value decomposition (GSVD) of a matrix.
This methodology aims to analyze several data sets of
variables that were collected on the same set of observations.
Originally, the comparisons were drawn from the compute
of the scalar product between the different tables.
In our approach, the condition is relaxing
allowing  the incorporation of different distances.
}
\examples{
{
data(Taraoceans)
pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
,as.data.frame(Taraoceans$pro.NOGs))
TaraOc_1<-scale(TaraOc[[1]])
Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
function(x){DataProcessing(x,Method='Compositional')})
colnames(Normalization[[1]])=pro.phylo
colnames(Normalization[[2]])=Taraoceans$GO
TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
TaraOc<-lapply(TaraOc,as.data.frame)
Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
}

}
\references{
\enumerate{
 \item Escoufier, Y. (1976). Operateur associe a un tableau de donnees.
  Annales de laInsee, 22-23, 165-178.
 \item Escoufier, Y. (1987). The duality diagram: a means
  for better practical applications. En P. Legendre & L. Legendre (Eds.),
  Developments in Numerical Ecology, pp. 139-156,
  NATO Advanced Institute, Serie G. Berlin: Springer.
 \item L'Hermier des Plantes, H. (1976). Structuration des
 Tableaux a Trois Indices de la Statistique. [These de Troisieme Cycle].
 University of Montpellier, France.
}
}
\author{
Laura M Zingatetti
}
