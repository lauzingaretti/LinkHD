#' Class \code{DistStatis} DistStatis S4 class (linkHD:Multiple Heterogeneous Dataset Integration)
#' Statis with Distance options implementation.

#' @section Features:
#' \enumerate{
#'   \item DistStatis (implements Statis method incorporating Distance options to integrate multiple heterogeneous datasets)
#'   \item Implement a LM (Linear Model) to variable selection
#'   \item Incorporate a method to variable clustering
#'   \item Incorporate some visualization tools: Compromise visualization, Relationship-visualization
#' }
#'
#' @section Fields:
#' \itemize{
#'  \item RV: Vectorial Correlation Matrix between studies.
#'  \item Inertia.RV: Inertia (\%) explained for all tables.
#'  \item Euclid.Im: Euclidean Image of all studies.
#'  \item Inertia.comp: Inertia (\%) explained for all dimensions of compromise matrix.
#'  \item Compromise.Coords: Projection of all observations in compromise (Coords).
#'  \item Compromise.Matrix: Compromise Matrix from statis methodology.
#'  \item RQO: Representation Quality of observations in compromise matrix.
#'  \item TableProjections: Projection of each table on Compromise configuration
#' }
#'
#'  @slot RV: Vectorial Correlation Matrix between studies.
#'  @slot  Inertia.RV: Inertia (\%) explained for all tables.
#'  @slot  Euclid.Im: Euclidean Image of all studies.
#'  @slot  Inertia.comp: Inertia (\%) explained for all dimensions of compromise matrix.
#'  @slot  Compromise.Coords: Projection of all observations in compromise (Coords).
#'  @slot  Compromise.Matrix: Compromise Matrix from statis methodology.
#'  @slot  RQO: Representation Quality of observations in compromise matrix.
#'  @slot  TableProjections: Projection of each table on Compromise configuration
#'
#' @section DistStatis-general-functions:
#' \describe{
#'  \item{DistStatis}{Getters for their respective slots.}
#' }
#'
#'  @author Laura M Zingaretti
#'
#' @examples
#' {
#' showClass('DistStatis')
#' }
#' @docType methods
#' @name DistStatis-class
#' @rdname DistStatis-Class
#' @export
#'
#'
#'
setClass(Class = "DistStatis", slots = c(Inertia.RV = "data.frame", RV = "matrix", Euclid.Im = "data.frame",
    Inertia.comp = "matrix", Compromise.Coords = "data.frame", Compromise.Matrix = "data.frame", RQO = "data.frame",
    TableProjections = "list"), prototype = prototype(Inertia.RV = data.frame(), RV = matrix(), Euclid.Im = data.frame(),
    Inertia.comp = matrix(), Compromise.Coords = data.frame(), Compromise.Matrix = data.frame(), RQO = data.frame(),
    TableProjections = list()))


setGeneric(name = "Inertia_RV", def = function(x) {
    standardGeneric("Inertia_RV")
})
setGeneric(name = "correl", def = function(x) {
    standardGeneric("correl")
})
setGeneric(name = "Euclid_Im", def = function(x) {
    standardGeneric("Euclid_Im")
})
setGeneric(name = "Inertia_comp", def = function(x) {
    standardGeneric("Inertia_comp")
})
setGeneric(name = "compromise_coords", def = function(x) {
    standardGeneric("compromise_coords")
})
setGeneric(name = "Compromise_matrix", def = function(x) {
    standardGeneric("Compromise_matrix")
})
setGeneric(name = "RQO", def = function(x) {
    standardGeneric("RQO")
})
setGeneric(name = "Trajectories", def = function(x) {
    standardGeneric("Trajectories")
})

#' @title Inertia_RV
#' @description Accessor to Inertia_RV from LinkData output.
#' @param x an object from DistSatis class.
#' @name Inertia_RV
#' @aliases  Inertia_RV,DistStatis-method
#' @docType methods
#' @rdname Inertia_RV-methods
#' @export
#' @return  Inertia_RV explained inertia for RV matrix from LinkData object
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Inertia_RV(Output)
#' }
setMethod("Inertia_RV", "DistStatis", function(x) (x@Inertia.RV))

#' @title correl
#' @description Accessor to RV (Vectorial correlation coefficient) from LinkData output.
#' @param x an object from DistSatis class.
#' @name correl
#' @aliases  correl,DistStatis-method
#' @docType methods
#' @rdname correl-methods
#' @export
#' @return  RV correlation coefficient for each input table to  LinkData function
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'correl(Output)
#' }
setMethod("correl", "DistStatis", function(x) (x@RV))
#' @title Euclid_Im
#' @description Accessor to the Observations Image Euclidean, i.e. the projections from LinkData output.
#' @param x an object from DistSatis class.
#' @name Euclid_Im
#' @aliases  Euclid_Im,DistStatis-method
#' @docType methods
#' @rdname Euclid_Im-methods
#' @export
#' @return  Euclid_Im Euclidean image of the input tables in LinData function.
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Euclid_Im(Output)
#' }
setMethod("Euclid_Im", "DistStatis", function(x) (x@Euclid.Im))

#' @title Inertia_comp
#' @description Accessor to explained inertia of compromise axis from LinkData output.
#' @param x an object from DistSatis class.
#' @name Inertia_comp
#' @aliases  Inertia_comp,DistStatis-method
#' @docType methods
#' @rdname Inertia_comp-methods
#' @export
#' @return  Inertia_comp explained inertia for Compromise matrix from LinkData object
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Inertia_comp(Output)
#' }
setMethod("Inertia_comp", "DistStatis", function(x) (x@Inertia.comp))

#' @title compromise_coords
#' @description Accessor to compromise coordinates from LinkData output.
#' @param x an object from DistSatis class.
#' @name compromise_coords
#' @aliases  compromise_coords,DistStatis-method
#' @docType methods
#' @rdname compromise_coords-methods
#' @export
#' @return  compromise_coords coordinates of observations in the compromise configuration from LinkData function
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'compromise_coords(Output)
#' }
setMethod("compromise_coords", "DistStatis", function(x) (x@Compromise.Coords))

#' @title Compromise_matrix
#' @description Accessor to Compromise Matrix from LinkData output.
#' @param x an object from DistSatis class.
#' @name Compromise_matrix
#' @aliases  Compromise_matrix,DistStatis-method
#' @docType methods
#' @rdname Compromise_matrix-methods
#' @export
#' @return  Compromise_matrix: Compromise matrix from LinkData object
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Compromise_matrix(Output)
#' }
setMethod("Compromise_matrix", "DistStatis", function(x) (x@Compromise.Matrix))
#' @title RQO
#' @description Accessor to RQO (% of individuals representation) from LinkData output.
#' @param x an object from DistSatis class.
#' @name RQO
#' @aliases  RQO,DistStatis-method
#' @docType methods
#' @rdname RQO-methods
#' @export
#' @return  RQO Representation Quality of the observations in the compromise configuration from LinkData object
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'RQO(Output)
#' }
setMethod("RQO", "DistStatis", function(x) (x@RQO))

#' @title Trajectories
#' @description Accessor to projections into the common configuration, i.e. compromise of each input table from LinkData output.
#' @param x an object from DistSatis class.
#' @name Trajectories
#' @aliases  Trajectories,DistStatis-method
#' @docType methods
#' @rdname Trajectories-methods
#' @export
#' @return  Trajectories contains a list of the projections
#' of each input table into the common configuration, i.e. the compromise
#' from LinkData object
#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo)
#',as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Trajectories(Output)
#' }
setMethod("Trajectories", "DistStatis", function(x) (x@TableProjections))

