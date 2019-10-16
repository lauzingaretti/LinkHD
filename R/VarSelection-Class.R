#' @title class VarSelection
#' @description  Class \code{VarSelection} S4 class (linkHD: integrating multiple heterogeneous datasets)
#' VarSelection is a class to perform variable selection from a DistStatis object.

#' @section Features:
#' \enumerate{
#'   \item class to perform variable selection using Linear Regression Biplot onto the Compromise-Subespace
#'   \item This method allow variable selection and classification
#'   }
#'
#'
#' @section Fields:
#' \itemize{
#'\item  Variables return all the selected variables (and the frecuency of selection).
#'\item  Coordinates represent the coordenates (Betas coeeficients on LM) of the selected variables.
#'\item  VarTable data.frame indicating which table selected variables come from.
#'\item  values data.frame contains the R2 or pvalue (fdr) of selected variables (it depends of the Crit used).
#' }
#'
#' @section Accesors:
#' \itemize{
#'\item  Variables return all the selected variables (and the frecuency of selection).
#'\item  Coordinates represent the coordenates (Betas coeeficients on LM) of the selected variables.
#'\item  VarTable dataframe indicating the table that each selected variable comes from.
#'\item  values data.frame which contains the R2 or pvalue (fdr) of selected variables (it depends of the Crit used).
#' }
#' @section VarSelection-class-general-functions:
#'\describe{
#'\item{print}{Generated basic output for VarSelection class}
#' }
#'
#'  @author Laura M Zingaretti
#' @references
#' \enumerate{

#' \item Gabriel, K. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika 58(3), 453--467.

#' \item Gower, J. & Hand, D. (1996). Biplots, Monographs on statistics and applied probability. 54. London: Chapman and Hall., 277 pp.

#'
#' }
#'
#' @examples
#' {
#' showClass('VarSelection')
#' }
#'
#' @name VarSelection-class
#' @rdname VarSelection-Class
#' @export


setClass(Class = "VarSelection", slots = c(Variables = "character", Coordinates = "data.frame", VarTable = "data.frame", 
    sign_values = "data.frame"), prototype = prototype(Variables = c(), Coordinates = data.frame(), VarTable = data.frame(), 
    sign_values = data.frame()))




setGeneric(name = "Variables", def = function(x) {
    standardGeneric("Variables")
})
setGeneric(name = "Var_coordinates", def = function(x) {
    standardGeneric("Var_coordinates")
})
setGeneric(name = "VarTable", def = function(x) {
    standardGeneric("VarTable")
})
setGeneric(name = "sign_values", def = function(x) {
    standardGeneric("sign_values")
})

#' @title Variables
#' @description Accessor to selected Variables from VarSelection output.
#' @param x an object from VarSelection class.
#' @name Variables
#' @aliases  Variables,VarSelection-method
#' @docType methods
#' @rdname Variables-methods
#' @export
#' @return  Variables list of selected variables from VarSelection object
#' @examples
#'{
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo),
#'as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,
#'Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Selection<-VarSelection(Output,TaraOc,Crit='Rsquare',perc=0.95)
#'Variables(Selection)
#'}
setMethod("Variables", "VarSelection", function(x) (x@Variables))

#' @title Var_coordinates
#' @description Accessor to the coordinates projections into the compromise configuration of the selected variables from VarSelection output.
#' @param x an object from VarSelection class.
#' @name Var_coordinates
#' @aliases  Var_coordinates,VarSelection-method
#' @docType methods
#' @rdname Var_coordinates-methods
#' @export
#' @return  Var_Coordinates, Coordinates of variables into the common configuration,
#' i.e. the compromise from LinkData function
#' @examples
#'{
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo),
#'as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,
#'Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Selection<-VarSelection(Output,TaraOc,Crit='Rsquare',perc=0.95)
#'Var_coordinates(Selection)
#'}
setMethod("Var_coordinates", "VarSelection", function(x) (x@Coordinates))


#' @title VarTable
#' @description Accessor to Table with the selected variables from VarSelection output.
#' @param x an object from VarSelection class.
#' @name VarTable
#' @aliases  VarTable,VarSelection-method
#' @docType methods
#' @rdname VarTable-methods
#' @export
#' @return VarTable data.frame with the name of input tables in
#' the LinkData function
#' @examples
#'{
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo),
#'as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,
#'Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Selection<-VarSelection(Output,TaraOc,Crit='Rsquare',perc=0.95)
#'VarTable(Selection)
#'}
setMethod("VarTable", "VarSelection", function(x) (x@VarTable))


#' @title sign_values
#' @description Accessor to R2 or p values of the selected variables from VarSelection output.
#' @param x an object from VarSelection class.
#' @name sign_values
#' @aliases  sign_values,VarSelection-method
#' @docType methods
#' @rdname sign_values-methods
#' @export
#' @return sign_values, data.frame with the R2 or FDR p-value
#' for each of the selected variables
#' @examples
#'{
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo),
#'as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,
#'Distance = c('ScalarProduct','Euclidean','Euclidean'))
#'Selection<-VarSelection(Output,TaraOc,Crit='Rsquare',perc=0.95)
#'sign_values(Selection)
#'}
setMethod("sign_values", "VarSelection", function(x) (x@sign_values))
