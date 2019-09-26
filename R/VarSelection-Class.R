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
#' }
#'
#' @section Accesors:
#' \itemize{
#'\item  Variables return all the selected variables (and the frecuency of selection).
#'\item  Coordinates represent the coordenates (Betas coeeficients on LM) of the selected variables.
#'\item  VarTable dataframe indicating the table that each selected variable comes from.
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
#' showClass("VarSelection")
#' }
#'
#' @name VarSelection-class
#' @rdname VarSelection-Class
#' @export


setClass(Class="VarSelection",
         slots=c(
           Variables="character",
           Coordinates="data.frame",
           VarTable='data.frame'),
         prototype=prototype(
           Variables=c(),
           Coordinates=data.frame(),
           VarTable=data.frame())
)


