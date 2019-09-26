#' @title Compromise-Plot
#' @description Plot a \code{CompromisePlot} of a DiStatis object
#'
#'
#' @param x DistStatis class object.
#' @param x_lab a character indicating x_label. Default is x.
#' @param y_lab a character indicating y_label. Default is y.
#' @param Name a character indicating plot title.
#' @param pchPoints pch for points in scatter plot.
#' @param colObs is a  character indicating the color for the observations. By Default is the QR (indicating the Quality of Representation of observations)
#' @param ... additional parameters from ggplot2 library
#'
#' @return plotted CompromisePlot/s of the component/s of the given DistStatis object.
#'
#' @author Laura M. Zingaretti
#'
#' @examples
#' {
#' \dontrun{
#' #'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,"Phylum"]
#'TaraOc<-list(Taraoceans$phychem,as.data.frame(Taraoceans$pro.phylo),as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),function(x){DataProcessing(x,Method="Compositional")})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c("phychem","pro_phylo","pro_NOGs")
#'TaraOc<-lapply(TaraOc,as.data.frame)
# Output<-LinkData(TaraOc,Scale =FALSE,Distance = c("ScalarProduct","Euclidean","Euclidean"))
#'CompromisePlot(Output) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#'                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
#'
#' }
#' }
#' @exportMethod CompromisePlot
#' @docType methods
#' @usage \S4method{CompromisePlot}{DistStatis}(x,x_lab=NULL, y_lab=NULL,
#' Name=NULL, pchPoints=2, colObs=NULL,...)
#' @name CompromisePlot
#' @rdname DistStatis-CompromisePlot
#' @aliases CompromisePlot,DistStatis-method
#' @import ggplot2
#'
setGeneric("CompromisePlot",def=function(x,x_lab=NULL, y_lab=NULL,
          Name=NULL, pchPoints=2,colObs=NULL,...){standardGeneric("CompromisePlot")})



setMethod(f="CompromisePlot", signature="DistStatis", definition=function(x,x_lab=NULL, y_lab=NULL,
         Name=NULL, pchPoints=2,colObs=NULL,...){
  ##Check that is at element is available
  if( class(x)!="DistStatis"){
    stop("CompromisePlot requires a DistStatis object")
  }

    PARACCIND<-x@Compromise.Coords
    colnames(PARACCIND)<-paste0("Dim",seq(1:ncol(PARACCIND)))
    if(is.null(x_lab)){
    x_lab= paste0("Dim 1( ",round(x@Inertia.comp[1,2],2)," %)")
    }
    if(is.null(y_lab)){
      y_lab= paste0("Dim 2( ",round(x@Inertia.comp[2,2],2)," %)")
    }

    if(is.null(colObs)){
      colObs<-x@RQO$`RQI(%)`
      QR<-colObs
    }else{
     QR<-colObs
    }

    ggplot2::ggplot(PARACCIND, aes(x=Dim1, y=Dim2)) + geom_point(size=pchPoints,aes(color=QR)) +
      labs(x = x_lab,y=y_lab) + ggtitle(Name)


})
