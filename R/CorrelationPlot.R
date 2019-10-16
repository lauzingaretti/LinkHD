#' @title Correlation-Plot
#' @description Plot a \code{CorrelationPlot} of a DistStatis object
#'
#' @param x an object from DistSatis class.
#' @param ... additional parameters from ggplot2 library
#'
#' @return correlation plot between tables from a  DistStatis object.
#'
#' @author Laura M. Zingaretti
#'
#' @examples
#' {
#'
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
#' CorrelationPlot(Output) +
#' theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#' panel.background = element_blank(),
#' axis.line = element_line(colour = 'black'))
#'
#'
#' }
#' @exportMethod CorrelationPlot
#' @inheritParams DistStatis
#' @docType methods
#' @usage \S4method{CorrelationPlot}{DistStatis}(x,...)
#' @name CorrelationPlot
#' @rdname DistStatis-CorrelationPlot
#' @aliases CorrelationPlot,DistStatis-method
#' @import ggplot2
#' @import ggpubr
#' @importFrom reshape2 melt
#' @import emmeans




setGeneric("CorrelationPlot", def = function(x, ...) {
    standardGeneric("CorrelationPlot")
})




setMethod(f = "CorrelationPlot", c("DistStatis"), definition = function(x, ...) {
    ## Check that is at element is available
    
    Z <- correl(x)
    diag(Z) <- 1
    m <- Z
    melted_cormat <- melt(m, preserve.na = FALSE)
    colnames(melted_cormat) <- c("x", "y", "RV")
    
    m <- melted_cormat
    
    p = ggplot2::ggplot(na.omit(m), aes(m[, 1], m[, 2]))
    
    p = p + geom_point(aes(size = abs(m[, 3]) * 1), color = "white")
    p = p + geom_point(aes(size = abs(m[, 3]), color = m[, 3]))
    p = p + scale_size_continuous(range = c(2, 10), name = "RV") + guides(size = FALSE) + scale_colour_gradient2()
    p + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + labs(x = "", 
        y = "") + theme(legend.title = element_blank())
    
})
