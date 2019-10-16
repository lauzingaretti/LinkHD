#' @title Global-Plot
#' @description this function outputs a plot from a DistStatis object.
#' The plot shows the projection of  the all common observation
#' onto each subspace used at the integration step
#'
#' @param x DistStatis class object.
#'
#' @return plotted GlobalPlot/s of the component/s of the given DistStatis object.
#'
#' @author Laura M. Zingaretti
#'
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
#'GlobalPlot(Output) +
#'theme(panel.grid.major = element_blank(),
#'panel.grid.minor = element_blank(),
#'panel.background = element_blank(),
#'axis.line = element_line(colour = 'black'))
#'
#'
#' }
#' @exportMethod GlobalPlot
#' @docType methods
#' @usage \S4method{GlobalPlot}{DistStatis}(x)
#' @name GlobalPlot
#' @rdname DistStatis-GlobalPlot
#' @inheritParams DistStatis
#' @aliases GlobalPlot,DistStatis-method
#' @import ggplot2
#' @importFrom gridExtra grid.arrange



setGeneric("GlobalPlot", def = function(x) {
    standardGeneric("GlobalPlot")
})

setMethod(f = "GlobalPlot", c("DistStatis"), definition = function(x) {
    ## Check that is at element is available
    
    p <- list()
    for (i in seq_along(Trajectories(x))) {
        p[[i]] <- ggplot2::ggplot(Trajectories(x)[[i]][, seq_len(2)], aes(x = Trajectories(x)[[i]][, 1], y = Trajectories(x)[[i]][, 
            2])) + geom_point(size = 2, aes(colour = "#000099")) + ggtitle(as.character(unique(Trajectories(x)[[i]]$Studies))) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.line = element_line(colour = "black")) + theme(legend.position = "None") + labs(x = "Dim 1", 
            y = "Dim 2") + scale_y_continuous(labels = scales::number_format(accuracy = 1e-04)) + scale_x_continuous(labels = scales::number_format(accuracy = 1e-04)) + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 45, 
                hjust = 1))
    }
    if (length(p) == 2) {
        grid.arrange(grobs = p, ncol = round(length(p)/2), nrow = length(p), newpage = TRUE)
    } else {
        grid.arrange(grobs = p, ncol = round(length(p)/2), nrow = round(length(p)/2), newpage = TRUE)
        
    }
    
})
