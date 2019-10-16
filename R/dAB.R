#'@title dAB
#'@description Function to estimate differential abundance (if nCluster in
#' LinkData function is at least 2). The function uses a
#'non parametric kruskal-wallis test follow up by corrected p-values.
#'The function is robust since it doesn't assume normality on data distribution.
#'This function calculates the differential abundance (at OTU level)
#'betweeen all the communities data
#'It is only used when CLusters (enterotypes-like) is activated in LinkData function.
#'The function takes into account the
#'compositional nature of the OTUs dataset.
#'The differential expression is an alternative way to perform variable selection
#'@param x is an object of DistStatis Class.
#'@param Data should be the same imput list than in LinkData object.
#' If you integrated microbial communities and other types of data,
#'please be careful: choose only the microbial communities as input to dab object!!!!
#'@param adjust.methods character, correction method.
#' Choose one between:  c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY',
#'  'fdr', 'none').
#'@param threshold fixed pre-defined threshold value, which is
#' referred to as the level of significance.
#'@export dAB
#'@return Diferentialb: a list with selected OTUs and their p-values.
#'
#'
#'
#'@author Laura M Zingatetti
#'
#'@references
#'{
#'\enumerate{
#'\item
#'\item Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis.
#' Journal of the American statistical Association, 47(260), 583-621.
#'\item Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical
#' and powerful approach to multiple testing.
#'  Journal of the Royal Statistical Society Series B 57, 289–300.
#'\item Wright, S. P. (1992). Adjusted P-values for simultaneous inference.
#' Biometrics 48, 1005–1013. (Explains the adjusted P-value approach.)
#'}
#'}

#' @examples
#' {
#'data(Taraoceans)
#'pro.phylo <- Taraoceans$taxonomy[ ,'Phylum']
#'TaraOc<-list(Taraoceans$phychem,
#'as.data.frame(Taraoceans$pro.phylo),as.data.frame(Taraoceans$pro.NOGs))
#'TaraOc_1<-scale(TaraOc[[1]])
#'Normalization<-lapply(list(TaraOc[[2]],TaraOc[[3]]),
#'function(x){DataProcessing(x,Method='Compositional')})
#'colnames(Normalization[[1]])=pro.phylo
#'colnames(Normalization[[2]])=Taraoceans$GO
#'TaraOc<-list(TaraOc_1,Normalization[[1]],Normalization[[2]])
#'names(TaraOc)<-c('phychem','pro_phylo','pro_NOGs')
#'TaraOc<-lapply(TaraOc,as.data.frame)
#'Output<-LinkData(TaraOc,Scale =FALSE,Distance =
#'c('ScalarProduct','Euclidean','Euclidean'),nCluster=3)
#'dAB(Output,Data=list(TaraOc[[2]]))
#' }
#'
#' @name dAB
#' @rdname dAB-dAB
#' @importMethodsFrom  MultiAssayExperiment assays
#'
dAB <- function(x, Data, adjust.methods = "BH", threshold = 0.05) {
    
    # adjust.method <- match.arg(adjust.methods)
    if (is(Data, "MultiAssayExperiment")) {
        Data <- as.list(assays(Data))
    }
    
    if (!is(x, "DistStatis")) {
        stop("x must be a DistStatis class object")
    }
    ## Obtain x name for future update
    nameObject <- deparse(substitute(x))
    
    if ("fit.cluster" %in% colnames(compromise_coords(x)) == FALSE) {
        stop("to perform this test, you should compute clusters when executed LinkData function")
    }
    methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    if (adjust.methods %in% methods == FALSE) {
        stop("Please choose a valid correction method")
    }
    
    if (threshold < 0 | threshold > 1) {
        stop("threshold value refer to the significance level and
             should be higher than 0 and smaller than 1, usually 0.05")
    }
    
    if (!is(Data, "list")) {
        stop("Data should be a list with at least one element.
             List data.frame should be communities data")
    }
    WWW <- compromise_coords(x)
    
    Comparison <- lapply(Data, function(x) {
        apply(x, 2, function(y) {
            kruskal.test(y, WWW$fit.cluster)[[3]]
        })
    })
    p_adjust <- lapply(Comparison, function(x) {
        p.adjust(x, method = adjust.methods)
    })
    Selected <- list()
    for (i in seq_along(p_adjust)) {
        Selected[[i]] <- p_adjust[[i]][p_adjust[[i]] < threshold]
    }
    
    names(Selected) <- names(Data)
    
    return(DiferentialAb = Selected)
    
}
