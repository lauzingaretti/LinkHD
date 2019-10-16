#' @title LinkData: multiple heterogeneous dataset integration
#' @description  Integrating multiple Heterogeneous Datasets
#' stored into a list. This function makes Statis using Distances options.
#' Statis is part of the PCA family and is based on singular value decomposition
#' (SVD) and the generalized singular value decomposition (GSVD) of a matrix.
#' This methodology aims to analyze several data sets of
#' variables that were collected on the same set of observations.
#' Originally, the comparisons were drawn from the compute
#' of the scalar product between the different tables.
#' In our approach, the condition is relaxing
#' allowing  the incorporation of different distances.
#'@param Data should be a list of dataframes or ExpressionSet data
#' with the same length of the number of tables to be integrate.
#' In each dataframe, the Observations (common elements on Statis)
#' should be in rows and the variables should be in columns.
#' Data also might be a MultiAssayExperiment object
#' from MultiAssayExperiment package, a software for
#' multi-omics experiments integration in Bioconductor.
#'@param Distance Vector indicating which distance (including scalar product)
#' should be applied to each study. If is missing,
#' the scalar product is used. The vector lenght must be equal
#' to the length of Data. Distance options: ScalarProduct, euclidean,
#' manhattan, canberra, pearson, pearsonabs, spearman,
#' spearmanabs, mahalanobis, BrayCurtis distance (please, use option Bray).
#' For binary data, the distance can be jaccard,
#' simple_matching, sokal_Sneath, Roger_Tanimoto, Dice,
#' Hamman, Ochiai, Phi_Pearson, 'Gower&Legendre.
#' Note that, use pre-processing option as compositional and
#' Euclidean is the same than use Aitchison distance for compositional data.
#'@param Center Logical. If TRUE, the data frame
#'is centered  by the mean. By default is FALSE.
#' If you have tables with different characteristics (continous phenotypes, frecuencies,
#'compositional data), we strongly recomendate normalize
#' datasets as a previous step through DataProcessing option.
#'@param Scale A logical value indicating whether the column vectors should be
#'standardized by the rows weight, by default is FALSE.
#'Note that all data into the list will be scaled.
#'If you don't need normalizing all data, you
#'could set this parameter as False and perform the normalization step
#'externally by using DataProcessing function.
#'If you have tables with different characteristics (continous phenotypes, frecuencies,
#'compositional data), we strongly recomendate normalize datasets
#' as a previous step through DataProcessing option.
#'@param CorrelVector Logical. If TRUE (default), the RV matrix is
#'computed using vectorial correlation, else
#'the Hilbert-Smith distance is used.
#'@param nCluster this variable indicates if common
#'elements on the dataset should be grouped (by default is zero, i.e. no-cluster).
#'@param cl_method categorical (pam or kmeans). pam is a robust
#' version of classical kmeans algorithm.
#'@return \item{LinkData}{DistStatis class object with the
#'corresponding completed slots according to the given model}
#'@export LinkData
#'@author Laura M Zingatetti
#'
#' @references
#' \enumerate{
#'  \item Escoufier, Y. (1976). Operateur associe a un tableau de donnees.
#'   Annales de laInsee, 22-23, 165-178.
#'  \item Escoufier, Y. (1987). The duality diagram: a means
#'   for better practical applications. En P. Legendre & L. Legendre (Eds.),
#'   Developments in Numerical Ecology, pp. 139-156,
#'   NATO Advanced Institute, Serie G. Berlin: Springer.
#'  \item L'Hermier des Plantes, H. (1976). Structuration des
#'  Tableaux a Trois Indices de la Statistique. [These de Troisieme Cycle].
#'  University of Montpellier, France.
#'}
#'
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
#' }
#'
#' @name LinkData
#' @rdname LinkData-LinkData
#' @import methods
#' @importFrom methods as is
#' @import scales
#' @importFrom cluster pam
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importMethodsFrom  MultiAssayExperiment assays

LinkData <- function(Data, Distance = c(), Center = FALSE, Scale = FALSE, CorrelVector = TRUE, nCluster = 0,
    cl_method = "pam") {

    ## Auxiliary functions Function for implementation the Statis Method with Distance options. The Data must be
    ## a list of data.frames or ExpressionSet check if elements of list are ESet data.
    if (!is(Data, "list") & !is(Data, "MultiAssayExperiment")) {
        stop("your dataset should be stored in a list containing data.frame or ExpressionSet object or a MultiAssayExperiment object ")
    }

    if (is(Data, "list")) {
        if (unique(lapply(Data, class)) == "ExpressionSet") {
            Data <- lapply(Data, function(x) {
                as(x, "data.frame")
            })
            Data <- lapply(Data, function(x) {
                t(x)
            })
            Data <- lapply(Data, as.data.frame)
        }
    }

    if (is(Data, "MultiAssayExperiment")) {
        Data <- assays(Data)
    }

    if (length(Data) < 2) {
        stop("you need a list of at least two data.frames")
    }

    if (nCluster < 0) {
        nCluster = 0
        message("warning: nCluster should be a positive integer")
    }



    # checking if data.frame have rows and columns names
    columns <- lapply(Data, colnames)
    rows <- lapply(Data, rownames)
    Names_Col <- lapply(columns, is.null)
    Names_Row <- lapply(rows, is.null)

    if (any(Names_Col == TRUE)) {
        stop("All your data.frame should include columns names")
    }

    if (any(Names_Row == TRUE)) {
        stop("All your data.frame should include rownames names")
    }

    if (is.null(names(Data))) {
        names(Data) <- seq_along(Data)
    }



    # rownames should be the same in all tables
    if (all(apply((x <- vapply(Data, rownames, character(nrow(Data[[1]])))), 2, function(y) +identical(y, x[,
        1]))) == FALSE) {
        stop("tables should be the same observations")
    }

    # Normalice Data by frec or by sd else use Normalize

    if (Scale == FALSE & Center == FALSE) {
        X <- lapply(Data, as.matrix)
        names(X) <- names(Data)
    } else {
        X <- Normalize(Data, scale = Scale, center = Center)
        names(X) <- names(Data)
    }



    # In--> number of observations.
    In <- nrow(X[[1]])
    # Create a Diagonal matrix with 1/n
    D <- diag(rep(1, In))/In

    # Z matrix is the centering matrix used to calculate distances, see Abdi (2007) see Abdi 2007 (Multiple
    # distances matrix)
    M1 <- matrix(rep(1, In)/In, ncol = In, nrow = 1)
    Ones <- matrix(rep(1, In), ncol = 1, nrow = In)
    Id <- diag(1, nrow = In, ncol = In)
    Z <- Id - Ones %*% M1


    # get data names put X row/col names using Data

    for (i in seq_along(Data)) {
        colnames(X[[i]]) <- colnames(Data[[i]])
        rownames(X[[i]]) <- rownames(Data[[i]])
    }


    # Calculate Distance for each table If not specify Distance or length(Distance)!=K, the scalar product is
    # used This list name (S) follows Abdi (2007)
    S = list()



    if (length(Distance) != length(Data)) {
        S <- lapply(X, ScalarProduct)
        message("The calculations were performed using the scalar product between the tables")
        Distance <- "scalar-product"
    }

    if (length(Distance) == length(Data)) {
        S <- lapply(seq_along(Distance), function(i) {
            ComputeDistance(X[[i]], Distance[i], Z)
        })
    }

    names(S) <- names(Data)

    ########################################################################## Vectorial correlation (SW)############################
    SW <- list()

    for (k in seq_along(S)) {
        wk <- as.matrix(S[[k]]) %*% D
        wk <- t(t(wk) %*% D)
        wk <- wk %*% t(wk)
        SW[[k]] <- wk
    }

    if (CorrelVector == TRUE) {
        sep <- matrix(unlist(SW), In * In, length(S))
        RV <- t(sep) %*% sep
        ak <- sqrt(diag(RV))
        RV <- sweep(RV, 1, ak, "/")
        RV <- sweep(RV, 2, ak, "/")
        dimnames(RV) <- list(names(S), names(S))
    } else {
        sep <- matrix(unlist(SW), In * In, length(S))
        RV <- t(sep) %*% sep
        dimnames(RV) <- list(names(S), names(S))
    }


    ######################################################################### INTER-STRUCTURE################################

    SvdRV <- svd(RV)
    # percentage of inertia explained by the dimensions
    InertiaExpRV <- c(((SvdRV$d)/sum(diag(SvdRV$d))) * 100)
    InertiaExpRV <- data.frame(InertiaExpRV)
    # inertia accumulated
    CumInertiaExpRV <- cumsum(InertiaExpRV)
    InertiaExpRV <- as.data.frame(cbind(SvdRV$d, InertiaExpRV, CumInertiaExpRV))
    colnames(InertiaExpRV) <- c("Value", "Inertia(%)", "Cumulative Inertia (%)")
    Dim <- c(paste("Dim", seq_along(S)))
    rownames(InertiaExpRV) <- Dim

    ### the following code is just to plot the Euclidean image of each table###

    cc <- SvdRV$u %*% sqrt(diag(SvdRV$d))
    if (any(cc[, 1] < 0))
        cc[, 1] <- -cc[, 1]
    ImSt <- as.data.frame(cc[, seq_len(2)])
    rownames(ImSt) <- names(S)
    colnames(ImSt) <- c("Dim1", "Dim2")

    ###### Calculate the cosine between studies#########
    A <- cc[, seq_len(2)] %*% t(cc[, seq_len(2)])
    M <- diag(1/sqrt(diag(A)))
    Cosin <- M %*% A %*% M
    Angulos <- acos(as.dist(Cosin))
    SqCos <- Cosin^2
    SqCos <- as.data.frame(SqCos)
    colnames(SqCos) <- names(S)
    rownames(SqCos) <- names(S)

    ##################################################### INTER-STRUCTURE:Compromise #############

    # The first eigenvector give the optimal weights to compute the compromise matrix.  Practically, the optimal
    # weights can be obtained by re-scaling these values such that their sum is equal to one. So the weights are
    # obtained by dividing each element of p1 by the sum of the elements of p1.


    sc <- sum(sqrt(diag(RV)))
    if (any(SvdRV$u[, 1] < 0)) {
        SvdRV$u[, 1] <- -SvdRV$u[, 1]
    }

    pit <- vapply(seq_len(nrow(RV)), function(i) {
        (sc/sqrt(SvdRV$d[1])) * SvdRV$u[i, 1]
    }, numeric(1))
    alphas <- pit/sum(pit)

    # WW is just the compromise matrix
    WW <- matrix(0, nrow = nrow(S[[1]]), ncol = ncol(S[[1]]))

    for (i in seq_along(S)) {
        WW <- WW + alphas[i] * S[[i]]
    }

    ### www data.frame compromise
    WWW <- as.data.frame(WW)
    rownames(WWW) <- as.matrix(rownames(Data[[1]]))
    colnames(WWW) <- as.matrix(rownames(Data[[1]]))


    ##################################################################################### Singular value Decomposition of Compromise matrix (WWW)################

    SvdComp <- svd(WW)
    InerComp <- c((SvdComp$d/sum(diag(SvdComp$d))) * 100)
    # InerComp porcentaje de inercia explicado por las dimensiones
    InerComp <- cbind(SvdComp$d, InerComp, cumsum(InerComp))
    Dim <- c(paste("Dim", seq_len(nrow(WW))))
    rownames(InerComp) <- Dim
    colnames(InerComp) <- c("Values", "Inertia", "Cumulative Inertia")


    ########### Projection of compromise' observations###########

    AP <- WW %*% SvdComp$u %*% diag((1/sqrt(SvdComp$d)))
    # AP for plot
    ProjObs <- as.data.frame(AP[, seq_len(min(4, ncol(AP)))])
    rownames(ProjObs) <- rownames(Data[[1]])
    A <- paste("Dim", seq_len(min(4, ncol(AP))))
    colnames(ProjObs) <- A

    ##################################################### Rendering Quality of the Observations (RQO,RQI)###
    RQITot <- as.data.frame((apply(AP[, seq_len(2)]^2, 1, sum)/apply(AP^2, 1, sum)) * 100)
    rownames(RQITot) <- rownames(S[[1]])
    colnames(RQITot) <- "RQI(%)"


    ###################################################### CLUSTER######################################
    if (nCluster > 0) {
        if (floor(nCluster) != nCluster) {
            stop("The number of cluster must to be a positive integer")
        }
        if (nCluster == 1) {
            nCluster == 2
            message("You must to ask by at least two clusters")
        }
        if (nCluster > nrow(ProjObs)) {
            stop("the number of clusters should be smaller than the number of observations")
        }

        # k-means over compromise matrix
        if (cl_method == "pam") {
            fit <- pam(ProjObs[, seq_len(2)], nCluster)
            Means_clusters = aggregate(ProjObs, by = list(fit$clustering), FUN = mean)
            # append cluster assignment
            ProjObs <- data.frame(ProjObs, fit$cluster)
        }
        if (cl_method == "km") {
            fit <- kmeans(ProjObs[, seq_len(2)], nCluster)
            Means_clusters = aggregate(ProjObs, by = list(fit$cl), FUN = mean)
            # append cluster assignment
            ProjObs <- data.frame(ProjObs, fit$cl)

        }
    }

    #################### INTRA-STRUCTURE################### The intrastructure step is a projection of the rows of each table of
    #################### the series into the multidimensional space of the compromise analysis.


    ### projection of each table on Compromise configuration

    Studies <- unique(names(S))
    Observations <- rownames(X[[1]])
    TableProjections <- lapply(seq_along(S), function(i) {
        TabProj(S[[i]], SvdComp, Studies[i], Observations)
    })
    names(TableProjections) <- names(S)
    # return projections

    ## Return of different slots
    .Object <- new("DistStatis", Inertia.RV = InertiaExpRV, RV = RV, Euclid.Im = ImSt, Inertia.comp = InerComp,
        Compromise.Coords = ProjObs, Compromise.Matrix = WWW, RQO = RQITot, TableProjections = TableProjections)

    validObject(.Object)
    return(.Object)
}








