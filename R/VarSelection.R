#'@title Variable Selection
#'@description Function to do variable selection using a Regression Biplot methodology.
#' This function calculates the regression biplot on the compromise matrix. Biplot can be
#' understood as the decomposition of a target matrix ($Y=XB$). Here, $Y$ is the matrix containing all
#'  variables taken into account in the analisis,$X$ is the matrix containing the explaining variables, i.e., the coordinates of compromise matrix
#'  and finally, $B$ are the regression coefficients to be estimated. Then, the method is interpreted as
#' a general linear regression into the $X$ matrix (${Y_hat}=X(X'X)^(-1)X'Y$) and the matrix $X(X'X)^(-1)X'$ is the projection matrix
#' onto the compromise configuration. We  use a classical linear model to obtain the regressors coefficients, however
#' the model could be extended and alternatives methods are able to use.
#' The quality of the regression biplot is measured using the proportion of explained variance
#' by each regression (adjusted r squared coefficient).
#'
#'@param x is an object of DistStatis Class.
#'@param Data should be a list of data.frame or ExpressionSet data with the same length of the number of tables to be integrate.
#'  In each dataframe, the Observations (common elements on Statis) should be in rows and the variables should be in columns. Data are the same
#'  data used to obtained the compromise configuration.It also can be a MultissayExperiment object, please check help of LinkData function and the package vignette.
#'@param intercept Logical. If is TRUE, the models with intercept are computed, else the intercept is zero.
#'@param model character. 'LM' for classical lm model. We've planned to implemening alternative models in the future.
#'@param Crit Character indicating the variable selection criteria.You could chose 'Rsquare' or 'p-val'.
#'@param perc The value of percentil that indicate how much data than are selected.
#'@param nDims Numeric that indicates the number of dimensions to use for do the model. Default is 2.
#'@param Normalize Logical. If is TRUE, the response variable in each model is normalized.
#'@return a \item{VarSelection}{VarSelection class with the
#'   corresponding completed slots
#'   according to the given model}
#'
#'
#'
#'@author Laura M Zingatetti
#'
#'@references
#'{
#'\enumerate{
#'\item
#'\item Gabriel, K. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika 58(3), 453--467.
#'\item Gower, J. & Hand, D. (1996). Biplots, Monographs on statistics and applied probability. 54. London: Chapman and Hall., 277 pp.
#'\item Greenacre, M. J. (2010). Biplots in practice. Fundacion BBVA.
#'}
#'}

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
#'}


#'@name VarSelection
#'@rdname VarSelection
#'@export
#'@import methods
#'@importMethodsFrom  MultiAssayExperiment assays
# x is an object comes from DistStatis class, intercept idicates if the model would be runned with/without
# intercept crit is the model selection criteria used to variable selection Rsquare/ p-value perc proportion
# of observations to select nDims N how many dims do you want to use as explained variables in the model?



VarSelection <- function(x, Data, intercept = FALSE, model = "LM", Crit = "Rsquare", perc = 0.9, nDims = 2,
    Normalize = FALSE) {
   #fix the problem with categories Dec 11, 2019

    ## Obtain x name for future update
    nameObject <- deparse(substitute(x))
    if (!is(x, "DistStatis")) {
        stop("x should be an object from DistStatis-class")
    }
    if (is(Data, "MultiassayExperiment")) {
        Data <- assays(Data)
    }

    WWW <- Compromise_matrix(x)

    ## Check parameters
    stopifnot(Crit[1] %in% c("Rsquare", "p-val"))
    if (is.null(perc) == TRUE) {
        perc = 0.95
    }
    if (Crit != "Rsquare" & Crit != "p-val") {
        stop("Non-valid selection criteria")
    }

    if (perc < 0 || perc > 1) {
        stop("Error in perc: Invalid percentil choice. Percentil have to be between 0 and 1")
    }
    if (nDims > ncol(WWW)) {
        nDims <- 2
        message("Arg nDims should be smaller than the number of observations")
    }
    if (nDims <= 1) {
        nDims <- 2
        message("Arg nDims must to be upper than 1")
    }
    if (model != "LM") {
        stop("We only had been implemented the standard Linear Model")
    }

    nDims <- round(nDims)

    SvdComp <- svd(as.matrix(WWW))

    AP <- as.matrix(WWW) %*% SvdComp$u %*% diag((1/sqrt(SvdComp$d)))
    # AP for plot
    ProjObs <- as.data.frame(AP[, seq_len(nDims)])
    rownames(ProjObs) <- rownames(WWW)
    colnames(ProjObs) <- paste("Dim", seq_len(nDims))


    Mat <- as.matrix(ProjObs)
    if (is(Data, "MultiAssayExperiment")) {
        Datos <- Data@ExperimentList@listData
    } else {
        Datos <- Data
    }



    if (all(apply((x <- vapply(Datos, rownames, character(nrow(Datos[[1]])))), 2, function(y) +identical(y,
        rownames(Datos[[1]])))) == FALSE) {
        stop("tables should be the same observations")
    }

    if (Normalize == TRUE) {
        Datos <- lapply(Datos, function(x) {
            scale(x, center = TRUE, scale = TRUE)
        })
    } else {
        Datos <- Datos
    }
    if (model == "LM") {
        s <- lapply(Datos, function(m) {
            apply(m, 2, function(y) LinModel(Mat, y, intercept = intercept))
        })
    }

    sn <- lapply(seq_along(s), function(i) {
        s[[i]][!vapply(s[[i]], is.null, logical(1))]
    })
    names(sn) <- names(Datos)



    # dataset contains all models (coordinates and R2, pval to each table)
    dataset <- lapply(sn, function(x) {
        as.data.frame(matrix(unlist(x), nrow = 4, byrow = FALSE))
    })

    # p val correction

    # dataset<-lapply(dataset,function(x){x[2,]<-x[2,]/ncol(x)})
    for (i in seq_along(dataset)) {
        colnames(dataset[[i]]) <- names(s[[i]][!vapply(s[[i]], is.null, logical(1))])
        dataset[[i]] <- rbind(dataset[[i]], rep(names(dataset)[i], ncol(dataset[[i]])))
        rownames(dataset[[i]]) <- c("R2", "pv", "x", "y", "Table")
        if (i == 1) {
            all_dat <- dataset[[i]]
        } else {
            all_dat <- cbind(all_dat, dataset[[i]])
        }
    }

    nvar <- sum(unlist(lapply(Datos, ncol)))
    # number of variable to select
    Tosel <- floor(nvar - nvar * perc)



    if (Crit == "Rsquare") {
        values <- unlist(lapply(dataset, function(x) {
            x[1, ]
        }))
        Nam<-names(values)
        values<-as.numeric(values)
        names(values)<-Nam
        R2 <- sort(values, decreasing = TRUE)
        Seleccionados <- R2[seq_len(Tosel)]
        indices <- which(values %in% Seleccionados)
        # indices<-order(values,decreasing = TRUE)[seq_len(Tosel)]
                if(any(Seleccionados<0)){
          Seleccionados[Seleccionados<0]=0
        }
    }

    if (Crit == "p-val") {
        values <- p.adjust(unlist(lapply(dataset, function(x) {
            x[2, ]
        })), method = "fdr")
        Nam<-names(values)
        values<-as.numeric(values)
        names(values)<-Nam
        # we could considerer add other methods
        pval <- sort(values)
        Seleccionados <- pval[seq_len(Tosel)]
        indices <- order(values)[seq_len(Tosel)]


    }

    ### look for the variables into each table
    Nam <- names(Seleccionados)
    index <- lapply(c(seq_along(dataset)), function(i) {
        startsWith(Nam, names(dataset)[i])
    })


    variables <- c()
    for (i in seq_along(dataset)) {
        variables <- c(variables, unlist(lapply(Nam[index[[i]]], function(m) {
            strsplit(m, paste0(names(dataset)[i], "."))[[1]][2]
        })))
    }

    if (Crit == "Rsquare") {
        Val <- all_dat[1, indices]
    }
    if (Crit == "p-val") {
        Val <- all_dat[2, indices]
    }

    .Object <- new("VarSelection", Variables = variables, Coordinates = all_dat[c(3:4), indices], VarTable = all_dat[5,
        indices], sign_values = Val)
    validObject(.Object)
    return(.Object)

}
