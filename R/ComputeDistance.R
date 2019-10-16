#' @import vegan
ComputeDistance <- function(X, Distance, Z = diag(1, nrow = dim(X))) {
    
    # 'distance funtions' to use: 1. ScalarProduct 2. euclidean 3. manhatan 4. canberra 5. pearson 6.
    # pearsonabs 7. spearman 8.  spearmanabs 9. mahalanobis binary Data 10. jaccard 11. sm 12.  sokal&Sneath 13.
    # RogeryTanimoto 14. Dice 15. Hamman 16. Ochiai 17. Sokal&Sneath 18. Phi de Pearson 19.  Gower&Leendre
    # compositional data
    Distance <- switch(tolower(Distance), scalarproduct = "ScalarProduct", euclidean = "euclidean", bray = "bray", 
        manhattan = "manhattan", canberra = "canberra", pearson = "pearson", pearsonabs = "pearsonabs", spearman = "spearman", 
        spearmanabs = "spearmanabs", mahalanobis = "mahalanobis", jaccard = as.numeric(1), simple_maching = as.numeric(2), 
        sokal = as.numeric(3), roger_tanimoto = as.numeric(4), dice = as.numeric(5), hamman = as.numeric(6), 
        ochiai = as.numeric(7), sokal2 = as.numeric(8), phi_pearson = as.numeric(9), gower_legendre = as.numeric(10))
    
    if (is.null(Distance)) {
        Distance <- "ScalarProduct"
    }
    
    if (Distance == "ScalarProduct") {
        if (compbin(X) == "Binary Data Imput") {
            stop("Your Data are Binary and the Distance chosen is for continuous data ")
        }
        S <- (X) %*% t(X)
    }
    if (Distance == "mahalanobis") {
        if (compbin(X) == "Binary Data Imput") {
            stop("Your Data are Binary and the Distance chosen is for continuous data ")
        }
        S <- as.matrix(maha(X))
    }
    if (Distance == "euclidean" || Distance == "manhattan" || Distance == "canberra" || Distance == "pearson" || 
        Distance == "pearsonabs" || Distance == "spearman" || Distance == "spearmanabs") {
        if (compbin(X) == "Binary Data Imput") {
            stop("Your Data are Binary and the Distance chosen is for continuous data ")
        }
        S <- -1/2 * (Z %*% as.matrix(distan(X, meth.dis = tolower(Distance), diag = TRUE, upper = TRUE)) %*% 
            t(Z))
    }
    if (Distance == "bray") {
        if (compbin(X) == "Binary Data Imput") {
            S <- -1/2 * (Z %*% as.matrix(vegan::vegdist(X, method = Distance, binary = TRUE)) %*% t(Z))
        } else {
            S <- -1/2 * (Z %*% as.matrix(vegan::vegdist(X, method = Distance, binary = FALSE)) %*% t(Z))
            
        }
    }
    
    if (Distance == 1 || Distance == 2 || Distance == 3 || Distance == 4 || Distance == 5 || Distance == 6 || 
        Distance == 7 || Distance == 8 || Distance == 9 || Distance == 10) {
        if (compbin(X) == "Binary Data Imput") {
            S <- -1/2 * (Z %*% as.matrix(dist.binary(X, method = as.numeric(Distance), diag = TRUE, upper = TRUE)) %*% 
                t(Z))
        }
        if (compbin(X) == "Non-Binary Data Imput") {
            stop(paste("You chose the ", Distance, "for binary data and the table  is not binary-like"))
        }
    }
    
    return(S)
}
