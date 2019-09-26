#' @import vegan
ComputeDistance<-function(X,Distance,Z=diag(1,nrow=dim(X))){

  #distance funtions to use
  #1. ScalarProduct
  #2. euclidean
  #3. manhatan
  #4. canberra
  #5. pearson
  #6. pearsonabs
  #7. spearman
  #8. spearmanabs
  #9. mahalanobis
  #binary Data
  #10. jaccard
  #11. sm
  #12. sokal&Sneath
  #13. RogeryTanimoto
  #14. Dice
  #15. Hamman
  #16. Ochiai
  #17. Sokal&Sneath
  #18. Phi de Pearson
  #19. Gower&Leendre
  #compositional data

  if (tolower(Distance)=="scalarproduct") Distance <- "ScalarProduct" else #(1)
    if (tolower(Distance)=="euclidean") Distance <- "euclidean" else #(2)
     if (tolower(Distance)=="bray") Distance <- "bray" else #(2bis)
      if (tolower(Distance)=="manhattan") Distance <- "manhattan" else #(3)
        if (tolower(Distance)=="canberra") Distance <- "canberra" else #(4)
          if (tolower(Distance)=="pearson") Distance <- "pearson" else #(5)
            if (tolower(Distance)=="pearsonabs") Distance <- "pearsonabs" else #(6)
              if (tolower(Distance)=="spearman") Distance <- "spearman" else #(7)
                if (tolower(Distance)=="spearmanabs") Distance <- "spearmanabs" else #(8)
                  if (tolower(Distance)=="mahalanobis") Distance <- "mahalanobis" else #(9)
                    if (tolower(Distance)=="jaccard") Distance <- as.numeric(1)  else #(10)
                      if (tolower(Distance)=="simple matching") Distance <- as.numeric(2)  else #(11)
                        if (tolower(Distance)=="sokal") Distance <- as.numeric(3) else #(12)
                          if (tolower(Distance)=="rogers&tanimoto") Distance <- as.numeric(4) else #(13)
                            if (tolower(Distance)=="dice") Distance <- as.numeric(5) else #(14)
                              if (tolower(Distance)=="Hamman") Distance <- as.numeric(6) else #(15)
                                if (tolower(Distance)=="Ochiai") Distance <- as.numeric(7) else #(16)
                                  if (tolower(Distance)=="sokal2") Distance <- as.numeric(8) else #(17)
                                    if (tolower(Distance)=="Phi-Pearson") Distance <- as.numeric(9) else #(18)
                                      if (tolower(Distance)=="Gower&Legendre") Distance <- as.numeric(10) else #(19)


                                        Distance <- "ScalarProduct"


                                      if(Distance=="ScalarProduct"){
                                        if (compbin(X)=="Binary Data Imput"){
                                          stop("Your Data are Binary and the Distance chosen is for continuous data ")

                                        }
                                        S <- (X)%*%t(X)
                                      }
                                      if(Distance=="mahalanobis"){
                                        if (compbin(X)=="Binary Data Imput"){
                                          stop("Your Data are Binary and the Distance chosen is for continuous data ")

                                        }
                                        S <- as.matrix(maha(X))
                                      }
                                      if(Distance=="euclidean" || Distance=="manhattan" ||Distance=="canberra"||Distance=="pearson" ||
                                         Distance=="pearsonabs"|| Distance=="spearman"|| Distance=="spearmanabs"){
                                        if (compbin(X)=="Binary Data Imput"){
                                          stop("Your Data are Binary and the Distance chosen is for continuous data ")
                                        }
                                        S <- -1/2*(Z%*%as.matrix(distan(X, meth.dis = Distance, diag = TRUE, upper = TRUE))%*%t(Z))
                                      }
                                      if(Distance=="bray"){
                                        if (compbin(X)=="Binary Data Imput"){
                                        S <- -1/2*(Z%*%as.matrix(vegan::vegdist(X, method=Distance,binary=TRUE))%*%t(Z))
                                          }else{
                                        S <- -1/2*(Z%*%as.matrix(vegan::vegdist(X, method=Distance,binary=FALSE))%*%t(Z))

                                        }
                                      }



                                      if(Distance==1|| Distance==2 ||Distance==3||Distance==4 ||
                                         Distance==5|| Distance==6|| Distance==7|| Distance==8 ||
                                         Distance==9 || Distance==10){
                                        if(compbin(X)=="Binary Data Imput"){
                                          S<- -1/2*(Z%*%as.matrix(dist.binary(X, method = as.numeric(Distance), diag = TRUE, upper = TRUE))%*%t(Z))
                                        }
                                        if(compbin(X)=="Non-Binary Data Imput"){
                                          stop(paste("You chose the ", Distance,"for binary data and the table  is not binary-like"))
                                        }
                                      }

                                      return(S)
}
