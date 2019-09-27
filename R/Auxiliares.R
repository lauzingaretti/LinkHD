#File Created by Laura M. Zingaretti
#Adapted from previous kimod package function
#Date April, 24 2019
#' @import stats

#geometric mean
g_mean <- function(Data) {
  if (any(na.omit(Data == 0))){
    0
  }
  else {exp(mean(log(c(Data)[which(is.finite(unlist(Data)) & Data > 0)])))
}
}


#Center_LR (Centered Log Ratio)
centerLR<-function (Data, b = exp(1))
{
  if (!is(Data,"data.frame") & !is(Data,"matrix")){
    stop("The dataset should be a matrix or data.frame")
  }
  #just if data.frame or matrix has only one column
  if (dim(Data)[2] == 1) {
    result <- list(Data.clr = Data, geom_m = rep(1, dim(Data)[1]))
  }

  else {
    gm <- apply(Data, 1, g_mean)
    Data.clr <- log(Data/gm, b)
    result <- list(Data.clr = Data.clr, geom = gm)
  }
    return(result)
}


#Normalize data (scale/center, both or one)

Normalize <-function(X,scale=TRUE,center=TRUE){
  if(!is(X,"list")){
    stop("Error:Object of type 'list' expected")
  }

  f2 <- function(v){
    if(!is.numeric(v)){
      stop("Error:Object of type 'numeric' expected")
    }
    row.w <- rep(1, length(v))/length(v)
    sqrt(sum(v * v * row.w)/sum(row.w))
  }

  # list data
  XC=list()
  # list Aux
  M=list()

  for (i in 1: length(X)){
    if (center==TRUE && scale==TRUE){
      XC[[i]] <- matrix(scale(X[[i]],center=TRUE,scale=FALSE),nrow=nrow(X[[i]]))
      M[[i]] <- apply(XC[[i]],2,f2)
      M[[i]][M[[i]]< 1e-08] <- 1
      XC[[i]]<- sweep(XC[[i]], 2, M[[i]], "/")
      XC[[i]] <- XC[[i]]*(as.numeric(1/sqrt(sum((XC[[i]]%*%t(XC[[i]]))^2))))
    }

    if (center==FALSE && scale==FALSE){
      XC[[i]] <- X[[i]]
    }
    if (center==TRUE && scale==FALSE){
      XC[[i]] =matrix(scale(X[[i]],scale=FALSE),nrow=nrow(X[[i]]))
    }
    if (center==FALSE && scale==TRUE){
      M[[i]] <- apply(X[[i]],2,f2)
      M[[i]][M[[i]]< 1e-08] <- 1
      XC[[i]]<- sweep(XC[[i]], 2, M[[i]], "/")
      XC[[i]]=XC[[i]]*(as.numeric(1/sqrt(sum((XC[[i]]%*%t(XC[[i]]))^2))))

}
  }

  return(XC)
}

#calculate Scalar product between configurations
ScalarProduct <-
  function(X,Row=TRUE){
    if(!is.matrix(X))
    {
      X <- as.matrix(X)
    }
    clases<-apply(X,2,class)
    if(any(clases!="numeric") && any(clases!="integer")){
      stop("The Scalar Product only should be used with numerical data")
    }
    if (Row==TRUE){
      Y <- X%*%t(X)
    }
    else{
      Y <- t(X)%*%X
    }
    return(Y)
  }

#Should be data try as Frecuencies?
#Then Data should be positives.
cia<-function(df){
  df <- as.data.frame(df)
  if (!is(df,"data.frame"))
    stop("data.frame expected")
  if (any(df < 0))
    stop("negative entries in table and your data should be frecuencies")
  if ((N <- sum(df)) == 0)
    stop("all frequencies are zero")
  df <- df/N

  row.w <- apply(df, 1, sum)
  col.w <- apply(df, 2, sum)
  df <- df/row.w
  df <- sweep(df, 2, col.w, "/") - 1
  if (any(is.na(df))) {
    fun1 <- function(x) {
      if (is.na(x))
        return(0)
      else return(x)
    }
    df <- apply(df, c(1, 2), fun1)
    df <- data.frame(df)
  }
  df
  return(df)
}


#########Distances computation#########

#functions to compute metrics

distan <- function(mat=NULL, meth.dis="euclidean",diag=FALSE,upper=FALSE) {

MEASURE <- c("euclidean","manhattan","canberra" ,"pearson", "pearsonabs", "spearman", "spearmanabs", "mahalanobis")
mea <- pmatch(meth.dis, MEASURE)

  if (is.na(mea)) stop("Error :Unknown Metric.")
  if (mea==1) DIS <- as.matrix(dist(mat, method="euclidean",diag,upper))
  if (mea==2) DIS <- as.matrix(dist(mat, method="manhattan",diag,upper))
  if (mea==3) DIS <- as.matrix(dist(mat, method="canberra",diag,upper))
  if (mea==4) DIS <- (1-cor(t(mat), method="pearson"))/2
  if (mea==5) DIS <- 1-abs(cor(t(mat), method="pearson"))
  if (mea==6) DIS <- (1-cor(t(mat), method="spearman"))/2
  if (mea==7) DIS <- 1-abs(cor(t(mat), method="spearman"))
  if (mea==8) DIS <- maha(mat)
  attr(DIS,"Metric")<-meth.dis

  return(DIS)
}

dist.binary <-
  function (df=NULL, method = NULL, diag = FALSE, upper = FALSE) {
    METHODS <- c("JACCARD S3", "SOCKAL & MICHENER S4", "SOCKAL & SNEATH S5",
                 "ROGERS & TANIMOTO S6", "CZEKANOWSKI S7", "GOWER & LEGENDRE S9", "OCHIAI S12", "SOKAL & SNEATH S13",
                 "Phi of PEARSON S14", "GOWER & LEGENDRE S2")
    if (!(inherits(df, "data.frame") | inherits(df, "matrix")))
      stop("df is not a data.frame or a matrix")
    df <- as.mat
    compbin<-function(X=NULL){
      if(!is(X,"data.frame") && !is(X,"matrix")){
        stop("invalid class of Object")
      }
      X<-na.omit(X)
      Unicos<-unlist(apply(X,1,unique))
      A<-c()
      if(max(Unicos)!=1){
        A<-c(A,1)
      }

      if(min(Unicos)!=0){
        A<-c(A,2)
      }

      M<-sort(unique(c(Unicos)))
      if(M[1]!=0 || M[2]!=1){
        A<-c(A,3)
      }

      if(length(A)==0){
        return("Binary Data Imput")
      }else{

        return("Non-Binary Data Imput")
      }
    }

    rix(df)
    if(!is.numeric(df))
      stop("df must contain  numeric values")
    if (any(df < 0))
      stop("non negative value expected in df")
    nlig <- nrow(df)
    d.names <- row.names(df)
    if(is.null(d.names))
      d.names <- 1:nlig
    nlig <- nrow(df)
    df <- as.matrix(1 * (df > 0))
    if (is.null(method)) {
      cat("1 = JACCARD index (1901) S3 coefficient of GOWER & LEGENDRE\n")
      cat("s1 = a/(a+b+c) --> d = sqrt(1 - s)\n")
      cat("2 = SOCKAL & MICHENER index (1958) S4 coefficient of GOWER & LEGENDRE \n")
      cat("s2 = (a+d)/(a+b+c+d) --> d = sqrt(1 - s)\n")
      cat("3 = SOCKAL & SNEATH(1963) S5 coefficient of GOWER & LEGENDRE\n")
      cat("s3 = a/(a+2(b+c)) --> d = sqrt(1 - s)\n")
      cat("4 = ROGERS & TANIMOTO (1960) S6 coefficient of GOWER & LEGENDRE\n")
      cat("s4 = (a+d)/(a+2(b+c)+d) --> d = sqrt(1 - s)\n")
      cat("5 = CZEKANOWSKI (1913) or SORENSEN (1948) S7 coefficient of GOWER & LEGENDRE\n")
      cat("s5 = 2*a/(2*a+b+c) --> d = sqrt(1 - s)\n")
      cat("6 = S9 index of GOWER & LEGENDRE (1986)\n")
      cat("s6 = (a-(b+c)+d)/(a+b+c+d) --> d = sqrt(1 - s)\n")
      cat("7 = OCHIAI (1957) S12 coefficient of GOWER & LEGENDRE\n")
      cat("s7 = a/sqrt((a+b)(a+c)) --> d = sqrt(1 - s)\n")
      cat("8 = SOKAL & SNEATH (1963) S13 coefficient of GOWER & LEGENDRE\n")
      cat("s8 = ad/sqrt((a+b)(a+c)(d+b)(d+c)) --> d = sqrt(1 - s)\n")
      cat("9 = Phi of PEARSON = S14 coefficient of GOWER & LEGENDRE\n")
      cat("s9 = ad-bc)/sqrt((a+b)(a+c)(b+d)(d+c)) --> d = sqrt(1 - s)\n")
      cat("10 = S2 coefficient of GOWER & LEGENDRE\n")
      cat("s10 =  a/(a+b+c+d) --> d = sqrt(1 - s) and unit self-similarity\n")
      cat("Select an integer (1-10): ")
      method <- as.integer(readLines(n = 1))
    }

    a <- df %*% t(df)
    b <- df %*% (1 - t(df))
    c <- (1 - df) %*% t(df)
    d <- ncol(df) - a - b - c

    if (method == 1) {
      d <- a/(a + b + c)
    }
    else if (method == 2) {
      d <- (a + d)/(a + b + c + d)
    }
    else if (method == 3) {
      d <- a/(a + 2 * (b + c))
    }
    else if (method == 4) {
      d <- (a + d)/(a + 2 * (b + c) + d)
    }
    # correction d'un bug signalé par Christian Düring <c.duering@web.de>
    else if (method == 5) {
      d <- 2*a/(2 * a + b + c)
    }
    else if (method == 6) {
      d <- (a - (b + c) + d)/(a + b + c + d)

    }
    else if (method == 7) {
      d <- a/sqrt((a+b)*(a+c))
    }
    else if (method == 8) {
      d <- a * d/sqrt((a + b) * (a + c) * (d + b) * (d + c))
    }
    else if (method == 9) {
      d <- (a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) *
                                  (d + c))
    }
    else if (method == 10) {
      d <- a/(a + b + c + d)
      diag(d) <- 1
    }
    else stop("Non convenient method")
    d <- sqrt(1 - d)
    # if (sum(diag(d)^2)>0) stop("diagonale non nulle")
    d <- as.dist(d)
    attr(d, "Size") <- nlig
    attr(d, "Labels") <- d.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
  }

#function to compute mahalanobis distance
maha <-function(df=NULL){
  if(!is(df,"data.frame") && !is(df,"matrix")){
    stop("invalid class of Object")
  }
  df <- data.frame(df)
  nlig <- nrow(df)
  d <- matrix(0, nlig, nlig)
  d.names <- row.names(df)
  fun1 <- function(x) {
    sqrt(sum((df[x[1], ] - df[x[2], ])^2))
  }
  df <- as.matrix(df)
  index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
  dfcov <- cov(df) * (nlig - 1)/nlig
  maha <- eigen(dfcov, symmetric = TRUE)
  maha.r <- sum(maha$values > (maha$values[1] * 1e-07))
  maha.e <- 1/sqrt(maha$values[1:maha.r])
  maha.v <- maha$vectors[, 1:maha.r]
  maha.v <- t(t(maha.v) * maha.e)
  df <- df %*% maha.v
  d <- c(unlist(apply(index, 1, fun1)))
  upper=FALSE
  diag=FALSE
  attr(d, "Size") <- nlig
  attr(d, "Labels") <- d.names
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  class(d) <- "dist"
  return(d)
}


compbin<-function(X=NULL){
  if(!is(X,"data.frame") && !is(X,"matrix")){
    stop("invalid class of Object")
  }
  X<-na.omit(X)
  Unicos<-unlist(apply(X,1,unique))

  if(length(Unicos)==2){
    if (sort(unique(Unicos))==c(0,1)){
      return("Binary Data Imput")
    }else{
      return("Non-Binary Data Imput")
   }
  }else{
    return("Non-Binary Data Imput")
 }

}

#####Projection of each table in Compromise Configuration

TabProj<-function(S,SvdComp,Studies,Observations){
  As<-S%*%SvdComp$u%*%diag(1/sqrt(SvdComp$d))
  As<-cbind(as.data.frame(As),rep(Studies,nrow(As)),Observations)
  colnames(As)<-c(paste0("CP",seq(1:(ncol(As)-2))),"Studies","Observations")
  return(As)
}

####get_upper_tri is just to get the upper correlation from a matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

####Variable selection using LM
LinModel<-function(Mat=NULL,X=NULL,intercept=FALSE){


  X1<-matrix(as.numeric(X),ncol=1)
  #X1<-matrix(scale(X1),ncol=1)
  L<-as.matrix(Mat)

  if(intercept==TRUE){
    Model<-lm(X1~L,  na.action=na.exclude)
  }
  if(intercept==FALSE){
    Model<-lm(X1~-1+L,  na.action=na.exclude)
  }
  Res<-summary(Model)

  R2<-Res$adj.r.squared
  #p value from f statistic!
  pval<-1-pf(Res$fstatistic[1],Res$fstatistic[2],Res$fstatistic[3])
  if(intercept==TRUE){
    XCoord<-Res$coefficients[2,1]
    YCoord<-Res$coefficients[3,1]
  }
  if(intercept==FALSE){
    XCoord<-Res$coefficients[1,1]
    YCoord<-Res$coefficients[2,1]
  }
  if (is.na(R2) & is.na(pval)){
    M2<-NULL
  }else{
  M2<-data.frame(c(R2,pval,XCoord,YCoord))
  colnames(M2)<-rownames(X)
  rownames(M2)<-c("R2-Adj","p-value","XCoord","YCoord")
  }
  return(M2)
}

TSSfunction = function(x){
  x/sum(x)
}

