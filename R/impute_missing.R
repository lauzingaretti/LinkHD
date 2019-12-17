#'@title impute_missing
#'@description This function imputes missing abundances by sampling the known values from each variable.
#' The values are sampled according to the distribution, i.e. the most frequent values have a better chance of being chosen.
#'
#'@param data Data frame with abundances values. Note that this function should only be used with the raw data (counts)
#'@export impute_missing
#'@return Imputed data frame
#'@author Laura M Zingatetti
#'@importFrom plyr count

#' @examples
#' {
#'# toy example. To simulate 13 missing columns with less than 50 \% of missing values in each.
#'
#'data('Ruminotypes')
#'Data<-Ruminotypes$`16_S`
#'#13 indicates the number of columns with missing values.
#'Columns<-sample(1:ncol(Data),13)
#'for (i in Columns){
#'n<-sample(1:30,1)
#'Data[sample(1:nrow(Data),n),i]<-NA
#'}
#'
#'A<-impute_missing(Data)
#'#check precision of imputed data
#'cor(A[,Columns[1]],Ruminotypes$`16_S`[,Columns[1]])
#'cor(A[,Columns[2]],Ruminotypes$`16_S`[,Columns[2]])
#'}
#'
#'
#'
#' @name impute_missing
#' @rdname impute_missing


impute_missing<-function(data){

  if(!is.data.frame(data)){
    stop("input data should be a data.frame")
  }

  if (!all(floor(data) == data, na.rm = TRUE)){
    stop("data should be integers, i.e. abundances itself")
  }

  if(any(apply(data,2,function(x) sum(is.na(x)))==nrow(data))){
    cat("the column/s has/have missing data to all the samples", which(apply(data,2,function(x) sum(is.na(x)))==nrow(data)))
    stop("We cannot impute data if all column is missing. Please, consider to delete that column")
  }


  if(any(apply(data,2,function(x) sum(is.na(x)))>=0.8*nrow(data))){
    cat("the column/s has/have missing more than 80/% of samples missing", which(apply(data,2,function(x) sum(is.na(x)))==nrow(data)))
    warning("You should delete that columns and re-run the function ")
  }


  if(all(apply(data,2,function(x) sum(is.na(x)))==0)){
    stop("No missing values were found. You can continue with the analysis")
  }

  B<- apply(data,2,function(x){ sample(count(na.omit(x))[,1], size=sum(is.na(x))
                                       , prob=count(na.omit(x))[,2]/sum(count(na.omit(x))), replace=TRUE)
  })



  if(any(unlist(lapply(B,function(x) length(x)))>0)){
    index<-which(unlist(lapply(B,function(x) length(x)))>0)
    dfz<-data[,index]
    BB<-B[index]

    Nom<-names(BB)

    for (i in (1:ncol(dfz))){
      dfz[is.na(dfz[,i]),][,i]<-BB[i]
    }
    data[,index]<-dfz
  }

  return("data"=data)
}
