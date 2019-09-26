#' @title DataProcessing
#' @description  function to Perform external datas' pre-processing.
#' This function allows an external pre-processing of the datasets including on the analysis in three ways: Standard, Compositional (centered log ratio) and
#' frequencies.
#'
#' @param Data a numeric data.frame.
#' @param Method character indicating the method used to Data preprocessing. If data are continous, use "Standard".
#' If Data are compositional, please use "Compositional" and clr (centered log-ratios functions) transformations are performed.
#' To compositional data, you also could use the option "TSS" Total Sum Scaling follow up bray (Bray-Curtis) in distance option.
#' The function also allows to processing frequencies- like data through "FreqNorm" option.
#' Note that when you use Compositional, we first sum 1 to all the counts (in order to performs the log transformation before).
#' @return a data.frame with normalized data.
#' @export DataProcessing
#' @author Laura M Zingatetti
#'
#'
#' @examples
#' {
#' data(Taraoceans)
#' Data<-Taraoceans$phychem
#' Data<-DataProcessing(Data,Method="Standard")
#' }
#

DataProcessing<-function(Data=NULL,Method="Standard"){

            if (is.null(Data)){
              stop("You should include a data.frame to standarized")
            }

            if(Method!="Standard" & Method!="Compositional" & Method!="FreqNorm" & Method!="TSS"){
              stop("You have to choose a valid Method to data standarization")
            }

            if(Method=="Standard"){
              Data<-scale(Data,scale=TRUE,center=TRUE)
            }


            if(Method=="FreqNorm"){
              Data<-cia(Data)
            }

            if(Method=="Compositional"){
            Data<-centerLR(Data+1)$Data.clr
            }

            if(Method=="TSS"){
            Data<- t(apply(Data+1, 1, TSSfunction))
               }

            return(Data)
            }
