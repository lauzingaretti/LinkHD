#'@title OTU2Taxa
#'@description This function aggregates OTUs into their taxonomic characteristics
#'(genus or level)and it  analyses the most significant selected genera into each table.
#' To each genera, the function returns the hypergeometric distribution function P(x>=X)
#'  to each count. The function also returns
#' filtered data  by  counts higher than one. In both cases, we implemented -log(p+0.05),
#' then a higher value means more significant, i.e., it is an enrichment genus or family.
#'@param Selection list or data frame from VarSelection or dAB function
#'@param TaxonInfo data.frame with taxonomic table associated to Data input.
#'For instance, if Data comes from 16_S level, TaxoInfo should be a data.frame with 16_S
#' associated taxonomic information.
#'Note that the first column of this table must have the OTUs ids.
#'@param tableName a character indicating the table name. For instance, if your data comes
#' from 16_S, this parameter should be '16_S'.
#' Note that, this argument must mutch with the names from the input list into LinkData function.
#'@param AnalysisLev It is a character indicating if data should be aggregate to genera
#' or family level
#'@export OTU2Taxa
#'@return List. The first element of this list contains all the selected taxa with
#'their associated value from the hyperg distribution -log(p+0.05); the second element
#'of this list have only taxas counting up to 1.
#'
#'
#'
#'@author Laura M Zingatetti
#'
#'@references
#'{
#'\enumerate{
#'\item
#'\item Da Wei Huang, B. T. S., & Lempicki, R. A. (2009). Bioinformatics enrichment tools:
#' paths toward the comprehensive functional analysis of large gene lists.
#' Nucleic acids research, 37(1), 1.
#'\item Zheng, Q., & Wang, X. J. (2008). GOEAST: a web-based software toolkit for
#' Gene Ontology enrichment analysis. Nucleic acids research, 36(suppl_2), W358-W363.
#'}
#'}
#' @examples
#' {
#'data('Ruminotypes')
#'Normalization<-lapply(list(Ruminotypes$`16_S`,Ruminotypes$Archaea,Ruminotypes$`18_S`),
#'function(x){DataProcessing(x,Method='Compositional')})
#'Dataset<-Normalization
#'names(Dataset)<-c('16_S','Archaea','18_S')
#'#Running LinkData
#'Output<-LinkData(Dataset,Distance=rep('euclidean',3),
#'Scale = FALSE,Center=FALSE,nCluster = 3)
#'Select_Var<-VarSelection(Output,Data=Dataset,Crit = 'Rsquare',perc=0.9)
#'SignTaxa<-OTU2Taxa(Selection=VarTable(Select_Var),
#'TaxonInfo=Ruminotypes$Taxa_16S,tableName='16_S',AnalysisLev = 'Family')
#'Selected<-SignTaxa$TotalUp1
#' }
#'
#' @name OTU2Taxa
#' @rdname OTU2Taxa-OTU2Taxa
#' @import stats
#' @importFrom reshape2 melt





OTU2Taxa <- function(Selection, TaxonInfo, tableName, AnalysisLev = "Genus") {

    if (!is.data.frame(Selection) && !is.list(Selection)) {
        stop("Selection should be a list or data.frame")
    }

    if (length(tableName) > 1) {
        stop("Please, you must to indicate only one table name")
    }

    if (("Genus" %in% AnalysisLev || "Family" %in% AnalysisLev) == FALSE) {
        stop("You should choose an unit to data aggregation.")
    }

    if (is(Selection, "data.frame")) {
        if (tableName %in% Selection == FALSE) {
            stop("wrong table name specification. Please,
                 make sure that you are using a right table name")
        }
        SelectedVar <- Selection[Selection %in% tableName]
        # number of variables from the table
        ntable <- ncol(SelectedVar)
        SelectedVar <- names(SelectedVar)
    }
    if (is(Selection, "list")) {
        if (tableName %in% names(Selection) == FALSE) {
            stop("wrong table name specification.
                 Please, make sure that you are using a right table name")
        }
        # changed names changed names
        if (is(Selection[[which(names(Selection) %in% tableName)]], "character")) {
            SelectedVar <- Selection[[which(names(Selection) %in% tableName)]]
            ntable <- length(SelectedVar)
        }
        if (is(Selection[[which(names(Selection) %in% tableName)]], "numeric")) {
            SelectedVar <- names(Selection[[which(names(Selection) %in% tableName)]])
            ntable <- length(SelectedVar)
        }
    }

    TaxonInfo <- as.data.frame(TaxonInfo)
    if (AnalysisLev %in% colnames(TaxonInfo) == FALSE) {
        stop("Your table should contains genus or family level")
    }

    TotaltoGenLev <- table(TaxonInfo[colnames(TaxonInfo) == AnalysisLev])
    if (!any(suppressWarnings(is.na(as.numeric(SelectedVar))))) {
        TotalSelected <- table(TaxonInfo[TaxonInfo[, 1] %in% as.numeric(SelectedVar), ][, colnames(TaxonInfo) ==
            AnalysisLev])
    } else {
        if (any(TaxonInfo[, 1] %in% SelectedVar)) {
            TotalSelected <- table(TaxonInfo[TaxonInfo[, 1] %in% SelectedVar, ][, colnames(TaxonInfo) == AnalysisLev])
        } else {
            TotalSelected = 0
        }
    }
    if (is(TotalSelected, "table")) {
        # hypergeometric test
        Test <- c()
        for (i in seq_along(TotalSelected)) {
            x = TotalSelected[i]
            m = TotaltoGenLev[names(TotaltoGenLev) %in% names(TotalSelected)[i]]
            n = sum(TotaltoGenLev[-which(names(TotaltoGenLev) %in% names(TotalSelected)[i])])
            k = sum(TotalSelected)
            Test <- c(Test, 1 - phyper(x, m = m, n = n, k = k))
        }

        Test <- sort(Test, decreasing = TRUE)
        # look for some empty genera or family
        if (length(which(nchar(trimws(names(Test))) == 0)) > 0) {
            names(Test)[which(nchar(trimws(names(Test))) == 0)] <- paste0("unknown_", seq(seq_along(which(nchar(trimws(names(Test))) ==
                0))), AnalysisLev)
        }

        Test2 <- Test[names(Test) %in% names(TotalSelected[which(TotalSelected > 1)])]
        Test <- (-log(Test + 0.05))
        Test2 <- (-log(Test2 + 0.05))
        return(list(TotalHyp = Test, TotalUp1 = Test2))
    } else {
        return(list(TotalHyp = c("selected variables are not annotated,please check"), TotalUp1 = c("selected variables are not annotated, please check!")))
    }

}
