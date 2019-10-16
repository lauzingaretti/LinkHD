#' @title Read_Data: a fast way to data reading.
#' @description  this function read all dataset in a folder and returns
#' list needed to Link_Data function input.
#'@param Path path to folder containing all dataset to integrate
#'
#'@return \item{List}{List including all dataset into the parent directory.
#' List names inherit the  names of the files}
#'@export Read_Data
#'@author Laura M Zingatetti
#'
#'
#' @examples
#'\dontrun{
#'Datos<-Read_Data('Path to parent folder',common_elements=1)
#'}
#'
#'
#'
#' @name Read_Data
#' @rdname Read_Data
#' @importFrom data.table fread
#' @importFrom rio convert
#'
Read_Data <- function(Path = "") {
    if (!dir.exists(Path)) {
        stop("Path is not valid")
    } else {
        # Create a vector of Excel files to read
        if (length(list.files(Path, pattern = "xls")) > 0 | length(list.files(Path, pattern = "xlsx")) > 0) {
            f1 = list.files(Path, pattern = "xlsx")
            f2 = list.files(Path, pattern = "xls")
            if (length(f1) > 0) {
                created <- mapply(convert, f1, gsub("xlsx", "csv", f1))
                unlink(f1)
            }
            if (length(f2) > 0) {
                created <- mapply(convert, f2, gsub("xls", "csv", f2))
                unlink(f2)
            }
        }
        Dat <- setdiff(list.files(Path, full.names = TRUE), list.dirs(Path, recursive = FALSE))

    }

    Datos <- list()
    for (i in seq_along(Dat)) {
        Datos[[i]] <- data.frame(fread(Dat[[i]], header = TRUE), row.names = 1, check.names = FALSE)
    }
    names(Datos) <- Dat
    return(Datos)

}




