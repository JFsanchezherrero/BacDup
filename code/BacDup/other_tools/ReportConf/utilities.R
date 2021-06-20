########################################################################################################################
## utilities.R
## created: 2012-05-10
## creator: Yassen Assenov

##### Modified by Alba Moya Garces - 2021-04-30
##### Original script https://github.com/epigen/RnBeads/blob/master/R/utilities.R
## ---------------------------------------------------------------------------------------------------------------------
## Collection of (mostly internal) helper constants and functions.
########################################################################################################################

## F U N C T I O N S ###################################################################################################
########################################################################################################################

## parameter.is.flag
##
## Checks if the provided parameter value is a flag.
## 
## @param value Value to be tested.
## @return \code{TRUE} if \code{values} is \code{TRUE} or \code{FALSE}, \code{FALSE} otherwise.
## @author Yassen Assenov
parameter.is.flag <- function(value) {
  is.logical(value) && length(value) == 1 && (!is.na(value))
}

########################################################################################################################
## validate.single
## Validates the given vector or list contains a single element. This function is used in validating function or method
## arguments.
##
## @param x          Value vector or list to validate.
## @param param.name Name of parameter or slot that is validated. This is used in the generation of failing message.
## @return Short message that encodes the result of the validation in the form of a \code{character}. It is either the
##         string \code{ok}, or a short phrase describing the divergence from the "single value assumption".
## @author Yassen Assenov
validate.single <- function(x, param.name = "x") {
  if (is.null(x) || length(x) == 0) {
    result <- paste("missing value for", param.name)
  } else if (length(x) > 1) {
    result <- paste("multiple values for", param.name)
  } else if (is.na(x)) {
    result <- paste("missing value for", param.name)
  } else {
    result <- "ok"
  }
  return(result)
}

########################################################################################################################

## write.line
##
## Writes a line of text to the specified text file. This function is used in the generation of HTML reports.
##
## @param txt    Character vector storing the text to be written. The elements of this vector are concatenated without
##               a separator.
## @param fname  Name of the file to write the text to.
## @param indent Indentation of the text, given as number of \code{TAB} characters.
## @param append Flag indicating if the line is to be appended to the text file. If this is \code{FALSE}, the file's
##               contents are overwritten.
## @author Yassen Assenov

write.line <- function(txt, fname, indent = 0, append = TRUE) {
  strprefix <- paste(rep("\t", times = indent), collapse = "")
  cat(strprefix, paste0(txt, collapse = ""), "\n", file = fname, sep = "",
      append = append)
}

########################################################################################################################

## validate.dir
## If there is a logger initialized, validates that the given directory exists.
##
## @param dname Name of directory to be validated.
## @author Yassen Assenov
validate.dir <- function(dname) {
  if (logger.isinitialized()) {
    logger.validate.file(dname, is.file = FALSE)
  }
}

########################################################################################################################

