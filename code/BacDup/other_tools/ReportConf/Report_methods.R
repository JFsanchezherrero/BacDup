########################################################################################################################
## Report-methods.R
## created: 2012-05-08
## creator: Yassen Assenov

##### Modified by Alba Moya Garces - 2021-04-15
##### Original script https://github.com/epigen/RnBeads/blob/369833af2810371f433909ad202dedaf1fa10e5d/R/Report-methods.R
## ---------------------------------------------------------------------------------------------------------------------
## Report class method definition.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

# Checks if the given strings are valid to be used as file names.
#
#  fname character} vector of candidate strings for names.
# @return logical} vector showing, for each element of fname}, if it meets the criteria for a vaild file
#         name.
# @author Fabian Mueller
# @noRd
is.valid.fname <- function(fname) {
  grepl("^[A-Za-z0-9._-]+$", fname)
}

########################################################################################################################

# Checks if the given path is valid to be used as a relative link (file or directory) in an HTML report.
#
#  dname Name of the file or directory.
# @return TRUE} if dname} meets the criteria for a file or directory to be easily referenced in an HTML
#         link, FALSE} otherwise.
# @author Yassen Assenov
# @noRd
is.valid.relative <- function(dname) {
  if (!grepl("^[A-Za-z0-9/\\._-]+$", dname)) {
    return(FALSE)
  }
  if (grepl("^/", dname)) {
    return(FALSE)
  }
  if (grepl("/$", dname)) {
    return(FALSE)
  }
  return(TRUE)
}

########################################################################################################################

# Creates the given directory. This function attempts to create all non-existant elements of the given path.
#
#  dname           Directory to be created.
#  accept.existing Flag indicating if an existing empty directory should be accepted. When this is TRUE}
#                        and dname} points to an existing empty directory, this function returns TRUE}.
#                        Hidden files (file names starting with "."} on Unix platforms) are ignored. Thus, if
#                        dname} already exists and contains hidden files only, it is treated as empty.
#  showWarnings    Flag indicating if warnings should be shown.
# @return                logical} indicating if the creation was successful.
#
# @seealso \link{dir.create}}
# @author Yassen Assenov
# @noRd
create.path <- function(dname, accept.existing = TRUE, showWarnings = TRUE) {
  if (file.exists(dname)) {
    if (accept.existing) {
      if (file.info(dname)[1, "isdir"]) {
        ## TODO: Should invisible files also be considered?
        return(length(dir(dname)) == 0)
      }
    }
    return(FALSE)
  }
  return(dir.create(dname, showWarnings = showWarnings, recursive = TRUE))
}

########################################################################################################################

# Completes the given HTML report by adding a footer notice and closing the <body>} and <html>} tags.
#
#  report Report instance to close.
#  rtype  Contents of the HTML file. It should be "report"} or "index"}.
# @author Yassen Assenov
# @noRd
complete.report <- function(report, report.type = "report") {
  
  ## Add list of references to the report
  if (length(report@references) != 0) {
    report <-  bdp.add.section(report, "References", NULL)
    reftexts <- report@references
    reftexts <- paste0("<a id=\"ref", 1:length(reftexts), "\">", reftexts, "</a>")
    bdp.add.list(report, as.list(reftexts), type = "o")
  }
  
  ## Close any open sub-sections and sections
  for (i in length(report@opensections):1) {
    if (report@opensections[i] != 0) {
      write.line('</div>', report@fname)
      report@opensections[i] <- 0L
    }
  }
  
  ## Add footer to the report
  
  
  ##TODO terminar de decidir qué poner
  
  write.line("\n<div id=\"copyright\">", report@fname)
  write.line("\t<div id=\"bacdup\">", report@fname)
  write.line(c("\tThis ", report.type, " was generated on ", format(Sys.time(), "%Y-%m-%d"), " by ",
               "<a href=\"http:\">BacDup</a> version "
               #,
               #paste(as.character(utils::packageVersion("RnBeads")), collapse = "."), "."
               ), report@fname)
  write.line("\t</div>", report@fname)
  # write.line("\t<div id=\"validlogo\">", report@fname)
  # write.line(c("\t\t<a href=\"http://validator.w3.org/check?uri=referer\">",
  #              "<img alt=\"Valid XHTML 1.1\" height=\"31\" src=\"http://www.w3.org/Icons/valid-xhtml11-blue\"",
  #              " width=\"88\" /></a>"), report@fname)
  # write.line(c("\t\t<a href=\"http://jigsaw.w3.org/css-validator/check/referer\">",
  #              "<img alt=\"Valid CSS\" height=\"31\" src=\"http://jigsaw.w3.org/css-validator/images/vcss-blue\"",
  #              " width=\"88\" /></a>"), report@fname)
  # write.line("\t</div>", report@fname)
  write.line("</div>", report@fname)
  
  ## Close all tags
  write.line("</body>", report@fname)
  write.line("</html>", report@fname)
  return(report)
}

########################################################################################################################

setMethod("show", "Report",
          function(object) {
            cat("Object of class Report - an XHTML Report in BacDup\n")
            infotext <- c("Report's file" = object@fname,
                          "Directory for configuration files" = object@dir.conf,
                          "Directory for data files" = object@dir.data,
                          "Directory for figure image files" = object@dir.fig,
                          "Directory for figure PDF files" = object@dir.pdfs,
                          "Directory for high-resolution image files" = object@dir.high,
                          "Number of references" = length(object@references))
            n <- max(nchar(names(infotext)))
            for (i in 1:length(infotext)) {
              cat(sprintf(paste0("%", n, "s: %s\n"), names(infotext)[i], infotext[i]))
            }
          }
)

########################################################################################################################

# rnb.initialize.reports -> cambio a bdp
#
# Creates a new directory to host HTML reports and copies the shared configuration files.
#
#  dir.reports       Directory to host report files. This must be a character} of length one that specifies
#                          a non-existent path, as this methods attempts to create it.
#  dir.configuration Subdirectory to host configuration files shared by the reports. This must be a
#                          character} of length one that gives location as a path relative to dir.reports}.
#                          Also, strong restrictions apply to the path name. See the description of the
#                          \link{createReport}} function for more details. This method creates the directory and
#                          copies configuration files that define cascading style sheet (CSS) definitions and
#                          Javascript functions used by the HTML reports.
# return TRUE if the report directory was successfully created and the configuration files were copied to the
#         specified location; FALSE otherwise.
#
#
# @author Yassen Assenov
# @export
bdp.initialize.reports <- function(dir.reports, dir.configuration = "configuration") {
  if (!(is.character(dir.reports) && length(dir.reports) == 1 && (!is.na(dir.reports[1])))) {
    stop("invalid value for dir.reports")
  }
  if (!(is.character(dir.configuration) && length(dir.configuration) == 1 && (!is.na(dir.configuration[1])))) {
    stop("invalid value for dir.configuration")
  }
  if (!is.valid.relative(dir.configuration)) {
    stop("invalid value for dir.configuration")
  }
  ## Remove trailing slash or backslash in dir.reports (if it exists)
  dir.reports <- sub("^(.+)[/\\\\]$", "\\1", dir.reports[1])
  ## Attempt to create dir.reports
  if (!create.path(dir.reports, FALSE, showWarnings = FALSE)) {
    return(FALSE)
  }
  ## Attempt to create dir.configuration
  dname <- normalizePath(file.path(dir.reports, dir.configuration), mustWork = FALSE)
  if (!create.path(dname, showWarnings = FALSE)) {
    return(FALSE)
  }
  ## Copy configuration files
    ##TODO crear una imagen BacDup
  cfiles <- c("arrow_down.png", "arrow_right.png", "pdf_active.png", "pdf_inactive.png", "report.css",
              "report.js")
  cfiles <- file.path("ReportConf/configuration", cfiles)
  return(all(file.copy(cfiles, dname)))
}

## M E T H O D S #######################################################################################################

setValidity("Report",
            function(object) {
              res <- validate.single(object@fname, "fname")
              if (res != "ok") { return(res) }
              res <- validate.single(object@dir.conf, "dir.conf")
              if (res != "ok") { return(res) }
              res <- validate.single(object@dir.data, "dir.data")
              if (res != "ok") { return(res) }
              res <- validate.single(object@dir.fig, "dir.fig")
              if (res != "ok") { return(res) }
              res <- validate.single(object@dir.pdfs, "dir.pdfs")
              if (res != "ok") { return(res) }
              res <- validate.single(object@dir.high, "dir.high")
              if (res != "ok") { return(res) }
              res <- validate.single(object@figures, "figures")
              if (res != "ok") { return(res) }
              res <- validate.single(object@tables, "tables")
              if (res != "ok") { return(res) }
              regex.fname <- "^[A-Za-z0-9/\\._-]+$"
              if (!grepl(regex.fname, basename(object@fname))) {
                return("fname is invalid")
              }
              if (!grepl(".+\\.(htm|html|xhtml|xml)$", tolower(basename(object@fname)))) {
                return("fname is invalid")
              }
              if (!is.valid.relative(object@dir.data)) {
                return("dir.data is invalid")
              }
              if (!is.valid.relative(object@dir.fig)) {
                return("dir.fig is invalid")
              }
              if (!is.valid.relative(object@dir.pdfs)) {
                return("dir.pdfs is invalid")
              }
              if (!is.valid.relative(object@dir.high)) {
                return("dir.high is invalid")
              }
              if (object@sections < 0) {
                return("sections must be a non-negative value")
              }
              if (object@figures < 0) {
                return("figures must be a non-negative value")
              }
              if (object@tables < 0) {
                return("tables must be a non-negative value")
              }
              TRUE
            }
)

########################################################################################################################

setMethod("initialize", "Report",
          function(.Object, fname, title, page.title, authors, dirs, init.configuration) {
            .Object@fname <- normalizePath(fname, mustWork = FALSE)
            ## Read or generate report-specific directory names
            dir.configuration <- dirs["configuration"]
            if (is.na(dir.configuration)) {
              dir.configuration <- "configuration"
            }
            if (!is.valid.relative(dir.configuration)) {
              stop("dir.configuration is invalid")
            }
            .Object@dir.conf <- dir.configuration
            dprefix <- sub("^(.*[\\\\/])*([^\\\\/]+)\\.[^\\.]+$", "\\2", fname)
            .Object@dir.data <- dirs["data"]
            if (is.na(.Object@dir.data)) {
              .Object@dir.data <- paste(dprefix, "data", sep = "_")
            }
            .Object@dir.fig <- dirs["fig"]
            if (is.na(.Object@dir.fig)) {
              .Object@dir.fig <- paste(dprefix, "images", sep = "_")
            }
            .Object@dir.pdfs <- dirs["pdfs"]
            if (is.na(.Object@dir.pdfs)) {
              .Object@dir.pdfs <- paste(dprefix, "pdfs", sep = "_")
            }
            .Object@dir.high <- dirs["high"]
            if (is.na(.Object@dir.high)) {
              .Object@dir.high <- .Object@dir.fig
            }
            validObject(.Object)
            if (!is.null(authors) && length(authors) != 0) {
              if (!all(grepl("^[A-Za-z ,.-]+$", authors))) {
                stop("invalid value for authors")
              }
              authors <- paste(authors, collapse = ", ")
            } else {
              authors <- NULL
            }
            
            ## Create directories
            dnames <- file.path(dirname(.Object@fname), c(.Object@dir.data, .Object@dir.fig, .Object@dir.pdfs, .Object@dir.high))
            dnames <- setdiff(unique(normalizePath(dnames, mustWork = FALSE)), dirname(.Object@fname))
            for (dname in dnames) {
              if (!create.path(dname)) {
                stop(paste("directory", dname, "cannot be created"))
              }
            }
            rm(dnames, dname)
            
            ## Validate configuration directory; create if it's missing
            if (init.configuration) {
              dname <- normalizePath(file.path(dirname(.Object@fname), dir.configuration), mustWork = FALSE)
              if (exists(dname)) {
                if (file.info(dname)[1, "isdir"] == FALSE) {
                  stop(paste(dname, "is not a directory"))
                }
              } else if (!create.path(dname)) {
                stop(paste("directory", dname, "cannot be created"))
              }
              ## Copy configuration files
              cfiles <- c("arrow_down.png", "arrow_right.png", "BacDup.png", "pdf_active.png", "pdf_inactive.png",
                         "report.css", "report.js")
              #cfiles <- "~/git/BacDup/developer/test_R/configuration"
              cfiles <- file.path("ReportConf/configuration", cfiles)
              if (!all(file.copy(cfiles, dname))) {
                stop(paste("configuration could not be initialized in", dname))
              }
              rm(dname, cfiles)
            }
            
            ## Create the HTML file and print its header
            dir.configuration <- ifelse(dir.configuration != ".", paste0(dir.configuration, "/"), "")
            wline <- function(txt, indent = 0) {
              write.line(txt, fname, indent = indent)
            }
            write.line("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", fname, append = FALSE)
            wline("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" \"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">")
            wline("<html xmlns=\"http://www.w3.org/1999/xhtml\">")
            wline("<head>")
            wline("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />", 1)
            if (!is.null(authors)) {
              wline(c("<meta name=\"author\" content=\"", authors, "\" />"), 1)
            }
            wline(c("<link rel=\"stylesheet\" href=\"", dir.configuration, "report.css\" type=\"text/css\" />"), 1)
            if (!is.null(page.title)) {
              wline(c("<title>", page.title[1], "</title>"), 1)
            }
            wline(c("<script src=\"", dir.configuration, "report.js\" type=\"text/javascript\"></script>"), 1)
            wline("</head>\n")
            wline(c("<body style=\"background:url(", dir.configuration, "BacDup.png) no-repeat 5px 25px;\">\n"))
            wline(c("<h1>", title, "</h1>\n"))
            .Object
          }
)

########################################################################################################################

# bdp.get.directory
#
# Gets the location of the given report-specific directory.
#
#  report   Report of interest.
#  dir      Type of directory to get. Must be one of "data"}, "images"}, "images-high"} or
#                 "pdfs"}.
#  absolute Flag indicating if the absolute path of the directory is to be returned. If this is FALSE},
#                 the directory name is returned relative to the report's HTML file location.
# 
# @seealso \linkS4class{Report}} for functions adding contents to an HTML report
# @author Yassen Assenov

bdp.get.directory <- function(report, dir = c("data", "images", "images-high", "pdfs"), absolute = FALSE) {
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (is.null(dir) || length(dir) == 0 || (!(dir[1] %in% c("data", "images", "images-high", "pdfs")))) {
    stop("invalid value for dir")
  }
  if (!parameter.is.flag(absolute)) {
    stop("invalid value for absolute; expected TRUE or FALSE")
  }
  dir <- dir[1]
  if (dir == "data") {
    result <- report@dir.data
  } else if (dir == "pdfs") {
    result <- report@dir.pdfs
  } else if (dir == "images") {
    result <- report@dir.fig
  } else { # dir == "images-high"
    result <- report@dir.high
  }
  if (absolute) {
    result <- normalizePath(file.path(dirname(report@fname), result), mustWork = FALSE)
  }
  return(result)
}

########################################################################################################################

# bdp.add.section
#
# Generates HTML code for a new section in the specified report.
#
#  report      Report to write the text to.
#  title       Section header. This must be a single-element character vector.
#  description Human-readable paragraph text of the section in the form of a character vector. Elements
#                    of this vector are concatenated without a separator to form the full description. Set this to
#                    NULL if the section does not (yet) contain text.
#  level       Section level as a single integer. It must be one of 1, 2 or 3,
#                    denoting section, subsection and sub-subsection, respectively.
#  collapsed   Flag indicating if the contents of this section is to be initially collapsed. Possible values are
#                    TRUE (the section is not visible), FALSE (default, the section is expanded) and
#                    "never" (the section cannot be collapsed or expanded).
#
# @author Yassen Assenov

bdp.add.section <- function(report, title, description, level = 1L, collapsed = FALSE) {
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!is.character(title) || validate.single(title) != "ok") {
    stop("invalid value for title")
  }
  if (!(is.null(description) || is.character(description))) {
    stop("invalid value for description")
  }
  if (is.numeric(level) && all(level == as.integer(level))) {
    level <- as.integer(level)
  }
  if (!(is.integer(level) && length(level) == 1 && 1 <= level[1] && level[1] <= 3)) {
    stop("invalid value for level")
  }
  if (is.character(collapsed)) {
    if (!(length(collapsed) == 1 && isTRUE(collapsed == "never"))) {
      stop('invalid value for collapsed; expected TRUE, FALSE or "never"')
    }
    hclass <- 'fixed'
  } else {
    if (!parameter.is.flag(collapsed)) {
      stop('invalid value for collapsed; expected TRUE, FALSE or "never"')
    }
    hclass <- ifelse(collapsed, 'collapsed', 'expanded')
  }
  ## Close previous sub-sections and sections if necessary.
  for (i in length(report@opensections):level) {
    if (report@opensections[i] != 0) {
      write.line('</div>', report@fname)
      report@opensections[i] <- 0L
    }
  }
  
  i <- report@sections <- report@sections + 1L
  attr.text <- ifelse(is.logical(collapsed), paste0(' onclick="toggleSection(\'', i, '\')"'), '')
  attr.text <- paste0(' class="', hclass, '" id="hsection', i, '"', attr.text)
  write.line(c("\n<h", (level + 1), attr.text, ">", title, "</h", (level + 1), ">"), report@fname)
  
  if (is.logical(collapsed)) {
    txt <- c('<div id="section', i, '" style="display:', ifelse(collapsed, 'none', 'block'), '">')
    write.line(txt, report@fname)
    
    ## Mark the newly open section
    report@opensections[level] <- i
  }
  if (!is.null(description)) {
    write.line(c("<p>", description, "</p>"), report@fname)
  }
  return(report)
}

########################################################################################################################

# bdp.add.paragraph
#
# Generates HTML code for a new paragraph in the specified report.
#
#  report          Report to write the text to.
#  txt             character} vector (or array) storing the text to be written. The elements of this vector
#                        are concatenated without a separator.
#  paragraph.class CSS class definition of the paragraph. This must be either NULL (default), 
#                         "centered" (This paragraph gives a formula or a short statement. Text is
#                               horizontally centered.
#                          "note" (This paragraph describes a note. Text is italic.)
#                          "task" (This paragraph describes a task. Text is bold and bright red.

# @author Yassen Assenov

bdp.add.paragraph <- function(report, txt, paragraph.class = NULL) {
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!is.character(txt)) {
    stop("invalid value for txt; expected character vector")
  }
  if (is.null(paragraph.class)) {
    paragraph.class <- ""
  } else {
    if (!(is.character(paragraph.class) && length(paragraph.class) == 1 && (!is.na(paragraph.class)))) {
      stop("invalid value for paragraph.class; expected single character")
    }
    paragraph.class <- paragraph.class[1]
    if (!(paragraph.class %in% c("centered", "note", "task"))) {
      stop("invalid value for paragraph.class; expected one of centered, note, task")
    }
    paragraph.class <- paste(" class=", paragraph.class, "", sep = "\"")
  }
  write.line(c("<p", paragraph.class, ">", txt, "</p>"), report@fname)
  return(invisible(report))
}

########################################################################################################################

# bdp.add.list
#
# Generates HTML code for a list in the specified report.
#
#  report Report to write the text to.
#  txt    Non-empty list of items to be written. An attribute named type}, if it exists, specifies the type
#               of the list. See the Details section for more information. Every item must be either a nested
#               list, denoting a sublist, or a character vector (or array), storing the text to be
#               written. Any other objects are coerced to a character type. Elements are concatenated without a
#               separator to form the text for a list item.
#  type   List type to be used for the list and/or its sublists in case the attribute type is not
#               specified.

# @details There are two ways to specify a list type: (1) setting a value for the attribute type of the list, or
#          (2) using the function's parameter type. The value of the function's parameter is used only for lists
#          and sublists that do not contain an attribute named type. The following types are supported:
#          "o" Ordered list using arabic numbers - 1, 2, 3, etc.
#          "u" Unordered list using bullet points.
#          
#          Note that every list type must be a one-element character vector containing one of the codes listed
#          above. Specifying any other value for list type results in an error.
#

# @author Yassen Assenov

bdp.add.list <- function(report, txt, type = "u") {
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!(is.list(txt) && length(txt) != 0)) {
    stop("invalid value for txt; expected non-empty list")
  }
  if (!(is.character(type) && length(type) == 1)) {
    stop("invalid value for type; character expected")
  }
  type <- tolower(type)
  ## TODO: Add support for more list types
  if (!(type %in% c("o", "u"))) {
    stop("invalid value for type; expected one of: o, u")
  }
  get.html <- function(txt, indent) {
    ltype <- attr(txt, "type")
    if (is.null(ltype)) { ltype <- type }
    if (ltype == "u") {
      tagtext <- "ul"
    } else if (ltype == "o") {
      tagtext = "ol"
    } else {
      stop("invalid list type; expected one of: o, u")
    }
    indenttext <- paste(rep("\t", times = indent), collapse = "")
    contents <- sapply(txt, function(litem) {
      if (is.list(litem)) {
        paste0("\t<li>\n", paste(get.html(litem, indent + 2), collapse = ""), "\n\t</li>\n")
      } else {
        paste0("\t<li>", paste(litem, collapse = ""), "</li>\n")
      }
    })
    paste0(indenttext, c(paste0("<", tagtext, ">\n"), contents, paste0("</", tagtext, ">")))
  }
  write.line(get.html(txt, 0L), report@fname)
  return(invisible(report))
}

########################################################################################################################

# bdp.add.table
#
# Generates HTML code for a table in the specified report.
#
#  report           Report to write the text to.
#  tdata            Matrix or data frame to be presented in HTML form. Column names, if present, are used to
#                         define table columns. If this table contains 0 (zero) rows or 0 columns, calling this
#                         function has no effect.
#  row.names        Flag indicating if row names should also be printed. If this parameter is TRUE and
#                         tdata defines row names, these are printed in the left-most column and are displayed
#                         as header cells. Keep in mind that data.frame always define row names.
#  first.col.header Flag indicating if all cells in the first column must be displayed as header cells. Note
#                         that, if both this parameter and row.names are TRUE, and tdata contains
#                         row names, the constructed HTML table will have 2 columns of header cells.
#  indent           Default indentation, in number of tabulation characters, to apply to HTML tags. This
#                         indentation is also applied to thead.
#  tag.attrs        Named character vector specifying the list of attributes to be set to the
#                         <table> element. Setting this to NULL or an empty character vector
#                         disables attributes.
#  thead            character vector storing a table header to include. This can, for example, be a
#                         character that defines column widths. Every element in this vector is written on a
#                         separate line, applying the indentation given by indent.
#  tcaption         Text to include as a caption below the table, or NULL if the table does not contain
#                         caption.
#  na               character to be used for printing NA values in the table. This parameter is
#                         not considered when printing thead or the table's column names.

# @author Yassen Assenov

bdp.add.table <- function(report, tdata, row.names = TRUE, first.col.header = FALSE, indent = 0,
                          tag.attrs = c("class" = "tabdata"), thead = NULL, tcaption = NULL, na = "<span class=\"disabled\">n/a</span>") {
  
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!(is.matrix(tdata) || is.data.frame(tdata))) {
    stop("invalid value for tdata; expected matrix or data.frame")
  }
  if (!parameter.is.flag(first.col.header)) {
    stop("invalid value for first.col.header; expected TRUE or FALSE")
  }
  if (is.numeric(indent) && all(indent == as.integer(indent))) {
    indent <- as.integer(indent)
  }
  if (!(is.integer(indent) && length(indent) == 1 && (!is.na(indent)) && 0 <= indent)) {
    stop("invalid value for indent; expected single non-negative integer")
  }
  if (!is.null(tag.attrs)) {
    if (!(is.character(tag.attrs) && length(tag.attrs) == length(names(tag.attrs)))) {
      stop("invalid value for tag.attrs; expected named character vector")
    }
    if (length(tag.attrs) != 0) {
      tag.attrs <- paste("", paste(names(tag.attrs), "=\"", tag.attrs, "\"", sep = "", collapse = " "))
    } else {
      tag.attrs <- ""
    }
  } else {
    tag.attrs <- ""
  }
  if (!is.null(thead)) {
    if (!(is.character(thead))) {
      stop("invalid value for thead; expected character")
    }
    if (length(thead) == 0) {
      thead <- NULL
    }
  }
  if (!is.null(tcaption)) {
    if (!is.character(tcaption)) {
      stop("invalid value for tcaption; expected character")
    }
  }
  if (!(is.character(na) && length(na) == 1)) {
    stop("invalid value for na; expected single character")
  }
  if (nrow(tdata) * ncol(tdata) == 0) {
    return(invisible(report))
  }
  
  ## Convert factors to characters
  if (is.data.frame(tdata)) {
    for (j in 1:ncol(tdata)) {
      if (is.factor(tdata[[j]])) {
        tdata[[j]] <- as.character(tdata[[j]])
      }
    }
  }
  
  wline <- function(txt, indentation = indent) {
    write.line(txt, report@fname, indent = indentation)
  }
  wline(c("<div class=\"restricted\"><table", tag.attrs, ">"))
  if (!is.null(thead)) {
    mapply(wline, thead)
  }
  startcolumn <- as.integer(ifelse(row.names & (!is.null(rownames(tdata))), 0, 1))
  if (!is.null(colnames(tdata))) {
    wline("<tr>")
    if (startcolumn == 0) {
      cnames <- c("", colnames(tdata))
    } else {
      cnames <- colnames(tdata)
    }
    for (cname in cnames) {
      wline(c("<td class=\"header\">", cname, "</td>"), indent + 1)
    }
    wline("</tr>")
    rm(cnames)
  }
  if (nrow(tdata) != 0) {
    for (i in 1:nrow(tdata)) {
      wline("<tr>")
      for (j in startcolumn:ncol(tdata)) {
        value <- base::ifelse(j == 0, rownames(tdata)[i], tdata[[i, j]])
        if (is.na(value)) { value <- na }
        cclass <- ifelse(j == 0 || (j == 1 && first.col.header), "header", "centered")
        wline(c("<td class=\"", cclass, "\">", value, "</td>"), indent + 1)
      }
      wline("</tr>")
    }
  }
  wline("</table></div>")
  if (!is.null(tcaption)) {
    wline(c("<p class=\"centered\"><span class=\"note\">", tcaption, "</span></p>"))
  }
  wline("", indentation = 0)
  return(invisible(report))
}

########################################################################################################################

# bdpreport.get.element.values
#
#  elements       Observed values of element names, in the form of a list} of character} vectors of
#                       equal length. The <i>i</i>-th element of this list should be extracted from the <i>i</i>-th
#                       table identifier or plot file name.
#  settings.count Number of variable elements for which descriptions are provided. The is the length of
#                       setting.names}.
# @return List of observed values per element; an empty list if no variation is observed (elements if of length
#         one). In case the number of variable elements does not match settings.count, the return value is
#         character containing an error message.
# @author Yassen Assenov

bdpreport.get.element.values <- function(elements, settings.count) {
  element.count <- range(sapply(elements, length))
  if (element.count[1] != element.count[2]) {
    return(NULL)
  }
  element.count <- element.count[1]
  if (element.count < settings.count) {
    return("too many entries in setting.names")
  }
  element.values <- lapply(1:element.count, function(i) { sort(unique(sapply(elements, "[", i))) })
  names(element.values) <- 1:length(element.values)
  i <- which(sapply(element.values, length) != 1)[1]
  if ((!is.na(i)) && settings.count < length(element.values) - i + 1) {
    return("too few entries in setting.names")
  }
  return(element.values)
}

########################################################################################################################

# bdpreport.create.settings.table
#
# Creates the HTML code for a table of selecting settings. This is used in figures and in listings of tables.
#
#  type              One of "fig" or "tab".
#  n                 Index of figure or table listing in the report.
#  element.values    List of all observed values for each element in the identifiers. This is extracted by
#                     bdpreport.get.element.values.
#  setting.names     List of element descriptors, as provided to the calling function.
#  selected.elements Element values that should be initially selected.
#  indent            Indentation, in number of tabulation characters, of the <table> tag.

# @author Yassen Assenov

bdpreport.create.settings.table <- function(type, n, element.values, setting.names, selected.elements, indent = 1L) {
  updatefun <- ifelse(type == "fig", "updateFigure", "updateTable")
  result <- character()
  wline <- function(txt, ind = 0L) {
    result <<- c(result, paste(c(paste(rep("\t", indent + ind), collapse = ""), txt), collapse = ""))
  }
  i <- length(element.values) - length(setting.names) + 1
  element.values <- element.values[i:length(element.values)]
  ## Validate setting.names
  missing.names <- sapply(1:length(element.values), function(i) {
    length(setdiff(element.values[[i]], names(setting.names[[i]]))) != 0
  })
  if (any(missing.names)) {
    return(NULL)
  }
  
  ## Create a table of settings
  wline("<table class=\"tabsettings\">")
  for (i in 1:length(element.values)) {
    si <- as.integer(names(element.values)[i])
    wline("<tr>")
    wline(c("<td>", names(setting.names)[i], "</td>"), 1)
    wline("<td>", 1)
    txt <- ifelse(length(element.values[[i]]) == 1, " disabled=\"disabled\"", "")
    txt <- c("<select", txt, " id=\"", type, n, "setting", si)
    wline(c(txt, "\" onchange=\"", updatefun, "(\'", type, n, "\')\">"), 2)
    s.options <- setting.names[[i]]
    s.options <- s.options[intersect(names(s.options), element.values[[i]])]
    for (j in 1:length(s.options)) {
      sv <- names(s.options)[j]
      txt <- ifelse(selected.elements[si] == sv, " selected=\"selected\"", "")
      wline(c("<option", txt, " value=\"", sv, "\">", s.options[j], "</option>"), 3)
    }
    wline("</select>", 2)
    wline("</td>", 1)
    wline("</tr>")
  }
  wline("</table>")
  wline("<br />")
  
  return(paste(result, collapse = "\n"))
}

########################################################################################################################

# bdp.add.tables
#
# Generates HTML code for a listing of tables (of which only one is visible at any moment) in the specified report.
#
#  report         Report to write the text to.
#  tables         Non-empty list of tables, each one represented by a data.frame or
#                       matrix. The names of this list are used as table identifiers; each one consists
#                       of elements separated by underscore character (_).
#  setting.names  List of table name element descriptors. Every variable elements in the table names must be
#                       included in this list.
#  selected.table Index of the table to be initially selected in this listing.
#  indent         Default indentation, in number of tabulation characters, to apply to every table.
#  ...            Other parameters passed to bdp.add.table.

# @author Yassen Assenov

bdp.add.tables <- function(report, tables, setting.names, selected.table = 1L, indent = 2L, ...) {
  
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!(is.list(tables) && length(tables) != 0 &&
        all(sapply(tables, function(x) { is.matrix(x) || is.data.frame(x) })))) {
    stop("invalid value for tables")
  }
  if (!(is.list(setting.names) && all(sapply(setting.names, is.vector)))) {
    stop("invalid value for setting.names")
  }
  if (is.numeric(selected.table) && all(as.integer(selected.table) == selected.table)) {
    selected.table <- as.integer(selected.table)
  }
  if (!is.integer(selected.table) && length(selected.table) == 1 && (!is.na(selected.table))) {
    stop("invalid value for selected.table")
  }
  if (selected.table < 1 || length(tables) < selected.table) {
    stop("invalid value for selected.table")
  }
  
  ## Identify name elements
  elements <- strsplit(names(tables), "_", fixed = TRUE)
  element.values <-  bdpreport.get.element.values(elements, length(setting.names))
  if (is.null(element.values)) {
    stop("inconsistent table names")
  }
  if (is.character(element.values)) {
    stop(element.values)
  }
  
  n <- report@tables + 1L
  wline <- function(txt, ind = 0L) {
    write.line(txt, report@fname, indent = ind)
  }
  
  ## Create a table of settings
  wline(c("<div class=\"figure\" id=\"tab", n, "figure\">"), 0)
  if (length(setting.names) != 0) {
    txt <- bdpreport.create.settings.table("tab", n, element.values, setting.names, elements[[selected.table]])
    if (is.null(txt)) {
      stop("missing plot file element descriptors")
    }
    wline(txt, 0)
  }
  
  ## Create tables
  for (i in 1:length(tables)) {
    txt <- c(" style=\"display:", ifelse(i == selected.table, "block", "none"), "\"")
    wline(c("<div id=\"tab", n, "_", names(tables)[i], "\"", txt, ">"), 1)
     bdp.add.table(report, tables[[i]], indent = indent, ...)
    wline("</div>\n", 1)
  }
  
  wline("</div>", 0)
  report@tables <- n
  return(report)
}

########################################################################################################################

# bdp.add.figure
#
# Generates HTML code for a figure in the specified report. A figure is a collection of images (plots), of which only
# one is visible at any given moment.
#
#  report         Report to write the text to.
#  description    Human-readable description of the figure. This must be a non-empty character vector. The
#                       elements of this vector are concatenated without a separator to form the full description.
#  report.plots   Object of type ReportPlot, or a list of such objects.
#  setting.names  List of plot file element descriptors. Every variable elements in the plot file names must be
#                       included in this list. Set this to empty list if no variable elements are present, that is, if
#                       the figure should present a single report plot.
#  selected.image Index of plot to be initially selected in the figure.

# @author Yassen Assenov

bdp.add.figure <- function(report, description, report.plots, setting.names = list(), selected.image = as.integer(1)) {
  
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!(is.character(description) && length(description) != 0)) {
    stop("invalid value for description")
  }
  if (inherits(report.plots, "ReportPlot")) {
    report.plots <- list(report.plots)
  }
  if (!(is.list(report.plots) && length(report.plots) != 0 && all(sapply(report.plots, inherits, "ReportPlot")))) {
    stop("invalid value for report.plots")
  }
  if (!(is.list(setting.names) && all(sapply(setting.names, is.vector)))) {
    stop("invalid value for setting.names")
  }
  if (is.numeric(selected.image) && all(as.integer(selected.image) == selected.image)) {
    selected.image <- as.integer(selected.image)
  }
  if (!is.integer(selected.image) && length(selected.image) == 1 && (!is.na(selected.image))) {
    stop("invalid value for selected.image")
  }
  if (selected.image < 1 || length(report.plots) < selected.image) {
    stop("invalid value for selected.image")
  }
  
  ## The following must be synchronizes with the values in report.css
  PLOT.IMAGE.WIDTH.MAX <- 850  ## maximum visible width, in pixels, of a figure plot image
  PLOT.IMAGE.HEIGHT.MAX <- 850 ## maximum visible height, in pixels, of a figure plot image
  
  iheight <- range(sapply(report.plots, slot, "height"))
  iwidth <- range(sapply(report.plots, slot, "width"))
  ilow <- range(sapply(report.plots, slot, "low.png"))
  if (iheight[1] != iheight[2] || iwidth[1] != iwidth[2] || ilow[1] != ilow[2]) {
    stop("inconsistent image sizes")
  }
  iheight <- iheight[1] * ilow[1]
  iwidth <- iwidth[1] * ilow[1]
  pdfs <- NULL
  if (all(sapply(report.plots, slot, "create.pdf"))) {
    pdfs <- unique(sapply(report.plots, slot, "dir.pdf"))
    if (length(pdfs) != 1) {
      pdfs <- NULL # inconsistent directory to store the PDF files; do not create a link
    }
  }
  ihigh <- range(sapply(report.plots, slot, "high.png"))
  ihigh <- ifelse(ihigh[1] == ihigh[2] && ihigh[1] > 0, ihigh[1], 0L)
  
  ## Get the file name elements
  fnames <- sapply(report.plots, slot, "fname")
  felements <- strsplit(fnames, "_", fixed = TRUE)
  element.values <-  bdpreport.get.element.values(felements, length(setting.names))
  if (is.null(element.values)) {
    stop("inconsistent figure file names")
  }
  if (is.character(element.values)) {
    stop(element.values)
  }
  
  fn <- report@figures + 1L
  wline <- function(txt, indent = 1) {
    write.line(txt, report@fname, indent = indent)
  }
  
  ## Create a table of settings if necessary
  if (length(setting.names) != 0) {
    txt <- bdpreport.create.settings.table("fig", fn, element.values, setting.names, felements[[selected.image]])
    if (is.null(txt)) {
      stop("missing plot file element descriptors")
    }
  } else {
    txt <- NULL
  }
  wline("<div class=\"figure\">", 0)
  if (!is.null(txt)) {
    wline(txt, 0)
  }
  
  ## Create image placeholder
  txt <- paste0(report@dir.fig, "/", fnames[selected.image], ".png")
  txt <- paste0("<img alt=\"Figure ", fn, "\" id=\"fig", fn, "image\" height=\"", iheight, "\" src=\"", txt,
                "\" width=\"", iwidth, "\" />")
  if (ihigh != 0) {
    ## Add link to a high-resolution image
    fimage.high <- ifelse(report@dir.high == report@dir.fig, "_high_resolution", "")
    fimage.high <- paste0(report@dir.high, "/", fnames[selected.image], fimage.high, ".png")
    txt <- paste0("<a href=\"", fimage.high, "\" id=\"fig", fn, "imagehigh\">", txt, "</a>")
  }
  if (iwidth > PLOT.IMAGE.WIDTH.MAX) {
    if (iheight > PLOT.IMAGE.HEIGHT.MAX) {
      txt <- paste0("<div class=\"restricted\">", txt, "</div>")
    } else {
      txt <- paste0("<div class=\"restrictedx\">", txt, "</div>")
    }
  } else if (iheight > PLOT.IMAGE.HEIGHT.MAX) {
    txt <- paste0("<div class=\"restrictedy\">", txt, "</div>")
  }
  wline(txt)
  
  ## Create figure title and description
  wline("<p class=\"caption\">")
  if (!is.null(pdfs)) {
    ## Create a link to a PDF file if applicable
    txt <- c("<span title=\"Open PDF\"><a class=\"figicon\" id=\"fig", fn, "pdf\" href=\"", report@dir.pdfs, "/",
             fnames[selected.image], ".pdf\"><img alt=\"Open PDF\" height=\"20\" onmouseout=\"pdficon(this, false)\" ",
             "onmouseover=\"pdficon(this)\" src=\"", report@dir.conf, "/pdf_inactive.png\" width=\"20\" /></a></span>")
    wline(txt, 2)
  }
  txt <- c("<a id=\"fig", fn, "caption\" onmousedown=\"showHide(\'fig", fn,
           "description\')\" onmouseout=\"doCursor()\" onmouseover=\"doCursor()\">Figure ", fn, "</a>")
  wline(txt, 2)
  wline("</p>")
  wline(c("<p class=\"figdescription\" id=\"fig", fn, "description\">", description, "</p>"))
  
  wline("</div>\n", 0)
  report@figures <- fn
  return(report)
}

########################################################################################################################

#  bdp.add.reference
#
# Adds a reference item to the given report.
#
#  report Report to add a reference item to.
#  txt    Text of the reference in the form of a non-empty character vector. The elements of this vector
#               are concatenated without a separator.

# @author Yassen Assenov

 bdp.add.reference <- function(report, txt) {
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!(is.character(txt) && length(txt) != 0)) {
    stop("invalid value for txt; expected character vector")
  }
  report@references <- union(report@references, paste(txt, collapse = ""))
  return(report)
}

########################################################################################################################

#  bdp.get.reference
#
# Creates a string that points to the given reference item in the specified report.
#
#  report Report that contains the reference to be cited.
#  txt    Text of the reference in the form of a non-empty character} vector. This reference must already
#               added to the report.
# @return Citation of the reference item (including a link) in the form of a one-element character vector. If
#         the specified reference item is not found in the report, this method returns an empty string.
#

# @author Yassen Assenov

 bdp.get.reference <- function(report, txt) {
  if (!inherits(report, "Report")) {
    stop("invalid value for report")
  }
  if (!(is.character(txt) && length(txt) != 0)) {
    stop("invalid value for txt; expected character vector")
  }
  result <- which(report@references == paste(txt, collapse = ""))
  if (length(result) != 0) {
    result <- paste0("[<a href=\"#ref", result, "\">", result, "</a>]")
  } else {
    result <- ""
  }
  return(result)
}

########################################################################################################################

if (!isGeneric("off")) {
  if (exists("off") && typeof(get("off")) == "closure") {
    setGeneric("off")
  } else {
    setGeneric("off", function(.Object,...) standardGeneric("off"))
  }
}

# off-methods
#
# Performs cleanup and/or other finishing activities and closes the specified device, connection, or document.
#
#  .Object       Object to be closed.
#  handle.errors Flag indicating if the method should attempt to catch and process errors (e.g. I/O errors)
#                      internally. Setting this to TRUE} does not guarantee that the method never stops with an
#                      error.
#

setMethod("off", "Report",
          function(.Object) {
            return(invisible(complete.report(.Object)))
          }
)

########################################################################################################################
########################################################################################################################

# createReport
#
# Creates a new report object.
#
#  fname              Single-element character vector denoting the name of the file to contain the HTML
#                           report. If this file already exists, it will be overwritten.
#  title              Title of the report in the form of a single-element character vector.
#  page.title         Web page title. This usually appears in the web browser's window title when the report is
#                           open. If specified, this must be a vector. Note that only the first element is used.
#  authors            Optional list of authors in the form of a character} vector. This list is included in
#                           the header of the generated HTML file. Note that author names can contain only Latin leters,
#                           space, dash (-), comma (,) or dot (.).
#  dirs               Location of the supporting directories, that is, paths that are expected to contain
#                           additional files linked to from the HTML report. See the Details section for a list
#                           of these directories.
#  init.configuration Flag indicating if the report configuration data should be initialized. If this parameter
#                           is TRUE, the method creates the respective directory and copies configuration files
#                           that define cascading style sheet (CSS) definitions and Javascript functions used by the
#                           HTML report. If such configuration files already exist, they will be overwritten. Since the
#                           aforementioned files can be shared by multiple reports, it is recommended that the
#                           configuration is initialized using the method bdp.initialize.reports, instead
#                           of setting this flag to TRUE.
# 
#
# @details
# If specified, the parameter dirs must be a character vector. The following names are read:
# 
#  "configuration" -- Directory that contains the auxilliary configuration files, such as style sheets and
#        Javascript files. If missing or NA}, the default value used is "configuration".
#  "data" -- Directory to contain the tables, lists and other generated data files that are linked to in
#        the HTML report. If missing or NA, the value used is formed from the file name fname (without
#        the extension) and the suffix "_data".
#  "fig" -- Directory to contain the low resolution PNG images shown in the HTML report. If missing or
#        NA, the value used is formed from the file name fname (without the extension) and the suffix
#        "_images".
#  "pdfs" -- Directory to contain the PDF images (if such are created). If not missing or NA}, the
#        value used is formed from the file name fname (without the extension) and the suffix "_pdf".
#  "high" -- Directory to contain the high resolution PNG images (if such are created). If missing or
#        NA}, the value used is the same as the fig} directory.
# 
# Any other elements, if present, are ignored. 
# See createReportPlot for creating image files.


# @author Yassen Assenov

createReport <- function(fname, title, page.title = "BacDup report", authors = NULL, dirs = NULL, init.configuration = FALSE) {
  if (!(is.character(fname) && length(fname) == 1 && (!is.na(fname[1])))) {
    stop("invalid value for fname")
  }
  if (!(is.character(title) && length(title) == 1 && (!is.na(title[1])))) {
    stop("invalid value for title")
  }
  if (!(is.null(page.title) || is.vector(page.title))) {
    stop("invalid value for page.title")
  }
  if (!(is.null(authors) || is.character(authors))) {
    stop("invalid value for authors")
  }
  if (is.null(dirs)) {
    dirs <- character(0)
  } else {
    if (!is.character(dirs)) {
      stop("invalid value for dirs")
    }
    dnames <- names(dirs)
    if (is.null(dnames)) { dnames <- character(0) }
    dnames <- intersect(dnames, c("configuration", "data", "fig", "pdfs", "high"))
    if (length(dnames) != 0 && any(is.na(dirs[dnames]))) {
      stop("missing values for dirs")
    }
  }
  if (!parameter.is.flag(init.configuration)) {
    stop("invalid value for init.configuration; expected TRUE or FALSE")
  }
  new("Report", fname, title, page.title, authors, dirs, init.configuration[1])
}