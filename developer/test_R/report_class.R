########################################################################################################################
## Report-class.R
## created: 2012-05-08
## creator: Yassen Assenov

##### Modified by Alba Moya Garces - 2021-04-15
##### Original script https://github.com/epigen/RnBeads/blob/369833af2810371f433909ad202dedaf1fa10e5d/R/Report-class.R
## ---------------------------------------------------------------------------------------------------------------------
## Report class definition 
########################################################################################################################

## C L A S S ###########################################################################################################

# Handler of a generated HTML report. Reports are initialized using the function createReport.

# Slots:
#    fname--Name of the file that contains the HTML report.
#    dir.conf--Directory that contains configuration files; usually shared between reports.
#    dir.data--Directory that contains the generated external lists and tables.
#    dir.fig--Directory that contains the generated figure image files.
#    dir.pdfs--Directory that contains the generated figure PDF files.
#    dir.high--Directory that contains the generated high-resolution image file.
#    sections--Number of sections and subsections currently added to the report.
#    opensections--Indices of currently active section and subsections.
#    figures--Number of figures currently added to the report.
#    tables--Number of selectable tables added to the report.
#    references--List of references to be added at the end of the report.

 # Methods and Functions:
 # 
 #  rnb.get.directory--Gets the location of a given report-specific directory.
 #  rnb.add.section--Generates HTML code for a new section in the report.
 #  rnb.add.paragraph--Generates HTML code for a new paragraph in the report.
 #  rnb.add.list--Generates HTML code for a list in the report.
 #  rnb.add.table--Generates HTML code for a table in the report.
 #  rnb.add.tables--Generates HTML code for a listing of tables in the report.
 #  rnb.add.figure--Generates HTML code for a figure in the report.
 #  rnb.add.reference--Adds a reference item to the report.
 #  off--Completes the HTML report by adding a reference section (if needed),
 #        a footer notice and closing the \code{<body>} and \code{<html>} tags.

setClass("Report",
         representation(fname = "character",
                        dir.conf = "character",
                        dir.data = "character",
                        dir.fig = "character",
                        dir.pdfs = "character",
                        dir.high = "character",
                        sections = "integer",
                        opensections = "integer",
                        figures = "integer",
                        tables = "integer",
                        references = "character"),
         prototype = prototype(fname = "",
                               dir.conf = "configutation",
                               dir.data = "data",
                               dir.fig = "images",
                               dir.pdfs = "images/pdf",
                               dir.high = "images/high_resolution",
                               sections = 0L,
                               opensections = rep(0L, 3L),
                               figures = 0L,
                               tables = 0L,
                               references = character()),
         package = "RnBeads")

## TODO ## Check deprecated arguments: representation, prototype, package