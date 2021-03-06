########################################################################################################################
## ReportPlot-class.R
## created: 2012-04-16
## creator: Yassen Assenov

##### Modified by Alba Moya Garces - 2021-04-26
##### Original script https://github.com/epigen/RnBeads/blob/369833af2810371f433909ad202dedaf1fa10e5d/R/Report-methods.R
## ---------------------------------------------------------------------------------------------------------------------
## ReportPlot class definition.
########################################################################################################################

## C L A S S ###########################################################################################################

# ReportPlot Class
#
# Information about the files created to store one generated plot in a report. Report plots are initialized using the
# function createReportPlot.
#
# Slots:

#   fname -- Relative file name. It does not include path or extension.
#   width -- Width of the image in inches.
#   height -- Height of the image in inches.
#   create.pdf -- Flag indicating if a PDF image is created.
#   low.png -- Resolution, in dots per inch, used for the figure image.
#   high.png -- Resolution, in dots per inch, used for the high-resolution image.
#   dir.pdf -- Directory that contains the generated PDF file.
#   dir.png.low -- Directory that contains the generated figure image file.
#   dir.png.high -- Directory that contains the generated high-resolution image file.
# 
#
# Methods and Functions:

# get.files -- Gets the list of all files that are planned to be generated,
#        or were already generated by the report plot.
# off -- Copies the figure to a PNG file (if needed) and closes the device
#        associated with the report plot.
# 
#
# @author Yassen Assenov

setClass("ReportPlot",
         slots= c(fname = "character", width = "numeric", height = "numeric",
                        create.pdf = "logical", low.png = "integer", high.png = "integer",
                        dir.pdf = "character", dir.png.low = "character", dir.png.high = "character"),
         prototype = prototype(fname = "temp", width = 7, height = 7,
                               create.pdf = TRUE, low.png = as.integer(100), high.png = as.integer(0),
                               dir.pdf = ".", dir.png.low = ".", dir.png.high = "."))



# ReportGgPlot Class
#
# Information about the files created to store one generated plot in a report. Report plots are initialized using the
# function createReportGgPlot. It inherits from the ReportPlot class and handling is
# analogous, except that it contains an additional slot to store a ggplot object.
#
# @section Slots:

#   ggp -- ggplot object to be printed.

# @author Fabian Mueller

setClass("ReportGgPlot",
         slots= c(ggp="ANY"),
         contains="ReportPlot",
         prototype = prototype(ggp=ggplot())
         )