


report.directory <- "BacDup_report"
#
source("utilities.R")
source("report_class.R")
source("Report_methods.R")
# if error: no existing definition for function ‘off’, please, clear workspace, restart R session and run again
require(ggplot2)
source("ReportPlot_class.R")
source("ReportPlot_methods.R")
# source("logger.R")
# source("options.R")
#
if (!bdp.initialize.reports(report.directory)) {
  stop(paste("Could not initialize reports in", report.directory))
}
## Initialize the report

report <- createReport(file.path(report.directory, "example.html"), "Example")

## Define references
rf1 <- "Autores. Titulo. Revista número (año) pá-ginas"
rf2 <- c("Autores. ", "Título. ",
         "revista número (año) pág-inas")
report <- bdp.add.reference(report, rf1)
report <- bdp.add.reference(report, rf2)

## Add report sections
txt <- c("This is example report ", bdp.get.reference(report, rf1), ". It is used in testing.")
report <- bdp.add.section(report, "Introduction", txt)
txt <- c("Some background knowledge is required ", bdp.get.reference(report, rf2), ".")
report <- bdp.add.section(report, "Background", txt, level = 2)

## Create plots


doplot <- function(type = "p", ...) {
  plot(x = c(0.4, 0.6, 0.8, 1), y = c(2, 8, 3, 9), type = type, ..., main = NA, xlab = expression(beta),
       ylab = "Measure")
}

report.plot1 <- createReportPlot("example_scatterplot_data_zoom", report, high.png = 200)
doplot(pch = 16, col = c("#000080", "#00FF00"))
off(report.plot1)

report.plot2 <- createReportPlot("example_lines_data_zoom", report, high.png = 200)
doplot(type = "l", lwd = 2, col = "#00FF00")
off(report.plot2)

report.plot3 <- createReportPlot("example_both_data_zoom", report, high.png = 200)
doplot(type = "b", lwd = 2, col = "#00FF00")
off(report.plot3)

report.plot4 <- createReportPlot("example_scatterplot_data_full", report, high.png = 200)
doplot(pch = 16, col = c("#000080", "#00FF00"), xlim = c(0, 1))
off(report.plot4)

report.plot5 <- createReportPlot("example_lines_data_full", report, high.png = 200)
doplot(type = "l", lwd = 2, col = "#00FF00", xlim = c(0, 1))
off(report.plot5)

report.plot6 <- createReportPlot("example_both_data_full", report, high.png = 200)
doplot(type = "b", lwd = 2, col = "#00FF00", xlim = c(0, 1))
off(report.plot6)

#logger.info("Generated plots")

## Generate report figure
report.plots <- list(report.plot1, report.plot2, report.plot3, report.plot4, report.plot5, report.plot6)
setting.names <- list(
  "plot type" = c("scatterplot" = "scatter plot", "lines" = "line plot", "both" = "line-and-point plot"),
  "values to visualize" = c("data" = "random data", "letters" = "random letters"),
  "methylation value range" = c("full" = "full", "zoom" = "zoomed in"))
description <- c("Example figure of four data points displayed in different plots.",
                 "The horizontal axis depicts methylation &beta; value, and the vertical axis represents a measurement.",
                 "All values and colors were selected randomly.")
description <- paste(description, collapse = " ")
report <- bdp.add.section(report, "The Figure", "Here comes <a href=\"#fig1image\">Figure 1</a>:")
report <- bdp.add.figure(report, description, report.plots, setting.names)

## Close the report
report <- bdp.add.section(report, "Summary", "The generation of this report was successful.")
off(report)
#logger.info("Closed the example report")

#logger.completed()

## Remove the generated files
#unlink(report.directory, recursive = TRUE)
