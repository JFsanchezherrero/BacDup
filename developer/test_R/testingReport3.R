#### --------------------------------------- ####
## Load packages
#### --------------------------------------- ####
source("utilities.R")
source("report_class.R")
source("Report_methods.R")
# if error: no existing definition for function ‘off’, please, clear workspace, restart R session and run again
require(ggplot2)
source("ReportPlot_class.R")
source("ReportPlot_methods.R")

source("BacDup.R")
# source("logger.R")
# source("options.R")
#### --------------------------------------- ####

#### --------------------------------------- ####
##    Set some variables
#### --------------------------------------- ####
## set some variables with argparse
## report_name = 
## working_dir =

## set working dir and get data
report.directory <- "BacDup_report2"
#setwd("/home/jfsanchez/git_repos/BacDup/developer/test_R")

#### --------------------------------------- ####
##      Get input data
#### --------------------------------------- ####

# Input information
## from given folder read info_annot.csv in folder:
## - report/input/info.csv
input_stats <- ""

#  Dup analysis
## from given folder read info_annot.csv in folder:
## - report/dups/info_annot.csv
dups_stats <- read.table(file="./info_annot.csv", header=TRUE, sep=',', row.names = 1)
#### --------------------------------------- ####


#### --------------------------------------- ####
##        Initialize the report
#### --------------------------------------- ####

# Control: Necessary? 
if (!bdp.initialize.reports(report.directory)) {
  stop(paste("Could not initialize reports in", report.directory))
}

report <- createReport(file.path(report.directory, "BacDup_report.html"), 
                       "BacDup summary report: Gene duplication analysis report")
#### --------------------------------------- ####

#### --------------------------------------- ####
## Define references: 
#### --------------------------------------- ####
## It may be possible to add in Bibtex format: check
## - add previous Bacterial duplications publications 
## - add link to github

rf1 <- "Autores. Titulo. Revista número (año) pá-ginas"
rf2 <- c("Autores. ", "Título. ",
         "revista número (año) pág-inas")
report <- bdp.add.reference(report, rf1)
report <- bdp.add.reference(report, rf2)
#### --------------------------------------- ####


#### --------------------------------------- ####
##        Introduction report section
#### --------------------------------------- ####
intro_description <- c("Add some useful description")
report <- bdp.add.section(report, "Introduction", intro_description)


#### --------------------------------------- ####
##        Input report section
#### --------------------------------------- ####
input_description <- c("Add some useful description")
report <- bdp.add.section(report, "Input data", input_description)

## create table with information from file: input_stats

#### --------------------------------------- ####
##        Duplicate analysis description report section
#### --------------------------------------- ####
dup_description <- c("Add some useful description")
report <- bdp.add.section(report, "Duplication analysis summary", dup_description)

## Create plots
##################################################
## Linear representation of duplicates per strain
##################################################
report.plot1 <- createReportPlot("ggline_genes", report, high.png = 200)
ggline_plot(dups_stats, 10, 'Strains', 'genes')
off(report.plot1)

report.plot2 <- createReportPlot("ggline_phage.count", report, high.png = 200)
ggline_plot(dups_stats, 10, 'Strains', 'phage_count')
off(report.plot2)

report.plot3 <- createReportPlot("ggline_total.prots", report, high.png = 200)
ggline_plot(dups_stats, 10, 'Strains', 'total_prots')
off(report.plot3)

report.plot4 <- createReportPlot("ggline_transpo", report, high.png = 200)
ggline_plot(dups_stats, 10, 'Strains', 'transpo')
off(report.plot4)

report.plot5 <- createReportPlot("ggline_hypothetical.coun", report, high.png = 200)
ggline_plot(dups_stats, 10, 'Strains', 'hypothetical_coun')
off(report.plot5)

report.plot6 <- createReportPlot("ggline_pseudo", report, high.png = 200)
ggline_plot(dups_stats, 10, 'Strains', 'pseudo')
off(report.plot6)

## regression
report.plot7 <- createReportPlot("lmplotter_genes", report, high.png = 200)
lm_plotter(dups_stats, "genes", "total_dupGroups", "genes")
off(report.plot7)

report.plot8 <- createReportPlot("lmplotter_phage.count", report, high.png = 200)
lm_plotter(dups_stats, "phage_count", "total_dupGroups", "phage_count")
off(report.plot8)

report.plot9 <- createReportPlot("lmplotter_total.prots", report, high.png = 200)
lm_plotter(dups_stats, "total_prots", "total_dupGroups", "total_prots")
off(report.plot9)

report.plot10 <- createReportPlot("lmplotter_phage.count", report, high.png = 200)
lm_plotter(dups_stats, "transpo", "total_dupGroups", "transpo")
off(report.plot10)

report.plot11 <- createReportPlot("lmplotter_hypothetical.coun", report, high.png = 200)
lm_plotter(dups_stats, "hypothetical_coun", "total_dupGroups", "hypothetical_coun")
off(report.plot11)

report.plot12 <- createReportPlot("lmplotter_pseudo", report, high.png = 200)
lm_plotter(dups_stats, "pseudo", "total_dupGroups", "pseudo")
off(report.plot12)
##################################################

##################################################
## Generate report figure
##################################################
report.plots <- list(report.plot1, report.plot2, report.plot3,
                     report.plot4, report.plot5, report.plot6,
                     report.plot7, report.plot8, report.plot9,
                     report.plot10, report.plot11, report.plot12)

setting.names <- list(
  "Plot type" = c("ggline" = "Lines", "lmplotter" = "Regression"),
  "Values to visualize" = c("genes" = "Genes", 
                            "phage.count" = "# Phages", 
                            "total.prots" = "Total proteins",
                            "hypothetical.coun" = "#hypothetical",
                            "transpo" = "transposase",
                            "pseudo" = "Pseudogenes"
                            )) 

description <- c("Add a useful description")
description <- paste(description, collapse = " ")
report <- bdp.add.figure(report, description, report.plots, setting.names)

#### --------------------------------------- ####
##            Close the report
#### --------------------------------------- ####
report <- bdp.add.section(report, "Summary", "The generation of this report was successful.")
off(report)
#logger.info("Closed the example report")

#logger.completed()

## Remove the generated files
#unlink(report.directory, recursive = TRUE)
