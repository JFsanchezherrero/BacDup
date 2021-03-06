if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)}
if (!require("tidyr")) {
  install.packages("tidyr")
  library(tidyr)}
if (!require("ggpubr")) {
  install.packages("ggpubr")
  library(ggpubr)}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)}

## plot regression
ggplotRegression <- function (fit, Xlabel_name, Ylabel_name) {
  
  #print(fit)
  #print(Xlabel_name)
  #print(Ylabel_name)
  
  ## https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) + 
    theme(plot.title = element_text(size=10)) +
    theme_classic() +
    labs(y = Ylabel_name) + labs(x = Xlabel_name)
  
}
## linear model representation
lm_plotter <- function(data, Xlabel_name, var1, var2, Ylabel_name="Number of duplicates") {
  formula_given <- paste(var1, "~", var2 )
  fit1 <- lm(formula_given, data)
  return(ggplotRegression(fit1, Xlabel_name, Ylabel_name))
}

## plot line of counts
ggline_plot <- function(dups_stats, sep, xlabel, var_name) {
  max_entries <- length(rownames(dups_stats)) 
  dups_stats['id'] <- rownames(dups_stats)
  dups_stats = dups_stats %>% select(., id, total_dupGroups, total_dups, all_of(var_name)) %>% 
    pivot_longer(., -id,names_to = "Line", values_to = "Value")
  
  p <- ggline(dups_stats, "id", y="Value", color="Line", 
              palette = "jco", add="jitter", point.size = 1) + 
    scale_x_discrete(breaks=seq(0,max_entries, as.integer(sep))) + ylim(0, NA) + 
    labs(y = "Counts") + labs(x = xlabel)
  return(p)
}
