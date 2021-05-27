library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)

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


## read input table
dups_stats <- read.table(file="/home/jfsanchez/git_repos/BacDup/developer/test_R/info_annot.csv", header=TRUE, sep=',', row.names = 1)

colnames(dups_stats)

ggline_plot(head(dups_stats),10, 'Strains', 'genes')
ggline_plot(dups_stats, 10, 'Strains', 'genes')
ggline_plot(dups_stats, 5, 'Strains', 'phage_count')
ggline_plot(dups_stats, 10, 'Strains', 'total_prots')

lm_plotter(dups_stats, "transpo", "total_dupGroups", "transpo_count")
lm_plotter(dups_stats, "genes", "total_dupGroups", "genes")
lm_plotter(dups_stats, "phage_count", "total_dupGroups", "phage_count")
lm_plotter(dups_stats, "hypothetical_coun", "total_dupGroups", "hypothetical_coun")
lm_plotter(dups_stats, "pseudo", "total_dupGroups", "pseudo")

# for each list of duplicate counts
test1 <- dups_stats[dups_stats$sample=="GCF_000010765.1",]
res <- as.numeric(unlist(strsplit(test1$list_dups[1], ":")))
summary(res)

gghistogram(res)

require(scales)
ggviolin(res) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(sides="l")


ggviolin(log(res, base = 10)) + scale_y_continuous(trans = log10_trans(), labels = trans_format("log10", math_format(10^.x))) + ylim(0, 2)


ggviolin(res, trim = TRUE) + scale_y_continuous(trans = "log10") + annotation_logticks(sides="l") 
ggboxplot(res, add="jitter") + scale_y_continuous(trans = "log10") + annotation_logticks(sides="l") 


## empirical cummulative density function
ggecdf(res)






ggdotchart_plot <- function(dups_df, sep, xlabel, var_name) {
  max_entries <- length(rownames(dups_df)) 
  dups_df['id'] <- rownames(dups_df)
  dups_df = dups_df %>% select(., id, total_dupGroups, total_dups, all_of(var_name)) %>% 
    pivot_longer(., -id,names_to = "Line", values_to = "Value")
  
  print(dups_df)
  
  p <- ggdotchart(dups_df, "id", y="Value", color="Line") + 
    labs(y = "Counts") + labs(x = xlabel)
  
  print(p)
  
    return(dups_df)
}


df <- ggdotchart_plot(dups_stats,10, 'Strains', 'genes')

df[df$id=="GCF_003017805.1",]

ggdotchart_plot(dups_stats, 10, 'Strains', 'total_prots')


library(data.table)
library(RColorBrewer)

dups_stats$sample <- rownames(dups_stats)
list_dups2 <- strsplit(dups_stats$list_dups, ":")
melt_list_dups <- data.frame(id = rep(dups_stats$sample, lapply(list_dups2, length)), 
           values = as.numeric(unlist(list_dups2)))


ggplot(melt_list_dups, aes(x = id, y=values))  + 
  geom_boxplot() + 
  xlab("Strains") + 
  ylab("Duplicates / group") + theme_classic( ) + #geom_jitter(color="red", size=0.6, alpha=0.9) +
  theme(axis.text.x  = element_text(size=10, angle = 45, hjust = 1))


ggplot(melt_list_dups, aes(x = id, y=values))  + 
  geom_boxplot() + 
  xlab("Strains") + 
  ylab("Duplicates / group") + theme_classic( ) + #geom_jitter(color="red", size=0.6, alpha=0.9) +
  #theme(axis.text.x  = element_text(size=10, angle = 45, hjust = 1)) 
  coord_flip()

ggplot(melt_list_dups, aes(x = id, y=values))  + 
  geom_bin2d(bins=50, color="black") + 
  xlab("Strains") + 
  ylab("Duplicates / group") + theme_classic( ) +
  theme(axis.text.x  = element_text(size=10, angle = 45, hjust = 1)) + 
  scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(5, "RdYlBu"))((10)))) 
  

ggplot(melt_list_dups, aes(x = id, y=values))  + 
  geom_count() +
  xlab("Strains") + 
  ylab("Duplicates / group") + theme_classic( ) + #geom_jitter(color="red", size=0.6, alpha=0.9) +
  theme(axis.text.x  = element_text(size=10, angle = 45, hjust = 1))


