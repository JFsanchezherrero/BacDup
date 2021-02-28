#https://github.com/molevol-ub/BacterialDuplicates/blob/master/scripts/R/BioCircos_plotter.R
#https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html
## load modules
#install.packages("BioCircos")
library(BioCircos)
#install.packages("colorspace")
library(colorspace)
library(dplyr)

BioCircos()

########################
## Functions
########################
#keep items without dots version
check_dot_refName <- function(name) {
  ifelse(grepl('.', name),return(unlist(strsplit(name,"\\."))[1]),
         return(name))
  }
## grepl() function searchs for matches of a string or string vector. 
## It returns TRUE if a string contains the pattern, otherwise FALSE
##
check_seqs <- function(list_items) {
  l <- NULL
  for (i in list_items) l<-c(l, check_dot_refName(i))
  return(unlist(l))
}
##
to_list <- function(df) {
  l <- list()
  for (i in rownames(df)) l[[ check_dot_refName(i) ]] <- df[i, ]
  return(l)
}
# ########
# tracks #
# ########

duplicate_group_tracks <- function(dataF) {
  ##########################
  ## for each duplicated group
  ##########################
  # chromosome, plasmids, diferents rec_id
  ## sort by rec_id & start position first
  dataF <- dataF[order(dataF$rec_id, dataF$start ),]

  ## sequence location
  seqs_vector <- check_seqs(as.character(dataF$rec_id))
  links_chromosomes_1 = head(seqs_vector, -1)
  links_chromosomes_2 = seqs_vector[-1]

  ## start position
  starts_vector <- dataF$start
  links_pos_1 = head(starts_vector)
  links_pos_2 = starts_vector

  links_labels = paste0("dup_",dataF$dup_id)

  ##########################
  ## create a track with the links
  ##########################
  DupTrack = BioCircosLinkTrack(trackname = paste0('myDupTrack_',links_labels),

                                opacities=0.3,

                                ## link position 1
                                gene1Chromosomes = links_chromosomes_1,
                                gene1Starts = links_pos_1,
                                gene1Ends = links_pos_1 + 1, #we set the end position of each link: start+1

                                ## link position 2
                                gene2Chromosomes = links_chromosomes_2,
                                gene2Starts = links_pos_2,
                                gene2Ends = links_pos_2 + 1, #we set the end position of each link: start+1

                                ## aesthetics
                                color = '#eeb242', #dataF$color[1], ## add the color from palette generated
                                displayAxis = T,

                                ## Labels: add name but do not show in track, only when mouse hovers
                                displayLabel = F,
                                labels = links_labels,

                                ## geneXNames
                                #gene1Names = c("test"),
                                #gene2Names = c("test2"),

                                maxRadius = 0.855, width = 0.5)

  ##
  return(DupTrack)
  }
##
parse_data <- function(bed_info) {
  ## add header
  colnames(bed_info) <- c("dup_id", "rec_id", "start", "end", "locus_tag", "product", "strand")

  # # ## create more informative:
  bed_info$Gene <- 'Gene: '
  cols <- c("Gene", "locus_tag")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = "" )

  cols <- c("new_gene_name", "product")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = ": " )

  cols <- c("new_gene_name", "strand")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = "; Strand: (" )

  cols <- c("new_gene_name", "dup_id")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = "); Duplicate Group: " )

  # create color palette
  rbPal <- colorRampPalette(rainbow_hcl(length(unique(bed_info$dup_id))))
  bed_info$color <- 'test'

  for (i in 1:length(unique(bed_info$dup_id)))
    bed_info[bed_info$dup_id==i,]$color <- rbPal(length(bed_info$dup_id))[i]

  return(bed_info)
}

add_dups <- function(bed_info) {
  ##########################
  ## add duplicate information
  ##########################
  ## connection between duplicates is done using BioCircosLinkTrack
  ## and duplicate genes are show using genes_tracklist

  ##########################
  ## add as many tracks as groups
  ##########################
  ## split by group
  bed_info_split <- split(bed_info, bed_info$dup_id)

  tracks_duplicates <- NULL
  for (i in bed_info_split)
    tracks_duplicates = tracks_duplicates + duplicate_group_tracks(i)

  return(tracks_duplicates)
}
##
add_genes_strand <- function(bed_info) {
  ##########################
  ## add genes per strand
  ##########################
  bed_info_strand <- split(bed_info, bed_info$strand)

  arcs_chromosomes_genes_strand_pos = check_seqs(as.character(bed_info_strand$`pos`$rec_id))  # Chromosomes on which the arcs should be displayed
  arcs_begin_genes_strand_pos = bed_info_strand$`pos`$start
  arcs_end_genes_strand_pos = bed_info_strand$`pos`$end

  ## Positive strand genes
  gene_pos_tracklist = BioCircosArcTrack('myGenePosStrandTrack',
                                         labels = as.character(bed_info_strand$`pos`$new_gene_name),
                                         arcs_chromosomes_genes_strand_pos,
                                         arcs_begin_genes_strand_pos,
                                         arcs_end_genes_strand_pos,
                                         minRadius = 0.94, maxRadius = 0.905, colors = "#ee42d4")
  pos_background_tracklist = BioCircosBackgroundTrack("myBackgroundTrack_posGenes",
                                                      minRadius = 0.90, maxRadius = 0.945,
                                                      borderColors = "black", borderSize = 0.3,
                                                      fillColors = "#FF000000")

  arcs_chromosomes_genes_strand_neg = check_seqs(as.character(bed_info_strand$`neg`$rec_id))  # Chromosomes on which the arcs should be displayed
  arcs_begin_genes_strand_neg = bed_info_strand$`neg`$start
  arcs_end_genes_strand_neg = bed_info_strand$`neg`$end

  ## Negative strand genes
  gene_neg_tracklist = BioCircosArcTrack('myGeneNegStrandTrack',
                                         labels=as.character(bed_info_strand$`neg`$new_gene_name),
                                         arcs_chromosomes_genes_strand_neg,
                                         arcs_begin_genes_strand_neg,
                                         arcs_end_genes_strand_neg,
                                         minRadius = 0.86, maxRadius = 0.89, colors = "#42eeb2")

  neg_background_tracklist = BioCircosBackgroundTrack("myBackgroundTrack_islandphat",
                                                      minRadius = 0.855, maxRadius = 0.895,
                                                      borderColors = "black", borderSize = 0.3,
                                                      fillColors = "#FF000000")

  tracks = gene_pos_tracklist + pos_background_tracklist + gene_neg_tracklist + neg_background_tracklist
  return(tracks)
}
# ##################
# # plot BioCircos #
# ##################
create_BioCircos <- function(seq_lengths, bed_info_file) {

  ## parse information
  myGenome= to_list(seq_lengths)
  myGenome
  bed_info <- parse_data(bed_info_file)

  ###############################
  ## add tracks
  ###############################
   
  ## duplicated gene connections
  tracks <- add_dups(bed_info)

  ## duplicated genes bands
  tracks <- tracks + add_genes_strand(bed_info)
  ##########################

  ##########################
  ## plot it
  ##########################
  BioCircos(genome = myGenome,
            tracklist = tracks,

            ## Aesthetics
            genomeFillColor = c('#ee5c42', '#42d4ee', '#42d4ee', '#42d4ee', '#42d4ee'), #"Spectral",
            chrPad = 0.1, ## distance between seqs
            displayGenomeBorder = T,

            ## Ticks
            genomeTicksScale = 1e+6, #genomeTicksLen = 10, #genomeTicksTextSize = 12,
            genomeTicksDisplay = TRUE,
            genomeTicksLen = 5,
            genomeTicksColor = "#000",
            genomeTicksTextSize = "0.6em",
            genomeTicksTextColor = "#000",

            ## Labels
            genomeLabelTextSize = 12,
            genomeLabelOrientation = 60,
            genomeLabelDy = 45,
            genomeLabelDx = 0.1,

            zoom = TRUE #TEXTModuleDragEvent = TRUE,
  )
}

# ################
# ## data folder
# ################
data_folder <- "./testing/"
# ## Klebsiella pneumoniae  ##
## strain KpC11
seq_lengths_kpc11 <- read.csv(paste0(data_folder, "Kpneumoniae/strain_KpC11/length.csv"),
                       header=FALSE, row.names=1)
bed_info_file_kpc11 <- read.table(paste0(data_folder, "Kpneumoniae/strain_KpC11/dup_annot.csv"), sep=",",
                                 header=TRUE)
bed_info_file_kpc11 <- bed_info_file_kpc11[,c("dup_id", "rec_id", "start", "end", "locus_tag", "product", "strand")]
# # keep only real duplicated groups
bed_info_file_kpc11 <- bed_info_file_kpc11 %>% group_by(dup_id) %>% filter(n() >1)
# 
create_BioCircos(seq_lengths = seq_lengths_kpc11,
                 bed_info_file = bed_info_file_kpc11)

## strain KpPF25
seq_lengths_KpPF25 <- read.csv(paste0(data_folder, "Kpneumoniae/strain_KpPF25/length.csv"),
                              header=FALSE, row.names=1)
bed_info_file_KpPF25 <- read.table(paste0(data_folder, "Kpneumoniae/strain_KpPF25/dup_annot.csv"), sep=",",
                                  header=TRUE)
bed_info_file_KpPF25 <- bed_info_file_KpPF25[,c("dup_id", "rec_id", "start", "end", "locus_tag", "product", "strand")]
# # keep only real duplicated groups
bed_info_file_KpPF25 <- bed_info_file_KpPF25 %>% group_by(dup_id) %>% filter(n() >1)
# 
create_BioCircos(seq_lengths = seq_lengths_KpPF25,
                 bed_info_file = bed_info_file_KpPF25)

## strain ATCC43816
seq_lengths_ATCC43816 <- read.csv(paste0(data_folder, "Kpneumoniae/strain_ATCC43816/length.csv"),
                               header=FALSE, row.names=1)
bed_info_file_ATCC43816 <- read.table(paste0(data_folder, "Kpneumoniae/strain_ATCC43816/dup_annot.csv"), sep=",",
                                   header=TRUE)
bed_info_file_ATCC43816 <- bed_info_file_ATCC43816[,c("dup_id", "rec_id", "start", "end", "locus_tag", "product", "strand")]
# # keep only real duplicated groups
bed_info_file_ATCC43816 <- bed_info_file_ATCC43816 %>% group_by(dup_id) %>% filter(n() >1)
# 

## strain BeachRanger
seq_lengths_BeachRanger <- read.csv(paste0(data_folder, "Kpneumoniae/strain_BeachRanger/length.csv"),
                                  header=FALSE, row.names=1)
bed_info_file_BeachRanger <- read.table(paste0(data_folder, "Kpneumoniae/strain_BeachRanger/dup_annot.csv"), sep=",",
                                      header=TRUE)
bed_info_file_BeachRanger <- bed_info_file_BeachRanger[,c("dup_id", "rec_id", "start", "end", "locus_tag", "product", "strand")]
# # keep only real duplicated groups
bed_info_file_BeachRanger <- bed_info_file_BeachRanger %>% group_by(dup_id) %>% filter(n() >1)
# 
create_BioCircos(seq_lengths = seq_lengths_BeachRanger,
                 bed_info_file = bed_info_file_BeachRanger)

