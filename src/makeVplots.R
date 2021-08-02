options(stringsAsFactors = F)
library(entropy)
library(ggplot2)
library(optparse)
library(RColorBrewer)
library(scales)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))
source('src/multiplot.R')

option_list = list(
  # Mandatory
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Signal file", metavar="character"),
  make_option(c("-o", '--output'), type="character", default=NULL, 
              help="output file handle", metavar="character"),
  make_option(c("-n", '--name'), type="character", default=NULL, 
              help="Motif name", metavar="integer"),
  
  # Optional
  make_option(c("--noplot"), action = 'store_true', default = F,
              help="Don't generate plot"),
  make_option(c("--split"), action = 'store_true', default = F, 
              help="Split plotting output in different files"),
  make_option(c("--step"), type = 'integer', default = 3, 
              help="Overlap window step (default = 3)"),
  make_option(c("--randomplot"), action = 'store_true', default = F,
              help="Also generate the randomized plot"),
  make_option(c("--maxfrags"), type = 'integer', default = 500000, 
              help="Maximum number of fragments (default = 500,000)"),
  make_option(c("--maxpoints"), type = 'integer', default = 250000, 
              help="Maximum points in v-plot (default = 250,000)"),
  make_option(c("--alpha"), type = 'numeric', default = 0.025, 
              help="Plot alpha (default = 0.025)"),
  make_option(c("--size"), type = 'numeric', default = 0.025, 
              help="Plot point size (default = 0.025)"),
  make_option(c("--ylim"), type = 'numeric', default = NA, 
              help="Y-axis limit value for enrichments (default = automatic)"),
  make_option(c("--ylim2"), type = 'numeric', default = 1000, 
              help="Y-axis limit value for V-plots (default = 1,000)"),
  make_option(c("--xrange"), type = 'numeric', default = 500, 
              help="X-axis range")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

# setwd('/lab/work/albanus/residenceTime_ver3')
# args$file <- 'tmp/vsignal/lab_gm12878/by_0.05/CTCF_known2.vsignal_allRegions.gz'
# args$file <- 'tmp/vsignal/buenrostro_rep1/by_0.05/AHR_1.vsignal_allRegions.gz'
# args$maxfrags <- 50000
# args$output <- 'test_vplot'
# args$name <- 'test_vplot'

# setwd('/lab/work/albanus/gm12878/tss_vplots')
# args$file <- 'vsignal/lab_gm12878.Active_TSS.vsignal.gz'
# args$output <- 'testTSS'
# args$name <- 'tss'
# args$xrange <- 2000

input    <- args$file
outfile  <- paste(args$output, '.out', sep = '')
plot_out <- paste(args$output, '.png', sep = '')
title    <- args$name
plot     <- !args$noplot
split    <- args$split
step     <- args$step
xrange   <- args$xrange

# Define which bins should be required for the information content enrichments
lim1 <- 25        # Motif-adjacent CIE peak
lim2 <- c(50, 70) # Motif proximal CIE peaks

# Minimum fragment size to include in the analyses
lowerFragmentSizeLimit <- 41

#
## Main functions
#

read_n_trim <- function(input, minSize = lowerFragmentSizeLimit){
  # Read data
  d <- read.delim(input)
  if(nrow(d) == 0){
    stop("Parsed vsignal file is empty. Closing.")
  }
  
  # Fetch number of reads in f-VICE regions before filtering
  fvice_reads_center <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim1)
  fvice_reads_side1  <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) >= lim2[1] & 
                               abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim2[2])
  fvice_reads_raw <<- fvice_reads_center + fvice_reads_side1
  
  # Trim to our region of interest
  d <- d[d$FragmentSize >= minSize,]
  maxSize  <<- max(d$FragmentSize)
  maxCoord <- max(d$FragmentMidpointDistanceToFeatureMidpoint)
  minCoord <- min(d$FragmentMidpointDistanceToFeatureMidpoint)
  coord_span <<- length(minCoord:maxCoord)
  
  # Fetch number of reads in f-VICE regions after filtering
  fvice_reads_center <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim1)
  fvice_reads_side1  <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) >= lim2[1] & 
                               abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim2[2])
  fvice_reads_trimmed <<- fvice_reads_center + fvice_reads_side1
  
  return(d)
}

# Randomize data
randomize <- function(d){
  fragdists <- d$FragmentMidpointDistanceToFeatureMidpoint
  fragdists <- sample(fragdists, replace = F)
  d$FragmentMidpointDistanceToFeatureMidpoint <- fragdists
  return(d)
}

preprocess <- function(d){
  # Transform fragment sizes and coordinates in factors to avoid 0-count sizes not being included
  minSize = min(d$FragmentSize)
  maxSize = max(d$FragmentSize)
  coords <- d$FragmentMidpointDistanceToFeatureMidpoint
  minCoord <- min(coords)
  maxCoord <- max(coords)
  d$FragmentSize <- factor(d$FragmentSize, levels = minSize:maxSize)
  d$FragmentMidpointDistanceToFeatureMidpoint <- factor(coords, levels = minCoord:maxCoord)
  
  # Split into positions
  d2 <- list(all = d)
  positions <- unique(d$FragmentPositionRelativeToFeatureMidpoint)
  for(i in positions){
    d2[[i]] <- d[d$FragmentPositionRelativeToFeatureMidpoint %in% i,]
  }
  
  # Make random set
  d2[['random']] <- randomize(d)
  
  # Make counts matrix
  for(i in 1:length(d2)){
    d2[[i]] <- table(d2[[i]]$FragmentMidpointDistanceToFeatureMidpoint, d2[[i]]$FragmentSize)
    d2[[i]] <- as.data.frame.matrix(d2[[i]])
  }
  
  # Test
  stopifnot(identical(d2$all, with(d2, L + O + R)))
  
  # Output
  return(d2)
}

compress_matrix <- function(dat){
  # Compress matrix in 5(y-axis) and 11(x-axis) bp windows to boost signal
  # (For some reason this is transposed - Past Ricardo has some explaining to do...)
  
  # comp_mat: rows=10bp position coordinate bins; cols=1bp fragment sizes
  # with 3 bp overlap (default) on the position bins only
  backstep <- 10 - step
  nrow_overlap <<- floor((nrow(dat) - backstep + 1) / step)
  comp_mat <- matrix(nrow = nrow_overlap, ncol = ncol(dat))
  comp_mat[1,] <- as.numeric(apply(dat[1:11,], 2, sum))
  end_previous <- 11; i <- 1
  positions <- 5
  while(end_previous < coord_span){
    i <- i + 1
    start <- end_previous - backstep
    end <- start + 10
    if(end >= coord_span){
      end <- coord_span
    }
    end_previous <- end
    # print(paste(i, start, end, sep = "-"))
    comp_mat[i,] <- as.numeric(apply(dat[start:end,], 2, sum))
    positions[i] <- start + 4
  }
  positions <- positions - xrange
  rownames(comp_mat) <- positions
  # print(paste(i, start, end, sep = "-"))
  # print(nrow_overlap)
  rm(start, end, end_previous, i)
  
  # # comp_mat2: rows=binned/overlapped position coordinate bins; cols=10bp fragment sizes bins
  # # with no overlap on the fragment bins
  # ncol_overlap <<- floor(ncol(comp_mat) / 10) + 1
  # comp_mat2 <- matrix(nrow = nrow_overlap, ncol = ncol_overlap)
  # comp_mat2[,1] <- apply(comp_mat[,1:10], 1, sum)
  # end_previous <- 10; i <- 1
  # while(end_previous < ncol(comp_mat)){
  #     i <- i + 1
  #     start <- end_previous + 1
  #     end <- start + 9
  #     if(end >= ncol(comp_mat)){
  #         end <- ncol(comp_mat)
  #     }
  #     end_previous <- end
  #     # print(paste(i, start, end, sep ="-"))
  #     comp_mat2[,i] <- as.numeric(apply(comp_mat[,start:end], 1, sum))
  # }
  # rownames(comp_mat2) <- positions
  # # print(paste(i, start, end, sep ="-"))
  # # print(ncol_overlap)
  # rm(start, end, end_previous, i, positions)
  
  return(comp_mat)
}

# Normalized entropy function
normentropy <- function(x){
  maxentro  <- entropy(rep(1, length(x)))
  normentro <- (maxentro - entropy(x)) / maxentro
  # if(is.na(normentro)){
  #     normentro <- 1
  # } else { }
  return(normentro)
}

## Subsample input fragment dataframe into required number of fragments
subsample <- function(frag_df, nfrags){
  if(!is.null(frag_df$FeatureNumber)){
    if(nrow(frag_df) > nfrags){
      counts <- as.data.frame(table(frag_df$FeatureNumber))
      counts$Var1 <- as.numeric(as.character(counts$Var1))
      counts <- counts[sample(1:nrow(counts), replace = F),]
      counts$cumulative <- cumsum(counts$Freq)
      counts$keep <- counts$cumulative <= nfrags
      to_keep <- counts$Var1[counts$keep]
    } else{
      to_keep <- unique(frag_df$FeatureNumber)
    }
    subsampled_data <- filter(frag_df, FeatureNumber %in% to_keep)
  } else {
    if(nrow(frag_df) > nfrags){
      to_keep <- sample(nrow(frag_df), nfrags, replace = F)
      subsampled_data <- frag_df[to_keep,]
    } else{
      subsampled_data <- frag_df
    }
  }
  return(subsampled_data)
}

#
## Run analyses
#

## Step 1: Read and pre-process data

raw_data   <- read_n_trim(input)
usedReads  <- nrow(raw_data)

subsampled_data <- subsample(raw_data, args$maxfrags)

processed_data     <- preprocess(raw_data)
processed_data_sub <- preprocess(subsampled_data)

## Step 2: Make compressed matrices

top_motifs         <- lapply(processed_data[c('L','R','O', 'all')], compress_matrix)
top_motifs_sub     <- lapply(processed_data_sub[c('L','R','O', 'all')], compress_matrix)
shuffled_data      <- compress_matrix(processed_data$random)
shuffled_data_sub  <- compress_matrix(processed_data_sub$random)

## Step 3: Calculate entropies

entro_pipeline <- function(top_significant){
  
  # Correct for diagonal bias
  zeroes <- apply(top_significant, 2, function(x) x == 0)
  my_shuffled_data <- shuffled_data
  my_shuffled_data[zeroes] <- 0
  
  t_significant_entro <- list()
  shuffled_entro <- list()
  norm_entro <- list()
  
  # Generate the vector of information contents/enrichment
  t_significant_entro$full <- apply(top_significant, 1, normentropy) / nrow(top_significant)
  shuffled_entro$full <- apply(my_shuffled_data, 1, normentropy) / nrow(my_shuffled_data)
  norm_entro$full <- log2(t_significant_entro$full / shuffled_entro$full)
  
  # Make information content/enrichment data frame
  posinfo <- data.frame(motif          = title,
                        signif_motifs  = t_significant_entro$full,
                        shuffled       = shuffled_entro$full,
                        normalized     = norm_entro$full,
                        pos = as.numeric(rownames(top_significant)))
  
  return(posinfo)
}

for(i in c('L','R','O', 'all')){
  out_variable  <- paste('entro_data_', i, sep = '')
  out_variable2 <- paste('entro_data_sub_', i, sep = '')
  assign(out_variable, entro_pipeline(top_motifs[[i]]))
  assign(out_variable2, entro_pipeline(top_motifs_sub[[i]]))
}; rm(i, out_variable)

posinfo <- cbind(entro_data_all[,1:3],
                 all_normalized     = entro_data_all$normalized,
                 left_normalized    = entro_data_L$normalized,
                 overlap_normalized = entro_data_O$normalized,
                 right_normalized   = entro_data_R$normalized,
                 position           = entro_data_all$pos)
posinfo <- as.data.frame(posinfo)

posinfo_sub <- cbind(entro_data_sub_all[,1:3],
                     all_normalized     = entro_data_sub_all$normalized,
                     left_normalized    = entro_data_sub_L$normalized,
                     overlap_normalized = entro_data_sub_O$normalized,
                     right_normalized   = entro_data_sub_R$normalized,
                     position           = entro_data_sub_all$pos)
posinfo_sub <- as.data.frame(posinfo_sub)

# Output information content enrichments
write.table(posinfo,  file = gzfile(paste(args$output, '.posinfo.gz',sep='')), sep = "\t",
            row.names = F, quote = F)

write.table(posinfo_sub,  file = gzfile(paste(args$output, '.posinfo_sub.gz',sep='')), sep = "\t",
            row.names = F, quote = F)

## Step 4: Calculate (non-normalized) f-VICE score

# Define which bins should be required for each track
fvice_cols <- cbind(abs(posinfo$position) <= lim1,                                      # Center peak
                  abs(posinfo$position) >= lim2[1] & abs(posinfo$position) <= lim2[2])  # Side peaks

# Fetch bins from positional information matrix and process signal
fetch_data <- function(posinfo, cols_to_get){
  subset_dat <- posinfo$all_normalized[cols_to_get]
  subset_dat[is.na(subset_dat)] <- 0 # Change NA's to 0
  subset_dat[subset_dat < 0] <- 0    # Set negative values to 0
  total_info <- mean(subset_dat)     # Returns the mean of the signal in the region
  return(total_info)
}
all_center <- fetch_data(posinfo, fvice_cols[,1])
all_sides  <- fetch_data(posinfo, fvice_cols[,2])
fvice     <- all_center + all_sides

all_center_sub <- fetch_data(posinfo_sub, fvice_cols[,1])
all_sides_sub  <- fetch_data(posinfo_sub, fvice_cols[,2])
fvice_sub <- all_center_sub + all_sides_sub

# Make f-VICE output
output <- data.frame(center = all_center,
                     side1  = all_sides,
                     fvice = fvice,
                     fvice_sub = fvice_sub,
                     nreads = usedReads,
                     nreads_sub = nrow(subsampled_data),
                     fvice_reads_raw = fvice_reads_raw,
                     fvice_reads_trimmed = fvice_reads_trimmed)

if(!is.null(raw_data$FeatureNumber)){
  output$nmotifs     <- max(raw_data$FeatureNumber)
  output$nmotifs_sub <- length(unique(subsampled_data$FeatureNumber))
}

write.table(output, file = outfile, sep = '\t', row.names = F, quote = F)


## Step 5: Plot information content per region

if(plot == T){
  
  # Subsample data for plotting
  d_to_plot <- subsample(subsampled_data, args$maxpoints)
  n_kept    <- nrow(d_to_plot)
  
  # ggplot-friendly data frames
  posinfo_melt <- entro_data_sub_all %>% 
    gather('type', 'value', 2:4) %>% 
    mutate(type = factor(type, 
                         levels = c("signif_motifs","shuffled","normalized"),
                         labels = c('High signal regions', 'Random', 'CIE'),
                         ordered = T)) %>% 
    mutate(category = ifelse(type %in% c('High signal regions', 'Random'), 
                             'Normalized information', 'CIE'))
  
  # Add dummy rows for manual axis limits (temporary)
  if(!is.null(args$ylim)){
    dummy <- head(posinfo_melt[posinfo_melt$category %in% "CIE",], 2)
    dummy$pos <- 510
    dummy$value <- c(args$ylim, args$ylim * -1)
    posinfo_melt <- rbind(posinfo_melt, dummy)
  }
  
  split_entro <- rbind(cbind(entro_data_sub_L, position = 'Left'),
                       cbind(entro_data_sub_O, position = 'Overlap'),
                       cbind(entro_data_sub_R, position = 'Right'))
  
  # Make plots
  p1 <- ggplot(posinfo_melt, aes(x = pos , y = value, color = type)) + 
    geom_line(na.rm = T) + 
    geom_hline(aes(yintercept = 0, linetype = category), posinfo_melt) + 
    labs(x = 'Fragment midpoint position relative to feature midpoint', 
         y = NULL, color = NULL) +
    facet_wrap(~ category, ncol = 1, strip.position = 'left', scales = 'free_y') +
    scale_x_continuous(breaks = seq(-xrange, xrange, 50), expand = c(0.0125, 0), 
                       limits = c(-xrange-5, xrange+5)) +
    scale_color_manual(breaks = c('High signal regions', 'Random'),
                       values = c("#E41A1C","#377EB8", 'black')) +
    scale_linetype_manual(breaks = c('CIE', 'Normalized information'), values = c(2, 0)) +
    guides(linetype = F) +
    theme_bw() + 
    theme(
      strip.background = element_blank(),
      strip.text   = element_text(size = 7),
      panel.border = element_rect(colour="gray90"),
      plot.title   = element_text(vjust= 1, face = 'bold'),
      axis.title.x = element_text(vjust = -0.5, size = 7),
      axis.title.y = element_text(vjust = 1, size = 6.5),
      axis.text.x  = element_text(size = 5, angle = 90, hjust = 1, vjust = .5),
      axis.text.y  = element_text(size = 5),
      legend.position = 'bottom'
    )
  # plot(p1)
  
  p1.1 <- ggplot(split_entro, aes(x = pos , y = normalized, color = position)) + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_line(na.rm = T) + 
    labs(x = NULL, y = '\nCIE\n', color = NULL) +
    scale_x_continuous(breaks = seq(-xrange, xrange, 50), expand = c(0.0125, 0)) +
    scale_colour_manual(values = c("blue", "red", "black")) +
    theme_bw() + 
    theme(
      strip.background = element_blank(),
      strip.text   = element_text(size = 7),
      panel.border = element_rect(colour="gray90"),
      plot.title   = element_text(vjust= 1, face = 'bold'),
      axis.title.y = element_text(vjust = 1, size = 6.5),
      axis.text.x  = element_blank(),
      # axis.text.x  = element_text(angle = 90, vjust = .5, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 5),
      legend.position = 'bottom'
    )
  # plot(p1.1)
  
  p2 <- ggplot(d_to_plot, aes(x     = FragmentMidpointDistanceToFeatureMidpoint, 
                              y     = FragmentSize, 
                              color = FragmentPositionRelativeToFeatureMidpoint)) + 
    geom_point(size = args$size, alpha = args$alpha) + 
    scale_x_continuous(expand = c(0.0125, 0)) +
    scale_y_continuous(expand = c(0.0125, 0)) +
    coord_cartesian(ylim = c(lowerFragmentSizeLimit, args$ylim2)) +
    scale_colour_manual(values = c("blue", "red", "black"), guide = F) +
    labs(title = title, y = '', x = '', 
         subtitle = paste0(n_kept, ' / ', nrow(raw_data), ' points plotted (', args$maxpoints,
                           ' requested)')) + 
    theme_bw(base_size = 10) + 
    theme(
      strip.background = element_rect(fill="gray90", colour=FALSE),
      panel.border  = element_rect(colour="gray90"),
      panel.spacing = unit(1, "lines"),
      plot.title    = element_text(vjust = 1, hjust = 0.5),
      axis.title.y  = element_text(vjust = 1, margin = margin(l = 8))
    )
  # plot(p2)
  
  layout <- matrix(c(1,1, 1,1, 1,1, 1,1,
                     2,2, 2,2, 2,2), ncol = 2, byrow = T)
  
  png(file = plot_out, width = 6, height = 8, units = 'in', res = 120)
  multiplot(layout = layout, p2, p1)
  dev.off()
  
  if(split == T){
    ggsave(p1, filename = paste(args$output, '.posinfo.pdf', sep = ''), device = cairo_pdf,
           width = 4, height = 3)
    ggsave(p2, filename = paste(args$output, '.vplot.png', sep = ''), width = 3.5, height = 3)
    
  }
}

# Plot randomized data
if(args$randomplot){
  d_shuf_to_plot <- randomize(d_to_plot)
  
  d_shuf_to_plot[,3] <- as.numeric(as.character(d_shuf_to_plot[,3]))
  d_shuf_to_plot[,4] <- as.numeric(as.character(d_shuf_to_plot[,4]))
  
  p_rand <- ggplot(d_shuf_to_plot, aes(x = FragmentMidpointDistanceToFeatureMidpoint, 
                                       y = FragmentSize)) + 
    geom_point(size = args$size, alpha = args$alpha) + 
    scale_x_continuous(expand = c(0.0125, 0)) +
    scale_y_continuous(expand = c(0.0125, 0)) +
    labs(title = title, y = 'Fragment size (bp)', x = 'Fragment distance to motif center (bp)') + 
    theme_bw(base_size = 10) + 
    theme(
      strip.background = element_rect(fill="gray90", colour=FALSE),
      panel.border  = element_rect(colour="gray90"),
      panel.spacing = unit(1, "lines"),
      plot.title    = element_text(vjust=1, hjust = 0.5),
      axis.title.y  = element_text(vjust=1)
    )
  # plot(p_rand)
  ggsave(p_rand, filename = paste0(args$output, '.random.png'), width = 6, height = 5)
}
