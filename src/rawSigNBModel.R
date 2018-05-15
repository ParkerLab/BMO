#!/usr/bin/env Rscript

# Raw Signal NB binomial modeling

#    This will take a raw signal distribution and fit it to NB to
# calculate probabilities scores.

# INPUT:
# 1) Motifs files
# 2) Factor name
# 3) Output file handle
# 4) Column where values are
# (mu) Manual NB mean
# (size) Manual NB overdispersion
# (TP) ChIP-seq true positives

library(optparse)
library(ggplot2)
library(scales)
library(MASS)
options(stringsAsFactors = F)

option_list = list(
  # Mandatory
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="raw signal file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file handle", metavar="character"),
  make_option(c("-c", "--colnum"), type="integer", default=NULL, 
              help="raw signal column on the first file", metavar="integer"),
  
  # Optional arguments
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="motif name", metavar="character"),
  make_option(c("--dir1"), type="character", default='motifsWithBigwigCoverage/', 
              help="directory with all read counts (include trailing slash)", 
              metavar="character"),
  make_option(c("--dir2"), type="character", default='motifsWithBigwigCoverage_outsidePeaks/', 
              help="directory with all read counts outside peaks (include trailing slash)", 
              metavar="character"),
  make_option(c("--colnum2"), type="integer", default=NULL, 
              help="raw signal column on the second file", metavar="integer"),
  make_option(c("-m", "--mu"), type="numeric", default=NULL, 
              help="manual NB mean", metavar="numeric"),
  make_option(c("-s", "--size"), type="numeric", default=NULL, 
              help="manual NB overdispersion", metavar="numeric"),
  make_option(c("-t", "--tp"), type="character", default=NULL, 
              help="true positives", metavar="character"),
  make_option(c("--plotFitted"), type="logical", default=T, 
              help="make fitted distribution plots", metavar="logical"),
#   make_option(c("--qvalcutoff"), type="numeric", default=7.957526e-09, 
#               help="manual q-value cutoff", metavar="numeric"),
  make_option(c("--qvalcutoff"), type="numeric", default=1e-7, 
              help="manual q-value cutoff", metavar="numeric"),
  make_option(c("--writeBed"), action = "store_true", default = FALSE,
            help="Write bed file with adjusted p-values"),


  # Comparison to other method (optional)
  make_option(c("--otherMethod"), type="character", default=NULL, 
              help="bed file with other method scores", metavar="character"),
  make_option(c("--threshold"), type="character", default=NULL, 
              help="other method threshold", metavar="character"),
  make_option(c("--otherName"), type="character", default=NULL, 
              help="name of other method", metavar="character"),
  make_option(c("--otherCol"), type="numeric", default=NULL, 
              help="other method's column on bed file", metavar="numeric"),
  
  # Look at ChIP-seq signal and motif distances (optional)
  make_option(c("--distanceChip"), type="character", default=NULL, 
              help="name of ChIP coverage and distances file", metavar="character"),
  make_option(c("--distanceCol"), type="numeric", default=9, 
              help="ChIP coverage column on bed file", metavar="numeric"),
  make_option(c("--chipCol"), type="numeric", default=10, 
              help="Summit disctance column on bed file", metavar="numeric"),
  
  # Graphical parameters (optional)
  make_option(c("--xlim"), type="numeric", default = 300, 
            help="histograms xlim", metavar="numeric"),
  make_option(c("--ylim"), type="numeric", default = NULL, 
            help="histograms ylim", metavar="numeric"),
  make_option(c("--xlim_p1"), type="numeric", default = NULL, 
              help="p1 histogram xlim", metavar="numeric")
)

opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser)
opt <- list()

# source('/home/albanus/scripts/source/rawSigNBModel/testing.R')

if (is.null(args$file)){
  print_help(opt_parser)
  stop("At least the first three arguments must be supplied (input, output, column).", call. = FALSE)
}

# Read data
d <- read.table(paste(args$dir1, args$file, sep = ''))
d2 <- read.table(paste(args$dir2, args$file, sep = ''))
outfile <- args$output
column <- args$colnum
if(!is.null(args$colnum2)){
  column2 <- args$colnum2
} else{
  column2 <- args$colnum
}


colnames(d)[column]  <- 'readCounts'
colnames(d2)[column2] <- 'readCounts'

## Parse optional arguments
# Name
if(!is.null(args$name)){
  factor <- args$name
} else {factor <- args$file}
# Evaluate performance
if(!is.null(args$tp)){
  tp <- read.table(args$tp)
  opt$evalPerf <- TRUE
} else {opt$evalPerf <- FALSE}
# Manual NB parameters
if(!is.null(args$mu) & !is.null(args$size)){
  mu <- as.numeric(args$mu)
  size <- as.numeric(args$size)
  opt$manual <- TRUE
} else {opt$manual <- FALSE}
# Compare to other method
if(!is.null(args$otherMethod) & !is.null(args$threshold) &
  !is.null(args$otherName) & !is.null(args$otherCol)){
  opt$compare <- TRUE
} else {opt$compare <- FALSE}

# Set NAs to 0 and convert to integer
if(!is.numeric(d[,column]) | !is.numeric(d[,column])){
    warning(paste(factor, 
                  ": One or more columns are not numbers. Forcefully converting and proceeding.", 
                  sep = ''))
    d[,column] <- as.numeric(d[,column])
    d2[,column] <- as.numeric(d2[,column])
}
d[is.na(d[,column]), column]   <- 0
d[,column] <- as.integer(round(d[,column],digits = 0))

d2[is.na(d2[,column2]), column2] <- 0
d2[,column2] <- as.integer(round(d2[,column2],digits = 0))


# Fit data into a Negative Binomial (complains a bit, but does its job)
n <- 15000

if(!opt$manual){  
  dist <- list(mus = NULL, sizes = NULL, logliks = NULL)
  for(i in 1:100){
    sampled <- sample(d2[,column2], n)
    fitted <- suppressWarnings(fitdistr(sampled,"negative binomial"))
    dist$sizes[i] <- fitted$estimate[1]
    dist$mus[i] <- fitted$estimate[2]
    dist$logliks[i] <- fitted$loglik
  }; rm(i, fitted)
  
  mu <- mean(dist$mus)
  muSD <- sd(dist$mus)
  size <- mean(dist$sizes)
  sizeSD <- sd(dist$sizes)
  loglik <- mean(dist$logliks)
  loglikSD <- sd(dist$logliks)
  
  negbin <- rnbinom(n = n, size = size, mu = mu)
  
  write.table(data.frame(factor = factor, size = size, mu = mu, loglik = loglik,
                         sizeSD = sizeSD, muSD = muSD, loglikSD = loglikSD), 
              file = paste(outfile,'.nb',sep=''), quote = F, row.names = F,
              sep = '\t', col.names = T)
} else {
  negbin <- rnbinom(n = n, size = size, mu = mu)
  sampled <- sample(d2[,column2], n)
}


# Use fitted or manual data to calculate signal probabilities
d$negbin <- pnbinom(q = d[,column], size = size, mu = mu, lower.tail = F)
d$log10negbin <- -log10(d$negbin)
d$qval <- p.adjust(d$negbin, 'BH')
d$log10qval <- -log10(d$qval)
d$p_by <- p.adjust(d$negbin, 'BY')
d$log10by <- -log10(d$p_by)

if(args$writeBed){
  out_data <- subset(d, select = -c(negbin,qval,log10qval,p_by))
  out_data <- out_data[,c(1:7,9,8)] # Keep column for compatibility
  outBed = paste(outfile,'.bed.gz',sep='')
  write.table(out_data, file = gzfile(outBed), 
              quote = F, row.names = F, sep = '\t', col.names = F)
}


# Evaluate performance
if(opt$evalPerf){
  # Define TPs
  d$tp <- paste(d$V1,d$V2,d$V3,sep="_") %in% paste(tp$V1,tp$V2,tp$V3,sep="_")
  # Sort by -log10 p-value, add index and cumulative fraction of TPs (Recall)
  d <- d[order(d$log10negbin, decreasing = T),]
  ntp <- sum(d$tp)
  ntn <- nrow(d) - ntp
  d$index <- 1:nrow(d)
  d$recall <- cumsum(d$tp)/ntp
  # Calculate Precision for each of the positions (TP / (TP + FP))
  d$ntp <- cumsum(d$tp)
  d$nfp <- cumsum(d$tp == F)
  d$precision <- d$ntp/(d$nfp+d$ntp)
  # Calculate FPR (FP / (TN + FP))
  d$tn <- ntn - d$nfp
  d$fpr <- d$ntp/(d$tn + d$nfp)
  
  # Calculate F1 scores
  d$f1 <- 2*(d$precision * d$recall)/(d$precision + d$recall)
  d$f1[is.na(d$f1)] <- 0
  
  maxf1qval <- d$qval[d$f1 == max(d$f1)]
  
  idxFdr <- c(sum(d$qval <= 1e-13),
              sum(d$qval <= 1e-12),
              sum(d$qval <= 1e-11),
              sum(d$qval <= 1e-10),
              sum(d$qval <= 1e-9),
              sum(d$qval <= 1e-8),
              sum(d$qval <= 1e-7),
              sum(d$qval <= 1e-6),
              sum(d$qval <= 1e-5), 
              sum(d$qval <= 1e-4), 
              sum(d$qval <= 0.001), 
              sum(d$qval <= 0.01), 
              sum(d$qval <= 0.05),
              sum(d$qval <= 0.1),
              sum(d$qval <= 0.2),
              sum(d$qval <= maxf1qval),
              sum(d$qval <= args$qvalcutoff),
              sum(d$p_by < 0.05))
  idxFdr[idxFdr == 0] <- 1
  idxFdr <- data.frame(index = idxFdr,
                       nReads = d[idxFdr, column],
                       log10negbin = d$log10negbin[idxFdr],
                       percentile = 1 - (d$index[idxFdr]/nrow(d)),
                       precision = d$precision[idxFdr],
                       recall = d$recall[idxFdr],
                       fpr = d$fpr[idxFdr],
                       fdr = factor(c('1e-13','1e-12','1e-11','1e-10','1e-9','1e-8','1e-7',
                                      '1e-6','1e-5','0.01%','0.1%','1%', '5%','10%','20%',
                                      'F1 max','Custom', 'Bonferroni'), 
                                    levels = c('1e-13','1e-12','1e-11','1e-10','1e-9','1e-8',
                                               '1e-7','1e-6','1e-5','0.01%','0.1%','1%',
                                               '5%','10%','20%','F1 max','Custom','Bonferroni'),
                                    ordered = T),
                       qval = c(1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,
                                0.001,0.01,0.05,0.1,0.2,maxf1qval,args$qvalcutoff,NA)
  )
  
  write.table(idxFdr, file = paste(outfile,'.out',sep=''), quote = F, row.names = F,
              sep = '\t', col.names = T)
  
  # Calculate peak statistics (assuming d2 is the subset of d outside peaks - might not be depending
  # on which dataset is used for d2)
  d$peaks <- paste(d$V1,d$V2,d$V3,sep="_") %in% paste(d2$V1,d2$V2,d2$V3,sep="_") == F
  
  d$categories <- ifelse(d$tp & d$peak, 'ChIP + Atac +',
                  ifelse(d$tp & !d$peak, 'ChIP + Atac -',
                  ifelse(!d$tp & d$peak, 'ChIP - Atac +',
                  ifelse(!d$tp & !d$peak, 'ChIP - Atac -',
                  'error'))))
  
  categories <- as.data.frame(as.matrix(table(d$categories)))
  categories <- rbind(
    as.data.frame(as.matrix(c(
              'All' = nrow(d),
              'ChIP -' = sum(!d$tp),
              'ChIP +' = sum(d$tp), 
              'Atac -' = sum(!d$peaks),
              'Atac +' = sum(d$peaks)))),
    categories)
  
  categories <- data.frame(factor = factor,
                           category = rownames(categories),
                           counts = categories$V1,
                           percentage = categories$V1 / nrow(d))
  
  # Get the percentages of the predictions that were correct
  
  categories2 <- as.data.frame(t(as.matrix(table(d$qval <= 0.01))), row.names = 'All')
  categories2 <- rbind(categories2,
                       as.data.frame.matrix(table(d$tp, d$qval <= 0.01), row.names = c('ChIP -', 'ChIP +'))
                       )
  categories2 <- rbind(categories2,
                       as.data.frame.matrix(table(d$peaks, d$qval <= 0.01), row.names = c('Atac -', 'Atac +'))
                       )
  categories2 <- rbind(categories2,
                       as.data.frame.matrix(table(d$categories, d$qval <= 0.01)))
  
  categories2 <- data.frame(factor = factor,
                            category = rownames(categories2),
                            calledTrue = categories2[,'TRUE'], 
                            calledFalse = categories2[,'FALSE'],
                            percentTrue = categories2[,'TRUE'] /apply(categories2,1,sum)
                            )
  categories <- merge(categories,categories2)
  
  # Output
  write.table(categories, file = paste(outfile,'.categories',sep=''), quote = F, row.names = F,
              sep = '\t', col.names = T)
  
  # Set categories to factor for plotting later
  d$categories <- factor(d$categories, levels = c('ChIP + Atac +', 'ChIP + Atac -', 
                                                  'ChIP - Atac +', 'ChIP - Atac -'), ordered = T)
}


## Plot ###

# Fitted distribution to sampled data
nbins <- ifelse(max(sampled) >= max(negbin), max(sampled), max(negbin))

if(args$plotFitted){
  p1 <- ggplot(NULL, aes(sampled)) + 
    geom_histogram(aes(fill = 'Sampled data', color = NULL), alpha = .5, bins = nbins) +
    geom_histogram(aes(negbin, fill = 'Negative binomial', color = NULL), alpha = .5, bins = nbins) +
    scale_fill_manual(name='Distributions', values = c('Sampled data' = 'black',
                                                       'Negative binomial' = 'red')
    ) +
    labs(x = 'Tags in motif vicinity', y = 'Number of motifs', title = factor) +
    theme_bw()
  if(!is.null(args$xlim_p1)){
    p1 <- p1 + coord_cartesian(xlim = c(0, args$xlim_p1))
  }
  ggsave(p1,filename=paste(outfile,'.fitted.pdf',sep=''), height=5, width=6)
  
  # Histogram of p-values
  p1.2 <- ggplot(d, aes(negbin)) + 
    geom_histogram(bins = 10, fill = 'white', color = 'black') + 
    theme_bw() +
    labs(x = 'p-value', y = 'Counts', 
         title = paste(args$name, ' p-values distribution', sep ='')) +
    scale_y_continuous(labels = comma) +
    theme(plot.title = element_text(face = 'bold', margin = margin(b=20,t=10)),
          axis.title.x = element_text(margin = margin(b=10,t=20)),
          axis.title.y = element_text(margin = margin(l=10,r=20)))
  ggsave(p1.2, filename=paste(outfile,'.p-histogram.pdf',sep=''), height=5.5, width=5.5)
}


if(opt$evalPerf){
  # Make PR plots of the relationships
  p2 <- ggplot(d, aes(x = recall, y = precision)) + 
    geom_line() +
    geom_vline(aes(xintercept = recall, color = fdr, lty = fdr), idxFdr) +
    geom_hline(aes(yintercept = precision, color = fdr, lty = fdr), idxFdr) +
    geom_point(aes(x = recall, y = precision, color = fdr), idxFdr, size = 3) +
    scale_x_continuous(labels = comma) +
    ylim(c(0,1)) +
    labs(x = 'Recall', y = 'Precision', color = 'FDR Threshold', lty = 'FDR Threshold', title = factor) +
    theme_bw()
  ggsave(p2,filename=paste(outfile,'.pr.pdf',sep=''), height=5, width=6.5)
  
  p2.1 <- ggplot(d, aes(x = index, y = f1)) + 
    geom_line() +
    scale_x_continuous(labels = comma) +
    labs(x = 'Motif index', y = 'F1 score', title = factor) +
    theme_bw()
  ggsave(p2.1,filename=paste(outfile,'.f1.pdf',sep=''), height=5, width=6.5)

  # Plot overall distributions
  p3 <- ggplot(d, aes(x = index, y = log10negbin)) +
    geom_line() +
    geom_vline(aes(xintercept = index, color = fdr, lty = fdr), idxFdr, lwd = 1) +
    geom_point(aes(x = index, y = log10negbin, color = fdr), idxFdr, size = 3) +
    scale_x_continuous(labels = comma) +
    labs(x = 'Ordered Motif Rank', y = expression('-log'[10]*'(p-value)'), 
         color = 'FDR Threshold', lty = 'FDR Threshold', title = factor) +
    theme_bw()
  ggsave(p3, file = paste(outfile, '.pdf', sep = ''), width = 5.5, height = 4)
  
  # Make categories plots
  source('/home/albanus/scripts/source/rawSigNBModel/plotCategories.R')
  # Compare to other method?
  if(opt$compare){
    source('/home/albanus/scripts/source/rawSigNBModel/compareMethods.R')
  }
  # Look at ChIP-seq signal and summit distances?
  if(!is.null(args$distanceChip)){
    source('/home/albanus/scripts/source/rawSigNBModel/summitDistances.R')
  }
}
