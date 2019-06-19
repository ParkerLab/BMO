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
  make_option(c("--qvalcutoff"), type="numeric", default=1e-7, 
              help="manual q-value cutoff", metavar="numeric"),
  make_option(c("--writeBed"), action = "store_true", default = FALSE,
            help="Write bed file with adjusted p-values"),
  
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


# Fit data into a Negative Binomial
n <- 10000
if (nrow(d2) < n) {
  warning("There are less than 10,000 motif instances outside peaks.")
  n <- nrow(d2)
}

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
