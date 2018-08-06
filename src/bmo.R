options(stringsAsFactors = F)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fitdistrplus))
suppressPackageStartupMessages(library(metap))
suppressPackageStartupMessages(library(parallel))

option_list = list(
    make_option(c("--f1"), type="character", default=NULL, 
                help="ATAC-seq NB file", metavar="character"),
    make_option(c("--f2"), type="character", default=NULL, 
                help="Co-occurring motifs file", metavar="character"),
    make_option(c("-n", "--name"), type="character", default=NULL, 
                help="Motif name", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="output directory", metavar="character"),
    make_option(c("-p", "--parallel"), type="integer", default=1, 
                help="Speed up things a bit with some parallelization", metavar="integer")
)
opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

# setwd("/lab/work/albanus/islet_BMO_redo")
# args$f1 <- "negative_binomials/gm12878/ALX1_1.bed.gz"
# args$f2 <- "raw_signal/gm12878/co-occurring/ALX1_1.bed"
# args$name <- "ALX1_1_test"
# args$output <- "tmp"
# args$parallel <- 4

# Print function
sys_print <- function(x){write(x, stdout())}

# Read data
sys_print("Reading data...")
df_nb     <- read.table(args$f1)
df_motifs <- read.table(args$f2)

# Fit NB distribution on co-occurring motifs
sys_print("Fitting distributions...")
co_occurring_motifs <- df_motifs$V8 - 1 # Make it zero-based for fitting NB
fitted <- fitdist(co_occurring_motifs, 'nbinom', discrete = T, method = 'mle')
parameters_fitted <- coef(fitted)

# Get p-values
pvals <- pnbinom(co_occurring_motifs, size = parameters_fitted['size'], 
                 mu = parameters_fitted['mu'], lower.tail = F)

# Make main data table
bmo_df <- data.frame(df_nb[,1:6], 
                     tag_score   = df_motifs$V7,
                     motif_score = df_motifs$V8,
                     tag_pval    = 10 ^ (-df_nb$V9),
                     motif_pval  = pvals)

# Add a pseudo-count for zero pvalues using R's maximum precision
if(min(bmo_df$tag_pval) == 0){
    sys_print("I found p=0 in at least one instance. Adding 2e-308 to all p-values.")
    bmo_df$tag_pval <- bmo_df$tag_pval + 2e-308
}

# Combine p-values with sumz (Liptak's method)
sys_print(paste("Calculating p-values with", args$parallel, "core(s)..."))
combine_p <- function(x){
    x <- as.numeric(x)
    combined_pval <- sumz(x)$p
    return(combined_pval)
}

cl <- makeCluster(args$parallel, type = 'FORK')
bmo_df$combined_p <- parApply(cl, bmo_df[,c("tag_pval", "motif_pval")], 1, combine_p)
stopCluster(cl)

# Adjust p-values
bmo_df$adj_pval <- p.adjust(bmo_df$combined_p, method = 'BY')
called_sig <- sum(bmo_df$adj_pval <= 0.05)

# Make output bed files
sys_print("Writing output...")
bmo_df$score <- -log10(bmo_df$adj_pval)
bmo_df$bound <- as.integer(bmo_df$adj_pval <= 0.05)

cols_to_print <- c(paste0("V", 1:4), "score", "V6")
bmo_bound <- bmo_df[bmo_df$bound == 1, cols_to_print]

outfile1 <- paste0(args$output, "/all/", args$name, ".all.bed")
outfile2 <- paste0(args$output, "/bound/", args$name, ".bound.bed")

dir.create(file.path(args$output, "bound"), showWarnings = FALSE)
dir.create(file.path(args$output, "all"), showWarnings = FALSE)

write.table(bmo_df, outfile1, quote = F, col.names = F, row.names = F, sep = "\t")
write.table(bmo_bound, outfile2, quote = F, col.names = F, row.names = F, sep = "\t")

sys_print("All done!")
