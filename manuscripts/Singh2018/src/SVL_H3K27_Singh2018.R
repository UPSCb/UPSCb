#' ---
#' title: "Abundance of H3K27me3 in SVL locus"
#' author: "Katja Stojkovič, Umeå Plant Science Center, Department of Forest Genetics and Plant Physiology, Swedish University of Agricultural Sciences"
#' date: "2018/5/7"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---


#' # Setup  
#' Load the libraries
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(genomeIntervals))
suppressPackageStartupMessages(library(stringr))

#' Change "data" with the path where you have stored bam files.  
#' Read file names...
bam_files <- list.files("data",
                        pattern = "SVL-locus.bam$",
                        full.names = TRUE)

#' ...and name the files
names(bam_files) <- str_extract(bam_files, "S[0-9]{2}-104?[A,B,C]-[K,H][0-9]{1,2}")

#' Select the samples
bam_files_H3K27 <- bam_files[grepl("K27", names(bam_files))]
bam_files_H3 <- bam_files[grepl("H3", names(bam_files))]

#' Construct sample information
samples_H3K27 <- data.frame(samples = c("S23-10A-K27",  "S24-10B-K27", "S25-10C-K27",  "S26-104A-K27", "S27-104B-K27", "S28-104C-K27"),
                       time_point = rep(c("10WSD", "4WC"), each = 3))

samples_H3 <- data.frame(samples = c("S33-10A-H3", "S34-10B-H3", "S35-10C-H3", "S36-104A-H3", "S37-104B-H3", "S38-104C-H3"),
                            time_point = rep(c("10WSD", "4WC"), each = 3))

#' Factors for correcting coverage for abundance differences defined earlier in the analysis
abundance.offset <- c(0.8017467, 0.8309579, 0.7561386, 0.9972086, 1.0115031, 1.1791633)
names(abundance.offset) <- c("10A", "10B",  "10C",  "104A", "104B", "104C")

#' Choose a colour palette
pal <- brewer.pal(12, "Paired")

#' # Calculating normalised H3K27 abundance  
#' Region of interest: SVL gene (Potrx004749g03714) ± 1 kb upstream/downstream region  
#'   
#' ## H3K27
#' Read in reads from .bam files  
gA_H3K27 <- lapply(as.list(bam_files_H3K27), function(bam) {
  readGAlignmentPairs(file = bam, 
                      index = bam,
                      param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isPaired = TRUE, hasUnmappedMate = FALSE))
  )
})

#' Compute coverage  
cov_gA_H3K27 <- lapply(gA_H3K27, FUN = coverage)

# extract covarage in region of interest
cov_gA_H3K27 <- lapply(cov_gA_H3K27, function(x) {
  x[GRanges("Potrx004749:18-8754")]
})

#' Correct coverage for abundance differences defined earlier in the analysis
cov_gA_corr_H3K27 <- mapply("/", cov_gA_H3K27, abundance.offset)
# turn into data frame
df_cov_H3K27 <- as.data.frame(lapply(cov_gA_corr_H3K27, function(x) as.numeric(unlist(x, use.names = FALSE))))


#' ## H3  
#' Read in reads from .bam files  
gA_H3 <- lapply(as.list(bam_files_H3), function(bam) {
  readGAlignmentPairs(file = bam, 
                      index = bam,
                      param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isPaired = TRUE, hasUnmappedMate = FALSE))
  )
})

#' Compute coverage  
cov_gA_H3 <- lapply(gA_H3, FUN = coverage)
# extract covarage in region of interest
cov_gA_H3 <- lapply(cov_gA_H3, function(x) {
  x[GRanges("Potrx004749:18-8754")]
})

#' Correct coverage for abundance differences defined earlier in the analysis
cov_gA_corr_H3 <- mapply("/", cov_gA_H3, abundance.offset)

# turn into data frame
df_cov_H3 <- as.data.frame(lapply(cov_gA_corr_H3, function(x) as.numeric(unlist(x, use.names = FALSE))))


#' ## H3K27/H3  
#' Divide abundance of H3K27 by H3  
# use log2 and add 1, not to have problems with value 0 
df_log_H3K27 <- log2(df_cov_H3K27 + 1)
df_log_H3 <- log2(df_cov_H3 + 1)

# subtract H3 from H3K27
df_cov_ratio <- df_log_H3K27 - df_log_H3

#' Calculate mean of three replicates in each time point  
df_cov_mean_ratio <- do.call(
  cbind,
  lapply(split.data.frame(t(df_cov_ratio),
                          samples_H3K27$time_point),
         colMeans))
df_cov_mean_ratio <- as.data.frame(df_cov_mean_ratio)


#' ## correct H3K27/H3 by the offset  
#' Calculate mean of H3 replicates in each time point  
df_log_mean_H3 <- do.call(
  cbind,
  lapply(split.data.frame(t(df_log_H3),
                          samples_H3$time_point),
         colMeans))
df_log_mean_H3 <- as.data.frame(df_log_mean_H3)

#' Calculate mean of H3 of all time points
log_mean_H3_timep <- apply(df_log_mean_H3, 1, mean)
#' Calculate offset from the mean for each time point
df_log_mean_H3_offset <- log_mean_H3_timep - df_log_mean_H3

#' Subtract calculated offset from mean H3K27/H3 values
df_log_mean_ratio_offset <- df_cov_mean_ratio - df_log_mean_H3_offset

#' # Graphs  
#' ## Bar plot  
#' Split region in bins  
# Length of the region of interest is 8737 nt (SVL gene ± 1 kb upstream/downstream region). 
# For the purpose of plotting, remove 1 nt, so the region is easily divided in bins of equal size.
df_cut <- df_log_mean_ratio_offset[-1, ]
length_region <- nrow(df_cut)
bin_size <- 546
split_region <- split(df_cut, ceiling((1:length_region)/bin_size))

#' Calculate mean and standard deviation  
# mean
mean_bin <- sapply(split_region, function(bin) {
  apply(bin, 2, mean)
})

# sd
sd_bin <- sapply(split_region, function(bin) {
  apply(bin, 2, sd)
})

stat_10WSD <- t(rbind(mean = data.frame(mean_bin)["10WSD", ], sd = data.frame(sd_bin)["10WSD", ]))
stat_4WC <- t(rbind(mean = data.frame(mean_bin)["4WC", ], sd = data.frame(sd_bin)["4WC", ]))

#' Plot  
# plot 10WSD
plot_mean <- barplot(stat_10WSD[ , "mean"], 
                     col = rgb(177, 89, 40, maxColorValue = 255, alpha = 70), 
                     ylim = c(-2.5, 1.0),
                     xaxt = "n",
                     ylab = "normalised H3K27me3 abundance")
arrows(x0=plot_mean,
       y0=stat_10WSD[ , "mean"]+stat_10WSD[ , "sd"],
       y1=stat_10WSD[ , "mean"]-stat_10WSD[ , "sd"], 
       angle=90,
       code=3,
       length=0.03,
       col = rgb(177, 89, 40, maxColorValue = 255, alpha = 255))

# add 4WC to the plot
plot_mean2 <- barplot(stat_4WC[ , "mean"], 
                      col = rgb(31, 120, 180, maxColorValue = 255, alpha = 50), 
                      xaxt = "n",
                      add = TRUE)
arrows(x0=plot_mean2,
       y0=stat_4WC[ , "mean"]+stat_4WC[ , "sd"],
       y1=stat_4WC[ , "mean"]-stat_4WC[ , "sd"], 
       angle=90,
       code=3,
       length=0.03,
       col = rgb(31, 120, 180, maxColorValue = 255, alpha = 200))

# add legend
legend("topright",
           legend = c("10WSD", "4WC"),
           fill = c(rgb(177, 89, 40, maxColorValue = 255, alpha = 70), rgb(31, 120, 180, maxColorValue = 255, alpha = 50)),
           border = "black",
           bty = "n")

#' ## Box plot
boxplot(df_log_mean_ratio_offset,
        names = c("10WSD", "4WC"),
        ylim = c(-3, 2),
        ylab = "normalised H3K27me3 abundance",
        xlab = "time point")

