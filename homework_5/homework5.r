suppressMessages({
 library("sleuth")
})

#set input and output dirs
datapath <- "/work/gene8940/lml38336/homework_5/kallisto" # NOTE: you need to modify this line to match the path made by your BASH script
resultdir <- "/work/gene8940/lml38336/homework_5" # NOTE: you need to modify this line to match the path made by your BASH script
setwd(resultdir)

#create a sample-to-condition metadata object
sample <- c("SRR5344681", "SRR5344682", "SRR5344683", "SRR5344684") # create vector of sample IDs
kallisto_dirs <- file.path(datapath, sample) # create vector of paths to kallisto output directories
condition <- c("WT", "WT", "dFNR", "dFNR") # create vector of treatments in same order as sample IDs
samples_to_conditions <- data.frame(sample,condition) # create dataframe associating treatments to sample IDs
samples_to_conditions <- dplyr::mutate(samples_to_conditions, path = kallisto_dirs) # add kallisto output paths to dataframe

# check that directories and metadata object are OK
# print(kallisto_dirs)
# print(samples_to_conditions)

# read data into sleuth_object, import bootstrap summaries and TPMs, and perform normalization/filtering steps
sleuth_object <- sleuth_prep(samples_to_conditions, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

# estimate parameters for the full linear model that includes the conditions as factors
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')

# estimate parameters for the reduced linear model that assumes equal transcript abundances in both conditions
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')

# perform likelihood ratio test to identify transcripts whose fit is significantly better under full model relative to reduced model
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')

# check that sleuth object is OK
models(sleuth_object)

#summarize the sleuth results for DE genes with q-values < 0.05
sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

#print the summary table for DE genes with q-values < 0.05
write.csv(x = sleuth_significant, file = "wt_vs_fnr_sleuth_q_0.05.csv", row.names = FALSE)

#visualize results for the 10 most significant DE genes
pdf(file="SleuthResults.pdf")
for(i in sleuth_significant$target_id[1:10]) {
  p1 <- plot_bootstrap(sleuth_object, i, units = "tpm", color_by = "condition")
  print(p1)
}
dev.off()

#Quit R
quit(save="no")
