setwd("~/Documents/UCLA/Bricei/ROHS")

# Read the samples file
samples <- read.table("samples", header = FALSE, colClasses = "character")

# Loop through each sample
for (i in 1:nrow(samples)) {
  sample <- gsub("Bede", "Bricei", (samples[i, 1]))
  # Construct the filename
  filename <- paste0("~/Documents/UCLA/Bricei/ROHS/", sample, "_rohsKB")
  
  # Read the data
  rohs_data <- read.delim(filename, header = FALSE)
  
  # Create a histogram
  hist_values <- rohs_data$V2 / 1000
  breaks <- 10
  x_label <- "ROHs in MB"
  
  hist_plot <- hist(hist_values, breaks = breaks, main = sample, xlab = x_label)
  
  # Save the plot in a file
  plot_filename <- paste0("hist_", sample, ".png")
  png(plot_filename)
  hist(hist_values, breaks = breaks, main = sample, xlab = x_label)
  dev.off()
}
