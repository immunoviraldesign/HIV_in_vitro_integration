# Upload library
library(circlize)
library(grid)

# clear variables
rm(list = ls())


                  
#Set as working directory -- Rename as appropriate for file location
setwd("~/Desktop/HIV_integration/Tables/")


#read in data table
#note that excel tables need to be converted to .txt files and with the headings removed for the first three files shown below for these commands to work

df = read.table("046_insertions_directional_new_calc_complex_CIGAR.txt", header = FALSE)

df = read.table("047_insertions_directional_new_calc_complex_CIGAR.txt", header = FALSE)

df = read.table("048_insertions_directional_new_calc_complex_CIGAR.txt", header = FALSE)

df = read.table("all_samples_insertions_directional_scaled_complex_CIGAR.txt", header = TRUE)

df = read.table("concensus_integration_sites.txt", header = TRUE)


# if plotting plasmid tracks, label columns for scaled counts from single sample df: eg. 046_insertions
colnames(df) = c("vector",	"position",	"scaling_position_F",	"scaling_position_R","F_scalar", "R_scalar",
                 "integration_counts_F",	"shifted_counts_F",	"merge_counts_F",	"scaled_counts_F",
                 "log2_merge_counts_F",	"log2_merge_scaled_F",	"integration_counts_R",	"shifted_counts_R",
                 "merge_counts_R",	"scaled_counts_R",	"log2_merge_counts_R",	"log2_merge_scaled_R",
                 "added_counts",	"log2_added_counts",	"log2_added_scaled",	"plasmid_sectors",	
                 "sectors",	"sector_position",	"log2_F_minus_R", "log2_F_minus_R_scaled")

# no need to label columns from "all_samples_insertions_directional_scaled_complex_CIGAR.txt" or "concensus_integration_sites.txt"




# Step 1: Initialise the chart giving factor and x-axis.
circos.clear()


# Step 2: Set track guidelines for main plasmid tracks
circos.par("track.height" = 2, "start.degree" = 90, "gap.degree" = 0)



# -OR- Set track guidelines for plasmid features
circos.par("track.height" = 2, "start.degree" = 90, "gap.degree" = 0, "cell.padding" = c(0, 0, 0, 0))



# Step 3: Initialize
# for main plasmid tracks
circos.initialize(factors=df$vector, x=df$position)

# -OR- for main plotting plasmid features
circos.initialize(factors=df$sectors, x=df$sector_position)



# Step 4: Build the plot regions.


# Use these commands if building plots for plasmid features.  It's easiest to plot all three rings 
# then delete the ones you don't need in illustrator


circos.trackPlotRegion(factors = df$sectors, y = df$merge_counts_F, panel.fun = function(x, y) {
  }, track.height = 0.05)
circos.trackPlotRegion(factors = df$sectors, y = df$merge_counts_F, panel.fun = function(x, y) {
}, track.height = 0.05)
circos.trackPlotRegion(factors = df$sectors, y = df$merge_counts_F, panel.fun = function(x, y) {
}, track.height = 0.05)


circos.text(30, 1, "NEO/KAN prom.", sector.index = "b", track.index = 1, facing = "bending.outside", niceFacing = TRUE,  cex = 1)
circos.text(400, 1, "NEO/KAN res.", sector.index = "d", track.index = 1, facing = "bending.outside", niceFacing = TRUE,  cex = 1)
circos.text(350, 1, "origin", sector.index = "f", track.index = 1, facing = "bending.outside", niceFacing = TRUE,  cex = 1)
circos.text(1, 1, "lac prom.", sector.index = "h", track.index = 1, facing = "bending.outside", niceFacing = TRUE,  cex = 1)
circos.text(30, 1, "MCS", sector.index = "j", track.index = 1, facing = "bending.outside", niceFacing = TRUE,  cex = 1)
circos.text(70, 1, "lacZ_a", sector.index = "l", track.index = 1, facing = "bending.outside", niceFacing = TRUE,  cex = 1)





# -OR- Use these commands for building main plasmid tracks.
# Use them one at a time in the following sequence
#1. choose one command beginning with circos.trackPlotRegion
#2. choose the corresponding circos.trackLines (from step 5)
#3. repeat 1 & 2 for each step as needed if plotting multitrack circles


# for all samples composite ln-scaled values with y-axis scale
circos.trackPlotRegion(factors = df$vector, y = df$ln_added_both, panel.fun = function(x, y) {
  circos.axis(major.at = c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500))
  circos.yaxis(side = "left")
}, track.height = .6)


# for integration concensus sites with y-axis scale plotted on single track
circos.trackPlotRegion(factors = df$vector, y = df$sum, panel.fun = function(x, y) {
  circos.axis(major.at = c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500))
  circos.yaxis(side = "left")
}, track.height = .6)



# for integration concensus sites with y-axis scale plotted on multi track
circos.trackPlotRegion(factors = df$vector, y = df$sum, panel.fun = function(x, y) {
  circos.axis(major.at = c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500))
  circos.yaxis(side = "left")
}, track.height = .25)

circos.trackPlotRegion(factors = df$vector, y = df$Fcount, panel.fun = function(x, y) {
  circos.yaxis(side = "left")
}, track.height = .25)

circos.trackPlotRegion(factors = df$vector, y = df$Rcount, ylim = c(0,7), panel.fun = function(x, y) {
  circos.yaxis(side = "left")
}, track.height = .25)



# for final multitrack LN scaled values with y-axis scale
circos.trackPlotRegion(factors = df$vector, y = df$ln_added_F, panel.fun = function(x, y) {
  circos.axis(major.at = c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500))
  circos.yaxis(side = "left")
}, track.height = .25)

circos.trackPlotRegion(factors = df$vector, y = df$ln_added_R, panel.fun = function(x, y) {
  circos.yaxis(side = "left")
}, track.height = .25)

circos.trackPlotRegion(factors = df$vector, y = df$ln_F_minus_R, panel.fun = function(x, y) {
  circos.yaxis(side = "left")
}, track.height = .25)






# Step 5: Add points



# for all samples composite ln-scaled values with y-axis scale
circos.trackLines(df$vector, df$position[order(df$position)], df$ln_added_both[order(df$position)], col = rgb(0.1,0.5,0.8,0.3), lwd=2, type="h")


# for integration concensus sites with y-axis scale plotted on multitrack
circos.trackLines(df$vector, df$position[order(df$position)], df$sum[order(df$position)], col = rgb(0.1,0.5,0.8,0.3), lwd=2, type="h")


# for integration concensus sites with y-axis scale plotted on multitrack
circos.trackLines(df$vector, df$position[order(df$position)], df$sum[order(df$position)], col = rgb(0.8,0,0,0.3), lwd=2, type="h", baseline = 0)
circos.trackLines(df$vector, df$position[order(df$position)], df$Fcount[order(df$position)], col = rgb(0.1,0.5,0.8,0.3), lwd=2, type="h")
circos.trackLines(df$vector, df$position[order(df$position)], df$Rcount[order(df$position)], col = rgb(0.7,0.3,0.9,0.3), lwd=2, type="h")


# for final multitrack LN scaled values with axis -- plotted in blue and purple and red
circos.trackLines(df$vector, df$position[order(df$position)], df$ln_added_F[order(df$position)], col = rgb(0.1,0.5,0.8,0.3), lwd=2, type="h")
circos.trackLines(df$vector, df$position[order(df$position)], df$ln_added_R[order(df$position)], col = rgb(0.7,0.3,0.9,0.3), lwd=2, type="h")
circos.trackLines(df$vector, df$position[order(df$position)], df$ln_F_minus_R[order(df$position)], col = rgb(0.8,0,0,0.3), lwd=2, type="h", baseline = 0)





                  