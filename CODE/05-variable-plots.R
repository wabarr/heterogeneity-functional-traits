library(ggplot2)
library(tidyr)
theme_set(theme_bw(9))
func <- read.csv(file = "DATA/full_heterogeneity_functional_dataset.csv", header=T)
func$log_nPatches <- log10(func$nPatches)
func$log_area <- log10(func$area)
func$CVwc <- func$stddev_WC/func$mean_WC
func$cv_WC <- func$stddev_WC/func$mean_WC
func$log_patchDensity <- log10(func$nPatches/func$area)

colIndices <- c(15, 17, 19, 20, 21, 34, 35, 36, 37, 39)

scale_200 <- ggpairs(
  dplyr::filter(func, analysis_scale==200), 
  columns = colIndices, 
  title="Spatial Grain = 200 m"
  )
ggsave("./FIGURES/SOM_Fig1_plots.pdf", plot = scale_200, width=9, height=9)

scale_5000 <- ggpairs(
  dplyr::filter(func, analysis_scale==5000), 
  columns = colIndices, 
  title="Spatial Grain = 5 km"
)
ggsave("./FIGURES/SOM_Fig2_plots.pdf", plot = scale_5000, width=9, height=9)