library(GGally)
library(tidyr)
theme_set(theme_bw(9))
func <- read.csv(file = "DATA/full_heterogeneity_functional_dataset.csv", header=T)
func$log_nPatches <- log10(func$nPatches)
func$log_area <- log10(func$area)
func$cv_WC <- func$stddev_WC/func$mean_WC
func$log_patchDensity <- log10(func$nPatches/func$area)
func$log_FRic <- log10(func$FRic)

colIndices <- which(
  colnames(func) %in% 
                c("log_area",
                  "patchRichness",
                  "log_patchDensity",
                  "cv_WC",
                  "nTrophic",
                  "nLocomotor",
                  "nbsp",
                  "FEve",
                  "log_FRic",
                  "FDiv"
                  )
                )

scale_200 <- ggpairs(
  dplyr::filter(func, analysis_scale==200), 
  columns = colIndices, 
  title="Spatial Grain = 200 m", 
  diag=list(continuous="barDiag")
  )
ggsave("./FIGURES/SOM_Fig1_plots.pdf", plot = scale_200, width=9, height=9)

scale_5000 <- ggpairs(
  dplyr::filter(func, analysis_scale==5000), 
  columns = colIndices, 
  title="Spatial Grain = 5 km",
  diag=list(continuous="barDiag")
)

ggsave("./FIGURES/Fig3_pairs_plots.pdf", plot = scale_5000, width=9, height=9)

scale_20000 <- ggpairs(
  dplyr::filter(func, analysis_scale==20000), 
  columns = colIndices, 
  title="Spatial Grain = 20 km",
  diag=list(continuous="barDiag")
)
ggsave("./FIGURES/SOM_Fig2_plots.pdf", plot = scale_20000, width=9, height=9)
