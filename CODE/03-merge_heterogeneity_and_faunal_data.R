library(dplyr)

hetero <- read.csv("DATA/heterogeneity_all_scales.csv")
WDPAids <- read.csv("DATA/WPDA_IDs.csv")
func <- read.csv("DATA/FunctionalDiversity_COMBINED_Rowan_Lintulaakso.csv")

fullDataset <- hetero %>%
  left_join(WDPAids, by=c("WDPAid" = "WDPA_ID")) %>%
  left_join(func) %>% 
  filter(patchRichness>0)


write.table(fullDataset, "DATA/full_heterogeneity_functional_dataset.csv", sep=",",row.names = FALSE)
