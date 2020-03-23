library(dplyr)
traits_occs <- read.table("./DATA/combined_lintulaakso_rowan_traits_occurrences.csv", header=T, sep = ",", stringsAsFactors = F)

#remove small creatures
traits_occs <- traits_occs %>% dplyr::filter(bodymass >= 500)


## format data for use in FD package - traits


traits <- 
  traits_occs %>%
  dplyr::select(taxon, locomotor_rowan, trophic_rowan, bodymass) %>%
  distinct()

class(traits) <- "data.frame"
row.names(traits) <- traits$taxon
traits$taxon <- NULL
traits$locomotor_rowan <- as.factor(traits$locomotor_rowan)
traits$trophic_rowan <- as.factor(traits$trophic_rowan)
# traits$bodymass <- log10(traits$bodymass)
# 

traits$trophic_rowan <- forcats::fct_collapse(traits$trophic_rowan,
  MF = c("MF","MF-B","MF-FG","MF-G"),
  B =  c("B"),
  F =  c("F", "FB", "FG", "FI", "FL"),
  O =  c("OM","OM-C","OM-FL","OM-I"),
  C =  c("C", "CB", "CI"),
  R =  c("R"),
  I =  c("I"),
  G =  c("G", "G-R")
  )

## format data for use in FD package - occurrences


occurrences <- 
  traits_occs %>% 
  dplyr::select(shortName, taxon)
occurrences <- table(occurrences$shortName, occurrences$taxon)
## using the object of class table causes errors
occ_DF <- as.data.frame(occurrences)
names(occ_DF) <- c("site", "taxon", "present")
occ_DF <- occ_DF %>% tidyr::spread(key=taxon, value=present)
rownames(occ_DF) <- occ_DF$site
occ_DF$site <- NULL


## format data for use in FD package - reorder traits to match occurrences


traits <- traits[names(occ_DF),]


## Test for completeness, and remove traits if not to make sure it doesn't go on

if(!length(complete.cases(traits))==nrow(traits)) rm(traits)



## Calculate distance based functional diversity metrics


fd_unordered_factors <- FD::dbFD(x=traits, a=occ_DF, corr="lingoes")
fd_unordered_factors

outputFRIC <- do.call(data.frame,fd_unordered_factors)
outputFRIC$shortName <- rownames(outputFRIC)

## count functional types
richness <- traits_occs %>% group_by(shortName) %>% 
  summarize(
    nLocomotor = length(unique(locomotor_rowan)),
    nTrophic = length(unique(trophic_rowan))
  )

outputFRIC <- left_join(outputFRIC, richness)

write.table(outputFRIC, file = "~/Dropbox/heterogeneity_funcric_JHE/JHE_revise_resubmit/heterogeneity-functional-traits/DATA//FunctionalDiversity_COMBINED_Rowan_Lintulaakso.csv", row.names = F, sep=",")


