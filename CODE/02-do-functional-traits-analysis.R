library(dplyr)
traits_occs <- read.table("./DATA/combined_lintulaakso_rowan_traits_occurrences.csv", header=T, sep = ",", stringsAsFactors = F)

#remove small creatures
traits_occs <- traits_occs %>% dplyr::filter(bodymass >= 500)


#simplify trophic levels
# traits_occs$trophic_rowan <- forcats::fct_collapse(traits_occs$trophic_rowan,
#                                               MF = c("MF","MF-B","MF-FG","MF-G"),
#                                               B =  c("B"),
#                                               F =  c("F", "FB", "FG", "FI", "FL"),
#                                               O =  c("OM","OM-C","OM-FL","OM-I"),
#                                               C =  c("C", "CB", "CI"),
#                                               R =  c("R"),
#                                               I =  c("I"),
#                                               G =  c("G", "G-R")
# )

## read in Kissling et al. 2014 MammalDiet traits 
mammalDIET <- read.csv("DATA/selected_mammalDiet_Traits_for_heterogenetiy_paper.csv", stringsAsFactors = FALSE)


## format data for use in FD package - traits

traits <- 
  traits_occs %>%
  left_join(mammalDIET, by=c("taxon" = "joinName")) %>%
  dplyr::select(taxon, locomotor_rowan, bodymass, Vertebrate, Invertebrate, Fish, Seed, Fruit, Root, Woody, Herbaceous, Nectar) %>%
  distinct()

convertLocomotorTraits <- function(rowID) {
  #convert categorical locomotor traits to a series of binary variables
  #operates on a single row of the traits dataframe
  row <- traits[rowID,]
  if(row$locomotor == 'AQ') {
    row$aquatic <- 1
    row$terrestrial <- 1
    row$arboreal <- 0
    row$fossorial <- 0
  }
  
  if(row$locomotor == 'T') {
    row$aquatic <- 0
    row$terrestrial <- 1
    row$arboreal <- 0
    row$fossorial <- 0
  }
  
  if(row$locomotor == 'A') {
    row$aquatic <- 0
    row$terrestrial <- 0
    row$arboreal <- 1
    row$fossorial <- 0
  }
  
  if(row$locomotor == 'F') {
    row$aquatic <- 0
    row$terrestrial <- 1
    row$arboreal <- 0
    row$fossorial <- 1
  }
  
  if(row$locomotor == 'TA') {
    row$aquatic <- 0
    row$terrestrial <- 1
    row$arboreal <- 1
    row$fossorial <- 0
  }
  return(row)
}

newTraitsDF <- do.call(rbind, lapply(1:nrow(traits),FUN=convertLocomotorTraits))
traits <-  dplyr::select(newTraitsDF, taxon, Vertebrate, Invertebrate, Fish, Seed, Fruit, Root, Woody, Herbaceous, Nectar, aquatic, terrestrial, fossorial, arboreal) 

#replace zeroes with 4s, to make an ordinal integer variable
traits$Vertebrate[traits$Vertebrate==0] <- 4
traits$Invertebrate[traits$Invertebrate==0] <- 4
traits$Fish[traits$Fish==0] <- 4
traits$Seed[traits$Seed==0] <- 4
traits$Fruit[traits$Fruit==0] <- 4
traits$Root[traits$Root==0] <- 4
traits$Woody[traits$Woody==0] <- 4
traits$Herbaceous[traits$Herbaceous==0] <- 4
traits$Nectar[traits$Nectar==0] <- 4

rm(newTraitsDF, mammalDIET)
class(traits) <- "data.frame"
row.names(traits) <- traits$taxon
traits$taxon <- NULL
#traits$locomotor_rowan <- as.factor(traits$locomotor_rowan)

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


fd <- FD::dbFD(x=traits, a=occ_DF, corr="lingoes")
fd

outputFRIC <- do.call(data.frame,fd)
outputFRIC$shortName <- rownames(outputFRIC)

## count functional types
richness <- traits_occs %>% group_by(shortName) %>% 
  summarize(
    nLocomotor = length(unique(locomotor_rowan)),
    nTrophic = length(unique(trophic_rowan)))

outputFRIC <- left_join(outputFRIC, richness)

write.table(outputFRIC, file = "~/Dropbox/heterogeneity_funcric_JHE/JHE_revise_resubmit/heterogeneity-functional-traits/DATA//FunctionalDiversity_COMBINED_Rowan_Lintulaakso.csv", row.names = F, sep=",")


