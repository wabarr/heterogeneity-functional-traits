GEEfiles <- list.files("DATA/", 
                       pattern="heterogeneity_results", 
                       full.names = TRUE)

extractValues <- function(string, Numeric=TRUE) {
  # this function processes the data to remove formatting from javascript objets, which is output from GEE
  # for example, the data may look like {labels=31}, where only the numeric portion is needed for analysis
  # this function uses regex to remove curly brackets and everything up to and including the equals sign
  strings <- gsub("\\{.+\\=|\\}","", string)
  if(Numeric) return(as.numeric(strings))
  return(strings)
}

readDataFrames <- function(filename){
  park <- read.csv(filename, header=TRUE)
  results <- data.frame(
    name=park$NAME,
    area=extractValues(park$area_sqkm),
    WDPAid=extractValues(park$WDPAID),
    patchRichness=extractValues(park$patchRichness),
    nPatches=extractValues(park$nPatches),
    stddev_WC=extractValues(park$stddev_WC),
    mean_WC=extractValues(park$mean_WC),
    analysis_scale=extractValues(park$analysis_scale)
  )
  return(results)
}

resultDFs <- lapply(GEEfiles, FUN=readDataFrames)
results <- do.call(rbind, resultDFs)
write.table(results, file="DATA/heterogeneity_all_scales.csv", row.names = FALSE, sep=",")
