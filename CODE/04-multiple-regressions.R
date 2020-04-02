library(tidyverse)
library(car)
func <- read.csv(file = "DATA/full_heterogeneity_functional_dataset.csv", header=T)
func$log_nPatches <- log10(func$nPatches)
func$log_area <- log10(func$area)
func[func==-Inf] <- NA
func$cv_WC <- func$stddev_WC/func$mean_WC
func$log_patchDensity <- log10(func$nPatches/func$area)
func$log_FRic <- log10(func$FRic)

## scale all variables to be used in models

func <- func %>% 
  group_by(analysis_scale) %>% 
  mutate(
    scale_log_patchDensity=scale(log_patchDensity),
    scale_FDiv=scale(FDiv), 
    scale_FEve=scale(FEve), 
    scale_log_FRic=scale(log_FRic),
    scale_CV_wc=scale(cv_WC),
    scale_patchRichness=scale(patchRichness), 
    scale_nbsp=scale(nbsp), 
    scale_log_area=scale(log_area),
    scale_nLocomotor=scale(nLocomotor),
    scale_nTrophic=scale(nTrophic)
    ) %>%
  ungroup()


nested_func <- 
  func %>% 
  group_by(analysis_scale) %>% 
  nest() %>%
  arrange(analysis_scale)

nested_func <- mutate(nested_func, 
                      FDiv_patchDensity_model = map(data, function(theDF) lm(scale_FDiv ~ scale_log_patchDensity + scale_nbsp + scale_log_area , data=theDF)),
                      FDiv_CVwc_model =  map(data, function(theDF) lm(scale_FDiv ~ scale_CV_wc + scale_nbsp + scale_log_area , data=theDF)),
                      FDiv_patchRichness_model =  map(data, function(theDF) lm(scale_FDiv ~ scale_patchRichness + scale_nbsp + scale_log_area , data=theDF)),
                      
                      FEve_patchDensity_model = map(data, function(theDF) lm(scale_FEve ~ scale_log_patchDensity + scale_nbsp + scale_log_area , data=theDF)),
                      FEve_CVwc_model =  map(data, function(theDF) lm(scale_FEve ~ scale_CV_wc + scale_nbsp + scale_log_area , data=theDF)),
                      FEve_patchRichness_model =  map(data, function(theDF) lm(scale_FEve ~ scale_patchRichness + scale_nbsp + scale_log_area , data=theDF)),
                      
                      FRic_patchDensity_model = map(data, function(theDF) lm(scale_log_FRic ~ scale_log_patchDensity + scale_nbsp + scale_log_area , data=theDF)),
                      FRic_CVwc_model =  map(data, function(theDF) lm(scale_log_FRic ~ scale_CV_wc + scale_nbsp + scale_log_area , data=theDF)),
                      FRic_patchRichness_model =  map(data, function(theDF) lm(scale_log_FRic ~ scale_patchRichness + scale_nbsp + scale_log_area , data=theDF)),
                      
                      # FDis_patchDensity_model = map(data, function(theDF) lm(scale(FDis) ~ scale_log_patchDensity + scale_nbsp + scale_log_area , data=theDF)),
                      # FDis_CVwc_model =  map(data, function(theDF) lm(scale(FDis) ~ scale_CV_wc + scale_nbsp + scale_log_area , data=theDF)),
                      # FDis_patchRichness_model =  map(data, function(theDF) lm(scale(FDis) ~ scale_patchRichness + scale_nbsp + scale_log_area , data=theDF)),
                      
                      nLocomotor_patchDensity_model = map(data, function(theDF) lm(scale_nLocomotor ~ scale_log_patchDensity + scale_nbsp + scale_log_area , data=theDF)),
                      nLocomotor_CVwc_model =  map(data, function(theDF) lm(scale_nLocomotor ~ scale_CV_wc + scale_nbsp + scale_log_area , data=theDF)),
                      nLocomotor_patchRichness_model =  map(data, function(theDF) lm(scale_nLocomotor ~ scale_patchRichness + scale_nbsp + scale_log_area , data=theDF)),
                      
                      nTrophic_patchDensity_model = map(data, function(theDF) lm(scale_nTrophic ~ scale_log_patchDensity + scale_nbsp + scale_log_area , data=theDF)),
                      nTrophic_CVwc_model =  map(data, function(theDF) lm(scale_nTrophic ~ scale_CV_wc + scale_nbsp + scale_log_area , data=theDF)),
                      nTrophic_patchRichness_model =  map(data, function(theDF) lm(scale_nTrophic ~ scale_patchRichness + scale_nbsp + scale_log_area , data=theDF))
                      )


extract_coefs <- function(mod) {
  #function to extract predictor coef and se plus covariate (nbsp), as well as predictor coef and se
  # plus p values and overall r^2
  stopifnot("lm" %in% class(mod))
  
  #turn the call into character vector so we can extract stuff from it
  call_ch <- as.character(mod$call)[2]
  #remove nbsp part because its always last and we don't need it
  call_ch <- gsub(" \\+ scale_nbsp \\+ scale_log_area", "", call_ch)
  #remove the scale_ part
  response_pred <-  gsub("scale_","", 
                         strsplit(call_ch, " ~ ")[[1]]
                    )

  mod_summary <- summary(mod)
  var_inf <- vif(mod)
  return(
    data.frame(
      response = response_pred[1],
      predictor = response_pred[2],
      predictor_coef = mod_summary$coefficients[2,1],
      predictor_se= mod_summary$coefficients[2,2],
      nbsp_coef = mod_summary$coefficients[3,1],
      nbsp_se = mod_summary$coefficients[3,2],
      log_area_coef = mod_summary$coefficients[4,1],
      log_area_se = mod_summary$coefficients[4,2],
      predictor_p = mod_summary$coefficients[2,4],
      nbsp_p = mod_summary$coefficients[3,4], 
      log_area_p = mod_summary$coefficients[4,4],
      r2_overall = mod_summary$r.squared,
      r2_drop_nbsp = summary(update(mod, ~ . - scale_nbsp, data=mod$model))$r.squared, 
      vif=sum(vif(mod)>5)
    )
  )
  
  }

results <- lapply(1:nrow(nested_func), FUN=function(i) {
  var_col_indices <- which(!names(nested_func) %in% c("analysis_scale", "data"))
  lapply(var_col_indices, FUN=function(j) {
    result <- extract_coefs(nested_func[i,j][[1]][[1]])
    result$analysis_scale <- nested_func[i,"analysis_scale"][[1]]
    return(result)
  })
})

results <- lapply(results, FUN=function(x){
  return(do.call(rbind, x))
})

results <- do.call(rbind, results)

levels(results$predictor) <- c("patchDensity", "CV woody cover", "patchRichness")

transparency <- 0.4
results$predictor_alpha <- as.numeric(results$predictor_p<=0.05)
results$predictor_alpha[results$predictor_alpha==0] <- transparency
#results$nbsp_alpha <- as.numeric(results$nbsp_p<=0.05)
#results$nbsp_alpha[results$nbsp_alpha==0] <- transparency
#results$log_area_alpha <- as.numeric(results$log_area_p<=0.05)
#results$log_area_alpha[results$log_area_alpha==0] <- transparency


theme_set(theme_bw(12))


forRegParams <- gather(results, key="standardized_regression_coefficient", value="value", predictor_coef, nbsp_coef, log_area_coef)
regParams <- ggplot(forRegParams) + 
  geom_segment(data=results, aes(x=analysis_scale/1000, xend=analysis_scale/1000, y=predictor_coef - predictor_se, yend=predictor_coef + predictor_se, alpha=predictor_alpha)) + 
  geom_segment(data=results, aes(x=analysis_scale/1000, xend=analysis_scale/1000, y=nbsp_coef - nbsp_se, yend=nbsp_coef + nbsp_se,alpha=predictor_alpha)) + 
  geom_segment(data=results, aes(x=analysis_scale/1000, xend=analysis_scale/1000, y=log_area_coef - log_area_se, yend=log_area_coef + log_area_se,alpha=predictor_alpha)) + 
  geom_point(aes(x=analysis_scale/1000, y=value, shape=standardized_regression_coefficient, color=standardized_regression_coefficient), size=1.8, alpha=rep(results$predictor_alpha,3)) + 
  facet_grid(predictor~response) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  labs(x="spatial grain of analysis (km)", y="value", color="standardized regression coefficient", shape="standardized regression coefficient") + 
  scale_x_continuous(trans = "log10") + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values=c("#1b9e77", "#d95f02","#7570b3"), labels=c("log10 area","nbsp", "predictor")) + 
  scale_shape(labels=c("log10 area","nbsp", "predictor")) +  
  scale_alpha(guide="none")
ggsave("FIGURES/Fig5_regression_parameters.pdf",plot = regParams,  height=5, width=8, ,useDingbats=F)


forR2Plot <- gather(results, key="r2", value="value", r2_overall, r2_drop_nbsp)
forR2Plot$r2 <- factor(forR2Plot$r2)
levels(forR2Plot$r2) <- c("without nbsp", "overall")

r2 <- ggplot(forR2Plot) + 
  geom_histogram(aes(x=value, fill=r2)) + 
  facet_grid(predictor~response) + 
  scale_fill_manual(values=c("#d95f02","#7570b3"), name=expression(R^2)) + 
  theme(legend.position = "bottom") 
  #labs(fill=expression("r^2"),parse=TRUE)
ggsave("FIGURES/Fig6_r2_values.pdf",plot = r2,  height=5, width=8,useDingbats=F)


outputTableLocation = paste0(baseDir,"/TABLES/SOM_Table3_regression_table.md")
unlink(outputTableLocation)
write(sprintf("### SOM Table 3 - Regression Results\n\nMultiple regression results at all spatial grains. Coefficient estimates are provided, as well as standard errors for the coefficients.  P-values for each coefficient are provided.  Two R^2 values are reported, the overall R^2, and the R^2 for the regression re-run with the nbsp covariate removed.\n\n"),file=outputTableLocation, append = FALSE)
writeDF <- function(x=.x, y=.y,file=outputTableLocation) {
  write(sprintf("### Analysis spatial grain = %s meters\n\n\n", y$analysis_scale), file=file, append=TRUE)
  write("\n\n", file=file, append=TRUE)
  write(knitr::kable(x, digits=3, format="pandoc"), file=file, append=TRUE)
  write("\n\n\\n", file=file, append=TRUE)
  }

results %>% 
  group_by(analysis_scale) %>% 
  arrange(predictor, response) %>% 
  select(-predictor_alpha, -vif) %>%
  group_walk(.f=writeDF)

baseDir <- getwd()
## this relies on pandoc which must be installed on your system https://pandoc.org/
system(
  sprintf("pandoc -s --reference-doc %s %s -o %s", 
          paste0(baseDir, "/TABLES/regression_table_template.docx"),
          outputTableLocation, 
          gsub("\\.md", "\\.docx",outputTableLocation)
          )
)
unlink(outputTableLocation)

