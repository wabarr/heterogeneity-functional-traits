library(tidyverse)
func <- read.csv(file = "DATA/full_heterogeneity_functional_dataset.csv", header=T)
func$log_nPatches <- log10(func$nPatches)
func$log_area <- log10(func$area)
func[func==-Inf] <- NA
func$cv_WC <- func$stddev_WC/func$mean_WC
func$log_patchDensity <- log10(func$nPatches/func$area)



nested_func <- 
  func %>% 
  group_by(analysis_scale) %>% 
  nest() %>%
  arrange(analysis_scale)
rm(func)
nested_func <- mutate(nested_func, 
                      FDiv_patchDensity_model = map(data, function(df) lm(scale(FDiv) ~ scale(log_patchDensity) + scale(nbsp) , data=df)),
                      FDiv_CVwc_model =  map(data, function(df) lm(scale(FDiv) ~ scale(cv_WC) + scale(nbsp) , data=df)),
                      FDiv_patchRichness_model =  map(data, function(df) lm(scale(FDiv) ~ scale(patchRichness) + scale(nbsp) , data=df)),
                      
                      FRic_patchDensity_model = map(data, function(df) lm(scale(FRic) ~ scale(log_patchDensity) + scale(nbsp) , data=df)),
                      FRic_CVwc_model =  map(data, function(df) lm(scale(FRic) ~ scale(cv_WC) + scale(nbsp) , data=df)),
                      FRic_patchRichness_model =  map(data, function(df) lm(scale(FRic) ~ scale(patchRichness) + scale(nbsp) , data=df)),
                      
                      FDis_patchDensity_model = map(data, function(df) lm(scale(FDis) ~ scale(log_patchDensity) + scale(nbsp) , data=df)),
                      FDis_CVwc_model =  map(data, function(df) lm(scale(FDis) ~ scale(cv_WC) + scale(nbsp) , data=df)),
                      FDis_patchRichness_model =  map(data, function(df) lm(scale(FDis) ~ scale(patchRichness) + scale(nbsp) , data=df)),
                      
                      nLocomotor_patchDensity_model = map(data, function(df) lm(scale(nLocomotor) ~ scale(log_patchDensity) + scale(nbsp) , data=df)),
                      nLocomotor_CVwc_model =  map(data, function(df) lm(scale(nLocomotor) ~ scale(cv_WC) + scale(nbsp) , data=df)),
                      nLocomotor_patchRichness_model =  map(data, function(df) lm(scale(nLocomotor) ~ scale(patchRichness) + scale(nbsp) , data=df)),
                      
                      nTrophic_patchDensity_model = map(data, function(df) lm(scale(nTrophic) ~ scale(log_patchDensity) + scale(nbsp) , data=df)),
                      nTrophic_CVwc_model =  map(data, function(df) lm(scale(nTrophic) ~ scale(cv_WC) + scale(nbsp) , data=df)),
                      nTrophic_patchRichness_model =  map(data, function(df) lm(scale(nTrophic) ~ scale(patchRichness) + scale(nbsp) , data=df))
                      )


extract_coefs <- function(mod) {
  #function to extract predictor coef and se plus covariate (nbsp), as well as predictor coef and se
  # plus p values and overall r^2
  stopifnot("lm" %in% class(mod))
  
  #turn the call into character vector so we can extract stuff from it
  call_ch <- as.character(mod$call)[2]
  #remove nbsp part because its always last and we don't need it
  call_ch <- gsub(" \\+ scale\\(nbsp\\)", "", call_ch)
  #remove the scale() part
  response_pred <- gsub("\\)", "", 
                    gsub("scale\\(","", 
                         strsplit(call_ch, " ~ ")[[1]]
                         )
                    )
  mod_summary <- summary(mod)
  return(
    data.frame(
      response = response_pred[1],
      predictor = response_pred[2],
      predictor_coef = mod_summary$coefficients[2,1],
      predictor_se= mod_summary$coefficients[2,2],
      covariate_coef = mod_summary$coefficients[3,1],
      covariate_se = mod_summary$coefficients[3,2],
      predictor_p = mod_summary$coefficients[2,4],
      covariate_p = mod_summary$coefficients[3,4], 
      overall_r2 = mod_summary$r.squared
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

results$predictor_alpha <- as.numeric(results$predictor_p<=0.05)
results$predictor_alpha[results$predictor_alpha==0] <- 0.4
results$covariate_alpha <- as.numeric(results$covariate_p<=0.05)
results$covariate_alpha[results$covariate_alpha==0] <- 0.4

theme_set(theme_bw(12))


forGG <- gather(results, key="variable", value="value", predictor_coef, covariate_coef, overall_r2)
regParams <- ggplot(forGG) + 
  geom_segment(data=results, aes(x=analysis_scale/1000, xend=analysis_scale/1000, y=predictor_coef - predictor_se, yend=predictor_coef + predictor_se, alpha=predictor_alpha)) + 
  geom_segment(data=results, aes(x=analysis_scale/1000, xend=analysis_scale/1000, y=covariate_coef - covariate_se, yend=covariate_coef + covariate_se,alpha=predictor_alpha)) + 
  geom_point(aes(x=analysis_scale/1000, y=value, shape=variable, color=variable), size=1.8, alpha=rep(results$predictor_alpha,3)) + 
  facet_grid(predictor~response) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  labs(x="spatial grain of analysis (km)", y="value") + 
  scale_x_continuous(trans = "log10") + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values=c("#1b9e77", "#d95f02","#7570b3"), labels=c("nbsp regression coef.", "overall R^2", "predictor coef.")) + 
  scale_shape(labels=c("nbsp regression coef.", "overall R^2", "predictor coef.")) +  
  scale_alpha(guide="none")
ggsave("FIGURES/Fig4_regression_parameters.pdf",plot = regParams,  height=5, width=8)


write.table(dplyr::select(results, -predictor_alpha, -covariate_alpha), file="./TABLES/regression_results_table.csv", row.names = F, sep=",")


