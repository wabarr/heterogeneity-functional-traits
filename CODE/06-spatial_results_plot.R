library(dplyr)
library(ggplot2)
library(ggsn)

theme_set(theme_bw(10))

africa <- read.table("./DATA/africa_for_ggplot.csv", sep=',', header=T)

africa_outline <- ggplot(africa, aes(x=long, y=lat)) + 
  geom_polygon(alpha=0.2) + 
  labs(x="longitude", y="latitude") + 
  annotate("text", x=-6.5, y=-25, label="N", color="grey80") + 
  annotate("text", x=-6.2, y=-20, label="|", color="grey80") + 
  annotate("point",x=-6.4, y=-19, size=2, pch=17, color="grey80") + 
  scalebar(africa, dist=1250, dist_unit="km", transform=TRUE, location = "bottomleft",height=0.04, st.size=2, box.fill=c("grey80", "white"),box.color="grey80", st.color="grey80", border.size=0.1)


sites <- read.csv("./DATA/full_heterogeneity_functional_dataset.csv") %>% filter(analysis_scale==1000, stddev_WC>0)

FRic <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=FRic), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="Functional Richness") 

#ggsave("FIGURES/FRic_allsites.pdf", width=12)

FDis <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=FDis), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="Functional Dispersion")

#ggsave("FIGURES/FDis_allsites.pdf", width=12)

FDiv <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=FDiv), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="Functional Divergence")

# #ggsave("FIGURES/FDiv_allsites.pdf", width=12)

patchDensity <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=log10(nPatches/area)), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="Patch Density (log10)", color="patchDensity") + 
  theme(legend.title = element_text(size=7))

patchRichness <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=patchRichness), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="Distinct Patch Types") + 
  theme(legend.title = element_text(size=5))

CV_wc <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=stddev_WC/mean_WC), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="CV Woody Cover") + 
  labs(color="CV_wc")

nTrophic <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=nTrophic), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="Trophic trait types") + 
  labs(color="nTrophic")

nLocomotor <- africa_outline + 
  geom_point(data=sites, aes(x=longitude, y=latitude, color=nLocomotor), size=4, alpha=0.7) + 
  coord_map() + 
  scale_color_distiller(palette = "Spectral") + 
  labs(title="Locomotor trait types") + 
  labs(color="nLocomotor") + 
  theme(legend.title = element_text(size=7))

pdf(file="./FIGURES/Fig2_spatial_plots.pdf", width = 8.5, height=11)
gridExtra::grid.arrange(FRic, FDiv, FDis, nLocomotor, nTrophic, CV_wc, patchDensity, patchRichness, ncol=2)
dev.off()
