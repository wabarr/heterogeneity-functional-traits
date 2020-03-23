library(ggplot2)

StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


theme_set(theme_classic(6))
## multivariate richness

set.seed(14)
n <- 30
forgg <- data.frame(trait1=rnorm(n,sd=1), 
           trait2=rnorm(n, sd=1), 
                    FRic=rep(c("High","Low"), each=n/2))

fric <- ggplot(forgg, aes(trait1, trait2)) + 
  geom_point(size=3) + 
  stat_chull(linetype="dashed", color="black", alpha=0.3) + 
  theme(legend.position = "none") + 
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="Functional richness (FRic)\nconvex hull volume") + 
  scale_fill_manual(values=c("#a4c045", "#045480")) + 
  scale_color_manual(values=c("#a4c045", "#045480")) 

## Functional Diversity (FDiv)


centroid_trait1 <- mean(forgg$trait1)
centroid_trait2 <- mean(forgg$trait2)
forgg$dist_from_centroid <- as.matrix(dist(data.frame(trait1=c(centroid_trait1,forgg$trait1), trait2=c(centroid_trait2,forgg$trait2))))[2:(nrow(forgg)+1),1]


fdis<-ggplot(forgg, aes(trait1, trait2)) + 
  geom_segment(aes(x=trait1, 
                   y=trait2, 
                   xend=centroid_trait1, 
                   yend=centroid_trait2), 
              size=0.2) + 
  geom_point(size=2) +
  annotate("point", x=centroid_trait1, y=centroid_trait2, size=3, color="#045480", shape=15) + 
  #annotate("text", x=2.05, y=-1.05, label=paste0("bar('dC')"), parse=TRUE, color="#a4c045", size=14) + 
  #scale_x_continuous(limits=c(-4, 2.5)) + 
  #scale_y_continuous(limits=c(-4, 2.5)) + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) + 
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="Functional dispersion (FDis)\nsummed distances from centroid\n  ")


chull_indices <- chull(forgg$trait1, forgg$trait2)
centroid_hull_trait1 <- mean(forgg$trait1[chull_indices])
centroid_hull_trait2 <- mean(forgg$trait2[chull_indices])
forgg$dist_from_hull_centroid <- as.matrix(dist(data.frame(trait1=c(centroid_hull_trait1,forgg$trait1), trait2=c(centroid_hull_trait2,forgg$trait2))))[2:(nrow(forgg)+1),1]


fdiv<-ggplot(forgg, aes(trait1, trait2)) + 
  stat_chull(linetype="dashed", color="black", alpha=0.3) +
  geom_segment(data=dplyr::filter(forgg, dist_from_hull_centroid > mean(dist_from_hull_centroid)),
               aes(x=trait1, 
                   y=trait2, 
                   xend=centroid_hull_trait1, 
                   yend=centroid_hull_trait2), 
               color="red", size=0.2, linetype="longdash") + 
  annotate("polygon",
           x=centroid_hull_trait1 + mean(forgg$dist_from_hull_centroid)* cos(seq(0,2*pi,length.out=100)),
           y=centroid_hull_trait2 + mean(forgg$dist_from_hull_centroid)* sin(seq(0,2*pi,length.out=100)), fill="white", color="#a4c045") +  
  geom_segment(data=dplyr::filter(forgg, dist_from_hull_centroid < mean(dist_from_hull_centroid)),
               aes(x=trait1, 
                   y=trait2, 
                   xend=centroid_hull_trait1, 
                   yend=centroid_hull_trait2), 
               size=0.2) + 
  geom_point(size=2) +
  annotate("point", x=centroid_hull_trait1, y=centroid_hull_trait2, size=3, color="#045480", shape=17) + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) + 
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="Functional divergence (FDiv)\nsummed deviations divided by\nsummed absolute deviations")
print(fdiv)

pdf("./FIGURES/Fig1_combined_fdiv_fric_fdis.pdf", width=6, height=3)
gridExtra::grid.arrange(fric,fdis, fdiv, ncol=3, nrow=1)
dev.off()
