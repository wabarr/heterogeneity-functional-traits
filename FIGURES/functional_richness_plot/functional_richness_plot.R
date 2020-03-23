library(ggplot2)

# ntax <- 1000
# noise <- 4
# 
# high_richness <- data.frame(
#   taxon=rep(c("A","B","C", "D"), each=ntax/4),
#   value=c(
#     rnorm(250, mean=-15, sd=noise),
#     rnorm(250, mean=-5, sd=noise),
#     rnorm(250, mean= 5, sd=noise),
#     rnorm(250, mean= 15, sd=noise)
#   ),
#   rich=rep("high richness", 1000)
#   )
# 
# low_richness <- data.frame(
#   taxon=rep(c("A","B","C", "D"), each=ntax/4),
#   value=c(
#     rnorm(250, mean=-7.5, sd=noise),
#     rnorm(250, mean=-2.5, sd=noise),
#     rnorm(250, mean= 2.5, sd=noise),
#     rnorm(250, mean= 7.5, sd=noise)
#   ),
#   rich=rep("low richness", 1000)
# )
# 
# richness <- do.call(rbind, list(high_richness, low_richness))
# 
# ggplot(richness, aes(x=value, fill=taxon)) + 
#   geom_density(alpha=0.8) + 
#   theme_bw(30) + 
#   labs(y="richness", x="trait value") + 
#   facet_grid(rich~.) + 
#   theme(
#     axis.ticks.y=element_blank(), 
#     axis.ticks.x=element_blank(),
#     axis.text.y=element_blank(),
#     axis.text.x=element_blank(), 
#     legend.position = "none")
# 
# ggsave("~/Dropbox/AncientEcoDatabase/AAPA_2018_poster/functional_richness_plot/univariate_richness_plot.jpg", dpi=300, width=8, height=6)


##note...this is straight outta the vignette
##https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html


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


theme_set(theme_classic(19))
## multivariate richness

set.seed(3)
forgg <- data.frame(trait1=c(rnorm(10,sd=2), rnorm(10, sd=1)), 
           trait2=c(rnorm(10,sd=2), rnorm(10, sd=1)), 
                    FRic=rep(c("High","Low"), each=10))

fric <- ggplot(forgg, aes(trait1, trait2)) + 
  geom_point(size=3) + 
  stat_chull(linetype="dashed", color="black", alpha=0.3) + 
  theme(legend.position = "none") + 
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="Functional Richness (FRic)\nvolume in ordinated space") + 
  scale_fill_manual(values=c("#a4c045", "#045480")) + 
  scale_color_manual(values=c("#a4c045", "#045480")) + 
  scale_x_continuous(limits=c(-4, 2.5)) + 
  scale_y_continuous(limits=c(-4, 2.5))
ggsave(plot=fric, filename="~/Dropbox/AncientEcoDatabase/AAPA_2018_poster/functional_richness_plot/bivariate_richness_plot.jpg", dpi=300, width=8, height=8)


## Functional Diversity (FDiv)


centroid_trait1 <- mean(forgg$trait1)
centroid_trait2 <- mean(forgg$trait2)
forgg$dist_from_centroid <- as.matrix(dist(data.frame(trait1=c(centroid_trait1,forgg$trait1), trait2=c(centroid_trait2,forgg$trait2))))[2:(nrow(forgg)+1),1]


fdiv<-ggplot(forgg, aes(trait1, trait2)) + 
  geom_segment(data=dplyr::filter(forgg, dist_from_centroid > mean(dist_from_centroid)),
                 aes(x=trait1, 
                   y=trait2, 
                   xend=centroid_trait1, 
                   yend=centroid_trait2), 
               color="red", size=1, linetype="dotted") + 
  annotate("polygon",
           x=mean(forgg$dist_from_centroid)* cos(seq(0,2*pi,length.out=100)),
           y=mean(forgg$dist_from_centroid)* sin(seq(0,2*pi,length.out=100)), fill="white", color="#a4c045") +  
  geom_segment(data=dplyr::filter(forgg, dist_from_centroid < mean(dist_from_centroid)),
               aes(x=trait1, 
                   y=trait2, 
                   xend=centroid_trait1, 
                   yend=centroid_trait2), 
              size=1) + 
  geom_point(size=3) +
  annotate("point", x=centroid_trait1, y=centroid_trait2, size=4, color="#045480") + 
  #annotate("text", x=2.05, y=-1.05, label=paste0("bar('dC')"), parse=TRUE, color="#a4c045", size=14) + 
  scale_x_continuous(limits=c(-4, 2.5)) + 
  scale_y_continuous(limits=c(-4, 2.5)) + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) + 
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="Functional Divergence (FDiv)\nsummed deviations from mean")

ggsave(plot=fdiv, filename="~/Dropbox/AncientEcoDatabase/AAPA_2018_poster/functional_richness_plot/FDiv_plot.jpg", dpi=300, width=8, height=8)







## Evenness

tree <- ape::mst(dist(forgg[,c("trait1","trait2")]))

edgeCount <- length(which(tree==1))
edgelist <- data.frame(from=rep(NA, edgeCount), to=rep(NA, edgeCount))

counter<-1
for(i in 1:nrow(tree)){
  for(j in 1:nrow(tree)){
    if(tree[i,j]==1){
      edgelist[counter,] <- c(i, j)
      counter <- counter + 1
    }
  }}

edgelist$trait1_from <- forgg$trait1[edgelist$from]
edgelist$trait2_from <- forgg$trait2[edgelist$from]
edgelist$trait1_to <- forgg$trait1[edgelist$to]
edgelist$trait2_to <- forgg$trait2[edgelist$to]

feve <- ggplot(forgg, aes(trait1, trait2)) +
  geom_segment(edgelist, 
               mapping=aes(x=trait1_from, y=trait2_from, xend=trait1_to, yend=trait2_to), 
               color="grey70") + 
  geom_point(size=3) + 
  scale_x_continuous(limits=c(-4, 2.5)) + 
  scale_y_continuous(limits=c(-4, 2.5)) + 
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="Functional Evenness (FEve)\nregularity on min. spanning tree")

ggsave(plot=feve, filename="~/Dropbox/AncientEcoDatabase/AAPA_2018_poster/functional_richness_plot/FEve_plot.jpg", dpi=300, width=8, height=8)

pdf("~/Dropbox/heterogeneity_funcric_JHE/JHE-figs/combined_fdiv_fric_feve.pdf", width=16, height=8)
gridExtra::grid.arrange(fric, fdiv,feve, ncol=3, nrow=1)
dev.off()
