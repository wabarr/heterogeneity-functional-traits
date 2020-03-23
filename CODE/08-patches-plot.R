require(gstat)
require(ggplot2)

set.seed(7)
size<-30
range<-1000
nclasses<-3
main<-NULL


theme_set(theme_minimal(8))

xy <- expand.grid(1:size, 1:size)
names(xy) <- c('x','y')

g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, range=range, model='Exp'), nmax=20)
yy <- predict(g.dummy, newdata=xy, nsim = 1)
names(yy) <- c("x","y","category")
yy$category <- cut(yy$category,nclasses,labels=F)


lowDensity <- 
  ggplot(data=yy) + 
  geom_tile(aes(x=x, y=y, fill=factor(category))) + 
  coord_fixed() + 
  theme(legend.position = "none", 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=c("#a4c045", "#045480", "sienna")) + 
  labs(title="A: patchRichness=3, low patchDensity")
print(lowDensity)




g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, range=range/100, model='Exp'), nmax=20)
yy <- predict(g.dummy, newdata=xy, nsim = 1)
names(yy) <- c("x","y","category")
yy$category <- cut(yy$category,nclasses,labels=F)

highDensity <- 
  ggplot(data=yy) + 
  geom_tile(aes(x=x, y=y, fill=factor(category))) + 
  coord_fixed() + 
  theme(legend.position = "none", 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=c("#a4c045", "#045480", "sienna")) + 
  labs(title="B: patchRichness=3, high patchDensity")
print(highDensity)

pdf(file = "./FIGURES/Fig3_patches.pdf", width=6, height = 3)
gridExtra::grid.arrange(lowDensity, highDensity, ncol=2)
dev.off()
