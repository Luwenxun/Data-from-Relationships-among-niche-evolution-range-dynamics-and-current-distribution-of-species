#Niche comparison
#Load the required packages for the analysis
library(raster)
library(ecospat)
library(ENMtools)
library("FactoMineR")
library(TPD)
#Niche comparison in G-space (using C. zawadskii and C. chanetii as an example)
#Load the required projections for the analysis
cza = raster("D:/cza.EM.current.tif")
cch = raster("D:/cch.EM.current.tif")
# Set habitat suitability value per raster cell below the maxSSS value to 0
cza_max <- calc(cza, function(x) { x[x<0.305] <- 0; return(x) })
cch_max <- calc(cch, function(x) { x[x<0.4418] <- 0; return(x) })
#Calculate niche breadth
cza_breadth <- raster.breadth(cza_max)
cch_breadth <- raster.breadth(cch_max)
#Calculate niche overlap
D<- raster.overlap(cch_max, cza_max)
#Niche comparison in E-space 
#Define the geographical backgrounds (using C. zawadskii as an example)
mcp <- function (xy) {
xy <- as.data.frame(coordinates(xy))
coords.t <- chull(xy[, 1], xy[, 2])
xy.bord <- xy[coords.t, ]
xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}
buffer.size <- 1
cza.occ<- read.delim("D:/cza-occ.csv", h=T, sep=",")
mcp.cza = mcp(cza.occ[,2:3])
mcp.cza.buffer <- gBuffer(mcp.cza, width = buffer.size)
proj4string(mcp.cza.buffer) <- proj4string(Clim_Cur)
bg.area.cza <- gIntersection(ChinaPoly, mcp.cza.buffer)
bg.raster.cza <- rasterize(bg.area.cza, Clim_Cur[[1]])
bg.points.cza <- randomPoints(mask = bg.raster.cza, n = 1000)
bg.points.cza <- cbind(bg.points.cza, extract(China_Clim, bg.points.cza))
bg.points.cza <- bg.points.cza[complete.cases(data.frame(bg.points.cza)), ]
write.csv(bg.points.cza, file = "bg-points-cza.csv")
#Niche comparison
#Calculate Fric (using C. zawadskii as an example)
#cza
cza_data<-read.delim("cza_climate.csv", header=TRUE, rowname=1, sep=",")
cza_pca <- PCA(cza_data[, 2:17], scale.unit = TRUE,  graph = TRUE)
cza.ind <- get_pca_ind(cza_pca)
cza.clim.pca <- cza.ind$coord
traits_cza <- cza.clim.pca[, 1:2]
sp_cza <- cza_data$ID
cza_sp <- TPDs(species = sp_cza, traits = traits_cza)
REND(TPDs = cza_sp)
#Calculate niche overlap
n.groups <- 9 
g.names <- c("cza", "cna", "cmo", "cch", "cma", "cor.cbs", "cor.sx", "chy.sc", "chy.sx")
g.codenames <- c("cza", "cna", "cmo", "cch", "cma", "cor.cbs", "cor.sx", "chy.sc", "chy.sx")
g.colors <- c("#FF6800", "#FFB500", "#FFE600", "#B9D400", "#2E8000", "#ABDDA4", "#9DCDD1","#5E75BC", "#00008B")
#Load the required presence and background locations for the analysis
all.spec.env <- read.delim("all.spec.env.csv", h=T, sep=",")
bg.cch <- read.delim("D:/bg-cch.csv", h=T, sep=",")
bg.cza <- read.delim("D:/bg-cza.csv", h=T, sep=",")
bg.cna <- read.delim("D:/bg-cna.csv", h=T, sep=",")
bg.cmo <- read.delim("D:/bg-cmo.csv", h=T, sep=",")
bg.cma <- read.delim("D:/bg-cma.csv", h=T, sep=",")
bg.chy.sc <- read.delim("D:/bg-chy-sc.csv", h=T, sep=",")
bg.chy.sx <- read.delim("D:/bg-chy-sx.csv", h=T, sep=",")
bg.cor.cbs <- read.delim("D:/bg-cor-cbs.csv", h=T, sep=",")
bg.cor.sx <- read.delim("D:/bg-cor-sx.csv", h=T, sep=",")
back.env <- list(bg.cza, bg.cna, bg.cmo, bg.cch, bg.cma, bg.cor.cbs, bg.cor.sx, bg.chy.sc, bg.chy.sx)
#Combine related files
all.back.env <- do.call(rbind.data.frame, back.env)
data.env <- rbind(all.spec.env, all.back.env)
#Conduct PCA
w <- c(rep(0, nrow(all.spec.env)), rep(1, nrow(all.back.env)))
pca.cal <- dudi.pca(data.env, row.w = w, center = TRUE, 
                    scale = TRUE, scannf = FALSE, nf = 2)
adtion <- cumsum(c(0, sapply(back.env, nrow)))
begnd <- nrow(all.spec.env)
scores.back <- list()
scores.spec <- list()
for (i in 1:n.groups) {
  scores.spec[[i]] <- pca.cal$li[row.sp[[i]], ]
  pos <- (begnd[1] + adtion[i] + 1) : (begnd[1] + adtion[i + 1])
  scores.back[[i]] <- pca.cal$li[pos, ]  
}
#Calculate occurence density for each species
R <- 100
z <- list()
for (i in 1:n.groups) {
  z[[i]] <- ecospat.grid.clim.dyn(total.scores.back,
                                  scores.back[[i]],
                                  scores.spec[[i]],
                                  R = R)
}
#Calculate Niche overlap and test similarity 
D <- matrix(nrow = n.groups, ncol = n.groups)
rownames(D) <- colnames(D) <- g.codenames
unfilling <- stability <- expansion <- sim <- D
for (i in 2:n.groups) {
  for (j in 1:(i - 1)) {
    x1 <- z[[i]]
    x2 <- z[[j]]
    # Niche overlap
    D[i, j] <- ecospat.niche.overlap (x1, x2, cor = TRUE)$D
    
    # Niche similarity 
    sim[i, j] <- ecospat.niche.similarity.test (x1, x2, rep,
                                                overlap.alternative = "higher")$p.D
    sim[j, i] <- ecospat.niche.similarity.test (x2, x1, rep,
                                                overlap.alternative = "higher")$p.D
    
    # Niche Expansion, Stability, and Unfilling
    index1 <- ecospat.niche.dyn.index (x1, x2, 
                                       intersection = NA)$dynamic.index.w
    index2 <- ecospat.niche.dyn.index (x2, x1,
                                       intersection = NA)$dynamic.index.w
    expansion[i, j] <- index1[1]
    stability[i, j] <- index1[2]
    unfilling[i, j] <- index1[3]
    expansion[j, i] <- index2[1]
    stability[j, i] <- index2[2]
    unfilling[j, i] <- index2[3]
  }
}
