# 02/12/2020

# Creation of shapes (polygons) from the segmented point cloud created by AMS3D (Computree platform from the ONF)

# packages
library(lidR)
library(rgdal)
library(data.table)


#### get AMS3D segmentation

cloudpoint <- read.table(file = "segmented_point_cloud.asc", header=T)

## optional: plot the point cloud with different colors for each crown
# convert to lidar point cloud
cloudpoint <- LAS(cloudpoint)
# plot
plot(cloudpoint, color = "ID", size=4,
     colorPalette = sample(rainbow(length(unique(cloudpoint@data$ID)))))

#### rasterize

cloudpoint <- cloudpoint[,c('X','Y','Z','ID')]
names(cloudpoint) <- c("x","y","z","cluster")
# order on z, to make sure the highest point of a raster cell will give its value to the raster cell
cloudpoint <- cloudpoint[order(cloudpoint$z),]

# be careful with the resolution of the raster (will depend on the density of the point cloud)
# if the point cloud has a lower resolution, increase the raster res
r0 <- raster(extent(cloudpoint), res=c(0.5,0.5))
rasterPoly <- rasterize(cbind(cloudpoint$x, cloudpoint$y), r0, field = cloudpoint$cluster, fun='last')

# majority filter
rasterPoly2 = focal(rasterPoly, w = matrix(1, nrow = 3, ncol = 3), fun = modal)

####  conversion to polygons
boundshp <- rasterToPolygons(rasterPoly2, n = 4, dissolve = T)
boundshp2 = boundshp[boundshp$layer!=0,]

#### correction of weird shapes (with holes or disjunct parts)

# to store the maximum height of a polygon
z90 <- zmax <- c()

# for each polygon:
for (j in 1:length(boundshp2$layer))
{
  new_poly = boundshp2[j,1]
  for (k in 1:length(new_poly@polygons[[1]]@Polygons))
  {
    
    # remove holes
    new_sp_nohole = SpatialPolygons( list( Polygons( list( new_poly@polygons[[1]]@Polygons[[k]] ),
                                                     ID = 1)))
    if (k == 1)
    {allnohole = new_sp_nohole}
    else
    {allnohole = rgeos::gUnion(allnohole,new_sp_nohole)}
  }
  
  # remove small parts of the crowns that are not connected to others
  if (length(allnohole@polygons[[1]]@Polygons) > 1)
  {
    new_sp1 = SpatialPolygons( list( Polygons( list( new_poly@polygons[[1]]@Polygons[[1]] ),
                                               ID = k)))
    area1 <- area(new_sp1)
    for (k in 2:length(allnohole@polygons[[1]]@Polygons))
    {
      new_sp2 = SpatialPolygons( list( Polygons( list( new_poly@polygons[[1]]@Polygons[[k]] ),
                                                 ID = k)))
      area2 <- area(new_sp2)
      
      if (area2 > area1)
      {
        area1 <- area2
        new_sp1 <- new_sp2
      }
    }
    allnohole <- new_sp1
  }
  
  # add dataframe
  allnoholedata = SpatialPolygonsDataFrame(allnohole, data = new_poly@data, match.ID=FALSE)
  
  # join polygons
  if (j==1)
  { bound_unique = allnoholedata }
  else
  { bound_unique = rbind(bound_unique,allnoholedata) }
  
  # add height
  z90 <- c(z90, quantile(cloudpoint$z[which(cloudpoint$cluster == new_poly@data$layer)],0.9))
  zmax <- c(zmax, max(cloudpoint$z[which(cloudpoint$cluster == new_poly@data$layer)]))
}

bound_unique@data$area_auto = sapply(slot(bound_unique, "polygons"), slot, "area")
bound_unique@data$z90 = z90
bound_unique@data$zmax = zmax

# minimum size for a polygon to be included in the SpatialPolygonsDataFrame
bound_unique2 = bound_unique[bound_unique$area_auto > 4.5,]

# projection
proj4string(bound_unique2) <- CRS("+proj=utm +zone=22 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# save
writeOGR(obj = bound_unique2 ,dsn = "directory",
         layer = "crowns_shapes",driver='ESRI Shapefile')


