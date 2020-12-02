# 02/12/2020

# Individual tree crown segmentation method based on ALS + high resolution RGB
# Refinement of the segmentation realized by AMS3D (Computree platform by ONF) with DeepForest (python package by Ben Weinstein)

# AMS3D is a segmentation method using the 3d lidar point cloud
# DeepForest draws bounding boxes around crowns on a RGB image

# here, results of both methods are combined to correct over-segmentation of the lidar point cloud (frequent problem with AMS3D) using RGB images.

# packages
library(rgdal)
library(lidR)


#### load data

# AMS3D results (crowns_shapes can be created by the code ShapesCreationFromAMS3D.R from the raw results of AMS3D)
shpAMS <- readOGR("crowns_shapes")
# create an ID for each AMS3D segment
shpAMS@data$IDAMS <- 1:length(shpAMS)
  
# boxes DeepForest
shpDeep <- readOGR("TrainedModel")
# compute area of each box
shpDeep@data$area_auto = sapply(slot(shpDeep, "polygons"), slot, "area")
# create an ID for each DeepForest segment
shpDeep@data$IDDeep <- 1:length(shpDeep)


##### merge AMS3D segments if they intersect with the same DeepForest box

## first, for each AMS3D segment intersecting by more than half of its area with a DeepForest bounding box, set idDeep


# idDeep of AMS3D segments set to NA
shpAMS@data$idDeep <- NA
  
for (i in 1:length(shpDeep))
{
  # for each bounding box of DeepForest
  bb <- shpDeep[i,]
  # get coordinates
  coords <- bb@polygons[[1]]@Polygons[[1]]@coords
  # change square box to octagon
  oct <- rbind(c(coords[1,1], coords[2,2]+ (coords[1,2] - coords[2,2])* 2/3),
               c(coords[1,1], coords[2,2]+ (coords[1,2] - coords[2,2])* 1/3),
               c(coords[3,1]+ (coords[2,1] - coords[3,1])* 2/3, coords[2,2]),
               c(coords[3,1]+ (coords[2,1] - coords[3,1])* 1/3, coords[2,2]),
               c(coords[3,1], coords[4,2]+ (coords[3,2] - coords[4,2])* 2/3),
               c(coords[3,1], coords[4,2]+ (coords[3,2] - coords[4,2])* 1/3),
               c(coords[4,1]+ (coords[5,1] - coords[4,1])* 1/3, coords[5,2]),
               c(coords[4,1]+ (coords[5,1] - coords[4,1])* 2/3, coords[5,2]),
               c(coords[1,1], coords[2,2]+ (coords[1,2] - coords[2,2])* 2/3))
  # create SpatialPolygons from the octoagon
  oct <- Polygons(list(Polygon(oct)), ID='A')
  oct <- SpatialPolygons(list(oct))
  
  # get the intersection with AMS3D segments
  inter <- raster::intersect(oct, shpAMS)
  
  # if some segments intersect:
  if (!is.null(inter))
  {
    # get the area of the intersections
    area_inter = sapply(slot(inter, "polygons"), slot, "area")
    # get ID for each segment intersecting
    idAMS <- inter@data$IDAMS
    # get AMS3D intersecting segments
    segAMS <- shpAMS[shpAMS$IDAMS %in% idAMS,]
    # compute % of segment included in the octagon
    ratio <- area_inter / segAMS@data$area_auto
    # get segments that have more than half of their area in the octagon
    ID_inside <- idAMS[ratio > 0.5]
    # change idDeep for these segments to i
    shpAMS@data$idDeep[ID_inside] <- i
  }
}
  
# optional: save before merging
writeOGR(shpAMS, "AMS_with_id_deep", driver='ESRI Shapefile')


## second, merge AMS3D segments if they have the same idDeep

# the merged polygons will be stored in "merged"
merged <- shpAMS

# get idDeep that have been attributed to AMS3D segments
idD <- unique(merged@data$idDeep)
idD <- idD[!is.na(idD)]

# for each idDeep:
for (i in idD)
{
  # create a new polygon with the corresponding segments 
  newPoly <- aggregate(merged[which(merged@data$idDeep == i),],dissolve=T)
  # create data frame for this new polygon
  dat <- data.frame(layer = 1, area_auto = area(newPoly), IDAMS = NA, idDeep = i)
  # create SpatialPolygonsDataFrame for this new polygon
  newPoly <- SpatialPolygonsDataFrame(newPoly, dat, match.ID = F)
  
  # remove oversegmented polygons
  merged <- merged[-which(merged@data$idDeep == i),]
  # add merged polygon
  merged <- rbind(merged, newPoly)
}

# write resulting SpatialPolygonsDataFrame
writeOGR(merged, "AMS_corrected_DF", driver='ESRI Shapefile' )
  
