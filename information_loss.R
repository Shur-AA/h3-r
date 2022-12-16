library(dplyr)
library(stars)
library(sf)
library(sp)
library(tidyverse)
library(akima)
library(tidyr)


options(scipen=999)

determine_epsg_code = function(longitude){
  return(paste(c('326', as.character(longitude %/% 6 + 31)), collapse=""))
}


# calculate raster cell area
get_rast_cellarea = function(raster_obj){
  # mean longitude of the image
  mean_lon = ((st_bbox(raster_obj)[1] + st_bbox(raster_obj)[3]) / 2)[[1]]
  metr_crs = st_crs(as.numeric(determine_epsg_code(mean_lon)))
  new_rast = st_transform(raster_obj, metr_crs)
  # dimension in meters on x and y
  x_dim = ((st_bbox(new_rast)[3] - st_bbox(new_rast)[1]) / dim(new_rast)['x'])[[1]]
  y_dim = ((st_bbox(new_rast)[4] - st_bbox(new_rast)[2]) / dim(new_rast)['y'])[[1]]
  return(x_dim * y_dim)
}


choose_h3_level = function(raster_cell_area){
  v = c(0:15, 4357449416078.392, 609788441794.134,
        86801780398.997, 12393434655.088, 1770347654.491,
        252903858.182, 36129062.164, 5161293.360, 737327.598,
        105332.513, 15047.502, 2149.643, 307.092,
        43.870, 6.267, 0.895)
  h3_cells_info = matrix(v, ncol = 2)
  closeness = abs(h3_cells_info[, 2] - raster_cell_area / 2)
  idx = match(min(closeness), closeness)
  return(h3_cells_info[idx, 1])
}



rpath1 = "./test_samples/ufa_igbp.tif"
#rpath1 = "./test_samples/etopo30.tiff"
#rpath2 = "./test_samples/ufa_glc.tif"
rast1 = read_stars(rpath1) # input 1-band raster file 1
#rast2 = read_stars(rpath2) # input 1-band raster file 2


tab1 = h3_raster_to_hex(rast1, choose_h3_level(get_rast_cellarea(rast1)))
tab2 = h3:::resample_up("majority", tab1$h3_ind,
                              as.integer(tab1$z))%>%
                                               as.data.frame()
tab22 = h3_raster_to_hex(rast1, choose_h3_level(get_rast_cellarea(rast1))-1)









