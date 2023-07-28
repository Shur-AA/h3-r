library(dplyr)
library(stars)
library(sf)
library(sp)
library(tidyverse)
library(akima)
library(tidyr)
library(readxl)

options(scipen=999)
rpath = "D:/3_Проекты/РФФИ сток/data/etopo15/africa_orange.tif"
rast = read_stars(rpath) # input 1-band raster file



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


initial_resolution = get_rast_cellarea(rast)
start_h3_l = choose_h3_level(initial_resolution)


tab = h3::h3_raster_to_hex(rast, start_h3_l)
c2 = "88ada22a35fffff"
k = h3:::flow_dir(tab$h3_ind, tab$z, c2)

t = readxl::read_excel("C:/Users/user/Downloads/fd.xlsx")
k = h3::h3_flow_acc(t$from, t$to)
write.csv(k, 'C:/Users/user/Downloads/fd.csv')




tab[is.na(tab)] = -9999
tab = tab1 %>% filter(!is.na(tab$z))

center = "88adaec893fffff"



local_raster_sum = h3::h3_simple_sum(tab1$h3_ind,
                                     tab1$z,
                                     tab2$h3_ind,
                                     tab2$z) %>%
  as.data.frame()

global_raster_stat = h3::h3_global_extremum(tab2$h3_ind,
                                            tab2$z,
                                            'max')
stat_type = 'max'
zonal_raster_sum = h3::h3_zonal_statistics(tab1$h3_ind,
                                           as.integer(tab1$z),
                                           tab2$h3_ind,
                                           tab2$z,
                                           stat_type,
                                           F) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(ind = rownames(.))
rownames(zonal_raster_sum) = seq(1, length(zonal_raster_sum$ind), 1)
colnames(zonal_raster_sum) = c('zone_code', stat_type, 'hex_ind')

focal_raster = h3::h3_simple_focal(tab1$h3_ind,
                                   tab1$z,
                                   stat_type, 2) %>%
  list() %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = 'indexes',
               values_to = 'values')
coords_focal = h3::h3_indexes_to_coords(focal_raster$indexes)%>%
  do.call(rbind, .) %>%
  as.data.frame()
coords_focal$ind = rownames(coords_focal)
colnames(coords_focal) = c('lon_hex', 'lat_hex', 'ind')
rownames(coords_focal) = seq(1, length(coords_focal$ind), 1)
coords_focal = cbind(coords_focal, focal_raster$values)
write.csv(k, "C:/Users/a.shurygina/Downloads/testfoc2_h3.csv")


hi = c('81263ffffffffff', '81267ffffffffff', '8126bffffffffff', '8126fffffffffff',
       '81273ffffffffff', '81277ffffffffff', '8127bffffffffff')
g = lapply(hi, function(X){h3_cell_azimuth(X)})


