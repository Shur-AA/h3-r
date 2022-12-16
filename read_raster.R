library(dplyr)
library(stars)
library(sf)
library(sp)
library(tidyverse)
library(akima)
library(tidyr)

options(scipen=999)
rpath1 = "./test_samples/ufa_igbp.tif"
#rpath1 = "./test_samples/etopo30.tiff"
rpath2 = "./test_samples/ufa_glc.tif"
rast1 = read_stars(rpath1) # input 1-band raster file 1
rast2 = read_stars(rpath2) # input 1-band raster file 2




tab1 = h3::h3_raster_to_hex(rast1, 8)
tab2 = h3::h3_raster_to_hex(rast2, 8)


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
#write.csv(coords_focal, "C:/Users/a.shurygina/Downloads/testfoc2_h3.csv")


hi = c('81263ffffffffff', '81267ffffffffff', '8126bffffffffff', '8126fffffffffff',
       '81273ffffffffff', '81277ffffffffff', '8127bffffffffff')
g = lapply(hi, function(X){h3_cell_azimuth(X)})

