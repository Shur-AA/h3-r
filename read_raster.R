library(dplyr)
library(stars)
library(sf)
library(sp)
library(tidyverse)
library(akima)
library(tidyr)

options(scipen=999)
rpath1 = "./test_samples/ufa_igbp.tif"
rpath2 = "./test_samples/ufa_glc.tif"
rast1 = read_stars(rpath1) # input 1-band raster file 1
rast2 = read_stars(rpath2) # input 1-band raster file 2



tab1 = h3::h3_raster_to_hex(rast1, 7)
tab2 = h3::h3_raster_to_hex(rast2, 7)

local_raster_sum = h3::h3_simple_sum(tab1$h3_ind,
                                      tab1$z,
                                      tab2$h3_ind,
                                      tab2$z) %>%
                                      as.data.frame()

global_raster_stat = h3::h3_global_extremum(tab2$h3_ind,
                                            tab2$z,
                                            'max')
stat_type = 'minority'
zonal_raster_sum = h3::h3_zonal_statistics(tab1$h3_ind,
                                          as.integer(tab1$z),
                                          tab2$h3_ind,
                                          tab2$z,
                                          stat_type) %>%
                                          do.call(rbind, .) %>%
                                          as.data.frame() %>%
                                          mutate(ind = rownames(.))
rownames(zonal_raster_sum) = seq(1, length(zonal_raster_sum$ind), 1)
colnames(zonal_raster_sum) = c('zone_code', stat_type, 'hex_ind')


