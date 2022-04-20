library(dplyr)
library(stars)
library(sf)
library(sp)
library(tidyverse)
library(akima)
library(tidyr)

options(scipen=999)


raster_to_points = function(raster_obj){
  #  get centers of raster cells
  rast_pnts = st_as_sfc(raster_obj, as_points = TRUE)
  rast_values = st_extract(raster_obj, rast_pnts)
  rast_pnt_grid = st_coordinates(rast_values) %>%
    as.data.frame() %>%
    cbind(as.data.frame(rast_values)[1])
  colnames(rast_pnt_grid) = c('lon_p', 'lat_p', 'value')
  return (rast_pnt_grid)
}

raster_to_bbox = function(raster_obj){
  box = st_bbox(raster_obj)
  envelope = box[c(1,3,2,4)]
  # polygon coords
  box_coords = data.frame(lon_p = c(as.numeric(envelope[1]),
                                    as.numeric(envelope[2]),
                                    as.numeric(envelope[2]),
                                    as.numeric(envelope[1]),
                                    as.numeric(envelope[1])),
                          lat_p = c(as.numeric(envelope[3]),
                                    as.numeric(envelope[3]),
                                    as.numeric(envelope[4]),
                                    as.numeric(envelope[4]),
                                    as.numeric(envelope[3])))
  return (box_coords)
}

linear_interpolation = function(rast_tab, hex_tab){
  rast_tab = as.data.frame(rast_tab)
  rast_tab[is.na(rast_tab)] = 0
  interpolated = interpp(x = rast_tab$lon_p,
                         y = rast_tab$lat_p,
                         z = rast_tab$value,
                         xo = hex_tab$lon_hex,
                         yo = hex_tab$lat_hex,
                         linear = TRUE,
                         duplicate = "mean") %>%
    as.data.frame()
  #   надо предусмотреть сохранение целых значений, если растр с кач. хар-ками
}

bilinear_interpolation = function(raster_obj, hex_tab){
  hex_pnts = st_as_sf(hex_tab,
                      coords = c("lon_hex", "lat_hex"),
                      crs = 4326)
  hex_values = st_extract(raster_obj,
                          hex_pnts,
                          bilinear = TRUE) %>%
    st_drop_geometry() %>%
    mutate(lon_hex = hex_tab$lon_hex) %>%
    mutate(lat_hex = hex_tab$lat_hex)
  colnames(hex_values) = c('z', 'x', 'y')
  return(hex_values)
}

raster_to_hex = function(rast_obj, h3_level, intrp_method = 'bilinear'){
  rast_pnts = raster_to_points(rast_obj)
  rast_extent = raster_to_bbox(rast_obj)
  hex_centers = h3:::hex_centers_inbbox(rast_extent$lon_p,
                                        rast_extent$lat_p,
                                        h3_level) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  hex_centers$ind = rownames(hex_centers)
  colnames(hex_centers) = c('lon_hex', 'lat_hex', 'ind')
  rownames(hex_centers) = seq(1, length(hex_centers$ind), 1)
  if (intrp_method == 'linear'){
    new = linear_interpolation(rast_pnts,  hex_centers)
  } else if (intrp_method == 'bilinear'){
    new = bilinear_interpolation(rast_obj,  hex_centers)
  }
  new$h3_ind = hex_centers$ind
  return(new)
}
