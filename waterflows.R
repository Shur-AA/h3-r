library(dplyr)
library(stars)
library(sf)
library(tidyverse)
library(akima)
library(tidyr)
library(readxl)

options(scipen=999)
rpath = "D:/3_Проекты/РФФИ сток/data/etopo15/africa_etopo60.tif"
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


initial_resolution = get_rast_cellarea(rast)
start_h3_l = choose_h3_level(initial_resolution)

# tiling raster
rast_extent = raster_to_bbox(rast)
hex_bnd_cntr = h3:::hex_boundary_inbbox(rast_extent$lon_p,
                                       rast_extent$lat_p,
                                       2, start_h3_l) %>%
  do.call(rbind, .) %>%
  as.data.frame()
hex_bnd_cntr$ind = rownames(hex_bnd_cntr)
rownames(hex_bnd_cntr) = seq(1, length(hex_bnd_cntr$ind), 1)


for (h in 728:nrow(hex_bnd_cntr)){
  print(h)
  # преобразуем координаты и делаем полигон
  ahex = hex_bnd_cntr[h,] %>% select(-ind)
  plg_m = matrix(c(ahex), ncol = 2, byrow = TRUE)
  ends = which(apply(plg_m, 1, function(x) identical(x, plg_m[1,])))
  a = max(ends) + 1
  if (a != 2){
    b = nrow(plg_m)
    plg_m = plg_m[-(a:b),] %>% as.data.frame()
  }else{
    plg_m = rbind(plg_m, plg_m[1,]) %>% as.data.frame()
  }
  plg_m$V1 = as.double(plg_m$V1)
  plg_m$V2 = as.double(plg_m$V2)
  plg_m = as.matrix(plg_m)
  coords = list(plg_m)
  pol = st_polygon(coords) %>% st_sfc()
  st_crs(pol) = st_crs(4326)

  # делаем буфер на 5,3 км, чтобы тайлы были внахлёст
  pol = st_transform(pol, crs = 3857) %>%
    st_buffer(dist = 6000) %>%
    st_transform(crs = 4326)


  tile = st_crop(rast, pol)


  tab = h3::h3_raster_to_hex(tile, start_h3_l)
  fd = h3::h3_flow_dir(tab$h3_ind, tab$z, hex_bnd_cntr$ind[h])

  if (length(fd) > 2){
    write.csv(fd, paste('C:/Users/user/Downloads/gidro/', 'fd', h, '.csv', sep = ''))
    fd = read.csv(paste('C:/Users/user/Downloads/gidro/', 'fd', h, '.csv', sep = ''))
    colnames(fd) = c('from', 'to')
    fa = h3::h3_flow_acc(fd$from, fd$to)
    write.csv(fa, paste('C:/Users/user/Downloads/gidro/', 'fa', h, '.csv', sep = ''))
  }
}





#tab = h3::h3_raster_to_hex(rast, start_h3_l)
#c2 = "88ada22a35fffff"
#c2 = "87ada22a0ffffff"
#k = h3:::flow_dir(tab$h3_ind, tab$z, c2)

#t = readxl::read_excel("C:/Users/user/Downloads/fd10.xlsx")
#k = h3::h3_flow_acc(t$from, t$to)
#write.csv(k, 'C:/Users/user/Downloads/fa10.csv')


#tab[is.na(tab)] = -9999
#tab = tab1 %>% filter(!is.na(tab$z))
