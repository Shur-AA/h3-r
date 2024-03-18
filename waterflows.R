library(dplyr)
library(stars)
library(sf)
library(tidyverse)
library(akima)
library(tidyr)
library(readxl)
library(RPostgres)


options(scipen=999)
rpath = "D:/3_Проекты/РФФИ сток/data/etopo15/africa_orange_problem2.tif"
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
#start_h3_l = 7


# tiling raster
rast_extent = raster_to_bbox(rast)
hex_bnd_cntr = h3:::hex_boundary_inbbox(rast_extent$lon_p,
                                       rast_extent$lat_p,
                                       3, start_h3_l) %>%
  do.call(rbind, .) %>%
  as.data.frame()
hex_bnd_cntr$ind = rownames(hex_bnd_cntr)
rownames(hex_bnd_cntr) = seq(1, length(hex_bnd_cntr$ind), 1)


for (h in c(5581)){
  print(h)
  h = 933
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

  # делаем буфер на 6 км, чтобы тайлы были внахлёст
  pol_plus = st_transform(pol, crs = 3857) %>%
              st_buffer(dist = 50000) %>%
              st_transform(crs = 4326)

  # получаем широкий растр
  tile = st_crop(rast, pol_plus)
  # конвертируем его в сетку
  tab = h3::h3_raster_to_hex(tile, start_h3_l)
  # выбираем исходным шестиугольником только те ячейки, чьи узлы в него попадают
  tab_pnt = st_as_sf(tab, coords = c('x', 'y'), crs = 4326)
  tab1 = tab_pnt[pol,]

  #  tab1 - стандартный шестиугольник; tab - расширенный шестиугольник


  fd_ext = h3:::fd_experiment(tab$h3_ind, tab$z)
  write.csv(fd_ext, paste('C:/Users/user/Downloads/gidro/', 'fd', h, '.csv', sep = ''))
  fd_ext = read.csv(paste('C:/Users/user/Downloads/gidro/', 'fd', h, '.csv', sep = ''))
  colnames(fd_ext) = c('from', 'to')
  fa = h3::h3_flow_acc(fd_ext$from, fd_ext$to)
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/', 'fa', h, '.csv', sep = ''))
  fa = read.csv(paste('C:/Users/user/Downloads/gidro/', 'fa', h, '.csv', sep = '')) %>%
    filter(`x` > 80)
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/', 'fa', h, '.csv', sep = ''))
}


# Эксперимент с расширяющейся областью (от ячейки 1010 l3)
start_tile = 1010
hex_dem = h3::h3_raster_to_hex(rast, 8)
hex_dem = select(hex_dem, -`x`, -`y`) %>%
          na.omit()
for (r in 101:150){
  this_tile = h3:::cell_vecinity_circle(hex_bnd_cntr[start_tile,]$ind, r) %>%
              data.frame() %>%
              left_join(hex_dem, by = c("." = "h3_ind"))
  colnames(this_tile) = c('h3_ind', 'z')


  fd_ext = h3:::fd_experiment(this_tile$h3_ind, this_tile$z)
  write.csv(fd_ext, paste('C:/Users/user/Downloads/gidro/expanding/', 'fd', r, '.csv', sep = ''))
  fd_ext = read.csv(paste('C:/Users/user/Downloads/gidro/expanding/', 'fd', r, '.csv', sep = ''))
  colnames(fd_ext) = c('from', 'to')
  fa = h3::h3_flow_acc(fd_ext$from, fd_ext$to)
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/expanding/', 'fa', r, '.csv', sep = ''))
  fa = read.csv(paste('C:/Users/user/Downloads/gidro/expanding/', 'fa', r, '.csv', sep = ''))
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/expanding/', 'fa', r, '.csv', sep = ''))
}


# Эксперимент со скользящим окном
hex_dem = h3::h3_raster_to_hex(rast, 8)
hex_dem = select(hex_dem, -`x`, -`y`) %>%
  na.omit()
for (c in hex_dem$h3_ind){
  this_tile = h3:::cell_vecinity_circle(c, 10) %>%
    data.frame() %>%
    left_join(hex_dem, by = c("." = "h3_ind"))
  colnames(this_tile) = c('h3_ind', 'z')

  fd_ext = h3:::fd_experiment(this_tile$h3_ind, this_tile$z)
  write.csv(fd_ext, paste('C:/Users/user/Downloads/gidro/walk/', 'fd_', c, '.csv', sep = ''))
  fd_ext = read.csv(paste('C:/Users/user/Downloads/gidro/walk/', 'fd_', c, '.csv', sep = ''))
  colnames(fd_ext) = c('from', 'to')
  fa = h3::h3_flow_acc(fd_ext$from, fd_ext$to)
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/walk/', 'fa_', c, '.csv', sep = ''))
  fa = read.csv(paste('C:/Users/user/Downloads/gidro/walk/', 'fa_', c, '.csv', sep = ''))
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/walk/', 'fa_', c, '.csv', sep = ''))
}




library(DBI)
library(deckgl)
#create connection object
con <- dbConnect(RPostgres::Postgres(),
                 user="postgres",
                 password="17041996",
                 host="localhost",
                 port=5432,
                 dbname="gidro")
dbListTables(con)   #list all the tables
cor = dbReadTable(con, "africa_fd")
#dbWriteTable(con, name, new_tab, row.names=F, overwrite=F (append=T))



#fatab = read.csv('C:/Users/user/Downloads/gidro/AfricaAcc.csv',
 #                skip = 30000000, header = F, nrows = 80000000)
fatab = select(tab, h3_ind, z)


for (num in 1:10){
  num = 150
  fa = read.csv(paste('C:/Users/user/Downloads/gidro/expanding/', 'fa', num, '.csv', sep = '')) %>%
      select(`X`, `x`)
  colnames(fa) = c('h3_ind', 'z')

  properties <- list(
    getHexagon = ~h3_ind,
    getFillColor = JS("d => [0, (1 - d.z / 100) * 255, 250]"),
    #getElevation = ~acc,
    #elevationScale = 2,
    getTooltip = "{{h3_ind}}: {{z}}"
  )

  deck <- deckgl(zoom = 7,
                 latitude=-29.7899884,
                 longitude=23.4763848) %>%
    add_h3_hexagon_layer(data = hex_dem, properties = properties) %>%
    add_control("Expand") %>%
    add_basemap()

  if (interactive()) deck
}


