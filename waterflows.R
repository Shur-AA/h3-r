library(dplyr)
library(stars)
library(sf)
library(tidyverse)
library(akima)
library(tidyr)
library(readxl)
library(RPostgres)
library(Matrix)
library(data.table)


options(scipen=999)
rpath = "D:/3_Проекты/РФФИ сток/data/etopo15/africa60_3.tif"
#rpath = "D:/3_Проекты/РФФИ сток/data/congo1.tif"
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
  #new_rast = raster_obj
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


# read common polygons

bas4 = st_read("D:/3_Проекты/РФФИ сток/data/todelete.geojson")
#%>%filter(HYBAS_ID == 1041213630)

bas_plus = st_transform(bas4, crs = 3857) %>%
            st_buffer(dist = 500) %>%
            st_transform(crs = 4326)

bastile = st_crop(rast, bas_plus)

tab = h3::h3_raster_to_hex(bastile, 8) %>% na.omit()






for (h in c(12, 13, 16)){
  print(h)
  h = 999
  start.time <- Sys.time()
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

  # делаем буфер на 5 км, чтобы тайлы были внахлёст
  pol_plus = st_transform(pol, crs = 3857) %>%
              st_buffer(dist = 200) %>%
              st_transform(crs = 4326)



  # получаем широкий растр
  tile = st_crop(rast, pol_plus)


  #rast = st_transform(rast, 4326)
  #tab = h3::h3_raster_to_hex(rast, 8) %>% na.omit()
  # конвертируем его в сетку
  tab = h3::h3_raster_to_hex(tile, start_h3_l) %>% na.omit()
  #
  # выбираем исходным шестиугольником только те ячейки, чьи узлы в него попадают
  tab_pnt = st_as_sf(tab, coords = c('x', 'y'), crs = 4326)
  tab1 = tab_pnt[pol,]

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken

  #  tab1 - стандартный шестиугольник; tab - расширенный шестиугольник

  fdem = h3:::fill_depr_jd(tab$h3_ind, tab$z)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken


  write.csv(fdem, paste('C:/Users/user/Downloads/gidro/', 'fd', h, '.csv', sep = ''))
  fdem = read.csv(paste('C:/Users/user/Downloads/gidro/', 'fd', h, '.csv', sep = ''))
  colnames(fdem) = c('from', 'to')
  fdem = filter(fdem, to != 'edge')
  write.csv(fdem, paste('C:/Users/user/Downloads/gidro/', 'fd', h, '.csv', sep = ''))
  fa = h3::h3_flow_acc(fdem$from, fdem$to)
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/', 'fa', h, '.csv', sep = ''))
  fa = read.csv(paste('C:/Users/user/Downloads/gidro/', 'fa', h, '.csv', sep = '')) %>%
    filter(`x` > 40)
  write.csv(fa, paste('C:/Users/user/Downloads/gidro/', 'fa', h, '.csv', sep = ''))

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken

}



fdem = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_nsa7.csv')
fdem = select(fdem, `X`, `x`)
colnames(fdem) = c('from', 'to')
fa = h3::h3_flow_acc(fdem$from, fdem$to)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_nsa7fa.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_nsa7fa.csv') %>%
  filter(`x` > 20)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_nsa7fa.csv')



fdem = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_fd.csv')
fdem = filter(fdem, to != 'edge')
fa = h3::h3_flow_acc(fdem$from, fdem$to)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_fa.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_fa.csv') %>%
  filter(`x` > 50)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_faf.csv')


fdem = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fd.csv')
fdem = filter(fdem, to != 'edge')
fa = h3::h3_flow_acc(fdem$from, fdem$to)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fa.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fa.csv') %>%
  filter(`x` > 50)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_faf.csv')



fdem = read.csv('C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_fd.csv')
fdem = filter(fdem, to != 'edge')
fa = h3::h3_flow_acc(fdem$from, fdem$to)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_fa.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_fa.csv') %>%
  filter(`x` > 50)
write.csv(fa, 'C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_faf.csv')


# negative buffering
buf_size = -800 * 7

nb = read.csv('C:/Users/user/Downloads/gidro/Congo/1041213640/fd333.csv')
coords = h3::h3_indexes_to_coords(nb$from) %>%
          do.call(rbind, .) %>%
          as.data.frame()
ch = st_as_sf(coords, coords = c('V1', 'V2'), crs = 4326) %>%
        st_union() %>%
        st_convex_hull() %>% st_transform(pol, crs = 3857) %>%
        st_buffer(dist = buf_size) %>%
        st_transform(crs = 4326)

hex_centers = h3:::hex_centers_inbbox(ch[[1]][[1]][,1],
                                      ch[[1]][[1]][,2],
                                      8) %>%
              do.call(rbind, .) %>%
              as.data.frame()
hex_centers$ind = rownames(hex_centers)
shortened = select(hex_centers, ind)

buf_res = inner_join(nb, shortened, by=c("from" = "ind"))

write.csv(buf_res, 'C:/Users/user/Downloads/gidro/Congo/1041213640/fd333b.csv')






# COTAT testing
library(dplyr)
library(stars)
library(sf)
library(tidyverse)
library(akima)
library(tidyr)
library(readxl)
library(RPostgres)
library(Matrix)
library(data.table)


fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 6, 250)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_cotat6.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_cotat3_90000.csv')
colnames(cfd) = c('from', 'to')
#vg = h3:::generalisation_verification(fd$from, fd$to, cfd$from, cfd$to, 10000)

fd = read.csv('C:/Users/user/Downloads/gidro/fd17.csv')
fd = select(fd, from, to)
dt = h3:::dren_tree(fd$from, fd$to)
write.csv(dt, 'C:/Users/user/Downloads/gidro/drent_17.csv')

fd = read.csv('C:/Users/user/Downloads/gidro/fd17.csv')
fd = select(fd, from, to)
cfd = h3:::vvrfra(fd$from, fd$to, 4)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/vvrfra_.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/vvrfra_.csv')
colnames(cfd) = c('from', 'to')

fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_fa.csv')
fa = select(fa, `X`, `x`)
colnames(fa) = c('ind', 'val')
cfd = h3:::nsa(fa$ind, fa$val, 6)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041213640/congo_1041213640_nsa6.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_nsa6.csv')
colnames(cfd) = c('from', 'to')


fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_fd.csv')
fd = select(fd, from, to)
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_cotat3_90000.csv')
cfd = select(cfd, `X`, `x`)
colnames(cfd) = c('from', 'to')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041174950/congo_1041174950_fa.csv')
fa = select(fa, `X`, `x`)
colnames(fa) = c('ind', 'val')
vg = h3:::generalisation_verification(fd$from, fd$to, cfd$from, cfd$to, fa$ind, fa$val, 10)

# DELETE

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156960/congo_1041156960_cotat4.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 7, 40)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041156960/congo_1041156960_cotat7.csv')

fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156960/congo_1041156960_fa.csv')
fa = select(fa, `X`, `x`)
colnames(fa) = c('ind', 'val')
cfd = h3:::nsa(fa$ind, fa$val, 7)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_nsa7.csv')


fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 6, 250)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_cotat6.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fa.csv')
fa = select(fa, `X`, `x`)
colnames(fa) = c('ind', 'val')
cfd = h3:::nsa(fa$ind, fa$val, 6)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_nsa6.csv')

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 5, 1700)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_cotat5.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fa.csv')
fa = select(fa, `X`, `x`)
colnames(fa) = c('ind', 'val')
cfd = h3:::nsa(fa$ind, fa$val, 5)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_nsa5.csv')

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 4, 12000)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_cotat4.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fa.csv')
fa = select(fa, `X`, `x`)
colnames(fa) = c('ind', 'val')
cfd = h3:::nsa(fa$ind, fa$val, 4)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_nsa4.csv')


fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156960/congo_1041156960_cotat4p.csv')
fd = select(fd, `X`, `x`)
colnames(fd) = c('from', 'to')
cfd = h3:::cotat(fd$from, fd$to, 3, 40)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041156960/congo_1041156960_cotat3p.csv')

fa = read.csv('C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_fa.csv')
fa = select(fa, `X`, `x`)
colnames(fa) = c('ind', 'val')
cfd = h3:::nsa(fa$ind, fa$val, 3)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041259940/congo_1041259940_nsa3.csv')








# cotat for night

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 7, 40)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat7_40.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat7_40.csv')
colnames(cfd) = c('from', 'to')
a = h3:::generalisation_verification(fd$from, fd$to, cfd$from, cfd$to, 10000)

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 6, 250)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat6_250.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat6_250.csv')
colnames(cfd) = c('from', 'to')
b = h3:::generalisation_verification(fd$from, fd$to, cfd$from, cfd$to, 10000)

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 5, 1700)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat5_1700.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat5_1700.csv')
colnames(cfd) = c('from', 'to')
c = h3:::generalisation_verification(fd$from, fd$to, cfd$from, cfd$to, 10000)

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 4, 12000)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat4_12000.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat4_12000.csv')
colnames(cfd) = c('from', 'to')
d = h3:::generalisation_verification(fd$from, fd$to, cfd$from, cfd$to, 10000)

fd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_fd.csv')
fd = select(fd, from, to)
cfd = h3:::cotat(fd$from, fd$to, 3, 90000)
write.csv(cfd, 'C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat3_90000.csv')
cfd = read.csv('C:/Users/user/Downloads/gidro/Congo/1041156950/congo_1041156950_cotat3_90000.csv')
colnames(cfd) = c('from', 'to')
e = h3:::generalisation_verification(fd$from, fd$to, cfd$from, cfd$to, 10000)











fdem = read.csv('C:/Users/user/Downloads/gidro/cotat11.csv')
colnames(fdem) = c('from', 'to')
fdem = filter(fdem, to != 'edge')
fa = h3::h3_flow_acc(fdem$from, fdem$to)
write.csv(fa,'C:/Users/user/Downloads/gidro/fa_11.csv')
fa = read.csv('C:/Users/user/Downloads/gidro/fa_11.csv') %>%
  filter(`x` > 1)
write.csv(fa,'C:/Users/user/Downloads/gidro/fa_11.csv')



write.csv(fdem, 'C:/Users/user/Downloads/gidro/vvv.csv')
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





// labeling watersheds (Step 3)

int w = 0; // number (label) of a watershed
std::unordered_map <std::string, int> watersheds;
std::vector<std::string> poor_points;
// go through all 'from' cells in current fd table
for (auto const & dir_pair : fdtab){

  int label = w;

  std::cout<<label<<std::endl;
  // has current cell exit of type 1 (clear case)?
    // if no exit
  if (dir_pair.second == "edge" ||
      dir_pair.second == "zero_drops" ||
      dir_pair.second == "undef"){
    // check if it is in poor points list
    if (std::find(poor_points.begin(), poor_points.end(),
                  dir_pair.first) != poor_points.end()){
      // if it is already in poor points list
      // take corresponding watershed label
      label = watersheds[dir_pair.first];
    }else{
      // if it is not in poor points list yet
      poor_points.push_back(dir_pair.first);
      watersheds[dir_pair.first] = ++w;
    }
  }else{
    std::vector<std::string> this_path; // path to poor point
    this_path.push_back(dir_pair.first);
    std::string toind = dir_pair.second; // where to flow
    std::map <std::string, int> control_loops; // table to count visits of each cell of a stream
    while(true){
      if (fdtab.find(toind) != fdtab.end()){
        this_path.push_back(toind);
        // checking loops
        control_loops[toind]++;
        bool flag = false;
        for (auto const & loop_pair : control_loops){
          // если посетили больше 1 раза, значит уже зациклились
          if (loop_pair.second > 1){flag = true;}
        }
        if (flag){
          std::cout<<"There was a loop!"<<std::endl;
          break;
        }else{
          toind = fdtab[toind];
        }
      }else{
        if (dir_pair.second != "edge" ||
            dir_pair.second != "zero_drops" ||
            dir_pair.second != "undef"){
          this_path.push_back(toind);
        }
        break;
      }
    }
    // now we have list with path ending in some poor point
    // check if the end is in the poor points list
    int pp_len = this_path.size();
    if (pp_len > 0){
      if (std::find(poor_points.begin(), poor_points.end(),
                    this_path[pp_len - 1]) != poor_points.end()){
        // if it is already in poor points list
        // take corresponding watershed label
        label = watersheds[this_path[pp_len - 1]];
      }else{
        // if it is not in poor points list yet
        poor_points.push_back(this_path[pp_len - 1]);
        watersheds[this_path[pp_len - 1]] = ++w;
        label = w;
      }
      for (int i; i < (pp_len - 1); i++){
        watersheds[poor_points[i]] = label;
      }
    }
  }
}





// check if it is in poor points list
if (std::find(poor_points.begin(), poor_points.end(),
              dir_pair.first) != poor_points.end()){
  // if it is already in poor points list
  // take corresponding watershed label
  label = watersheds[dir_pair.first];
}else{
  // if it is not in poor points list yet
  poor_points.push_back(dir_pair.first);
  watersheds[dir_pair.first] = ++w;
}



// pre-post processing for flat areas which are next to lpp
// and may not outpour via them
for (auto const & acell : w_lowest_pp){
  // go through all lpp cells
  // get the cell's neighbors
    std::vector<std::string> this_vicinity = cell_vecinity(acell.second, 1);
    for (auto const & vic_ind : this_vicinity){
      // if the cell is from flat (zd) area and has same height with current cell
      if (std::find(need_postprocess1.begin(), need_postprocess1.end(), vic_ind)
            != need_postprocess1.end() && geotab[vic_ind] == geotab[acell.second]){
        fdtab1[vic_ind] = acell.second;
        np_copy.erase(std::remove(np_copy.begin(),
                                  np_copy.end(), vic_ind),
                                  np_copy.end());
        need_postprocess1.erase(std::remove(need_postprocess1.begin(),
                                            need_postprocess1.end(), vic_ind),
                                            need_postprocess1.end());
      }
    }
  }
