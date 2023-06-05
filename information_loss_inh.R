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

entropy = function(intab, val_col){
  intab = intab %>% filter(!is.na(intab[,val_col]))
  cells_count = length(intab[,val_col])
  entr = intab %>% group_by(intab[,val_col]) %>%
    summarise(val_cnt = n()) %>%
    mutate(pij = val_cnt / (cells_count * 1.0)) %>%
    mutate(log2pij = log2(pij)) %>%
    mutate(pijlog2pij = pij * log2pij)
  entr = sum(entr$pijlog2pij)
  return(-entr)
}


# calculates Kullback–Leibler divergence between two vectors of variables
# note that v1 is reality and v2 is a model
KL_divergence = function(v1, v2){
  v1 = v1[!is.na(v1)] %>% as.data.frame()
  v2 = v2[!is.na(v2)] %>% as.data.frame()
  colnames(v1) = c('val')
  colnames(v2) = c('val')
  v1_count = length(v1$val)
  v2_count = length(v2$val)
  v1p = v1 %>% group_by(val) %>%
    summarise(val_cnt = n()) %>%
    mutate(pij = val_cnt / (v1_count * 1.0))
  v2p = v2 %>% group_by(val) %>%
    summarise(val_cnt = n()) %>%
    mutate(pij = val_cnt / (v2_count * 1.0))
  vp = full_join(v1p, v2p, by='val')
  vp[is.na(vp)] = -1
  vp = vp %>% mutate(kl = pij.x * log2(pij.x / pij.y))
  return(sum(filter(vp, !is.nan(vp$kl))$kl))
}


# mutual information in 3 notation
# https://jinjeon.me/post/mutual-info/
# https://medium.com/swlh/a-deep-conceptual-guide-to-mutual-information-a5021031fad0
# (X - reality, Y - model)
# tables must have columns 'h3_ind' and 'z'
mut_inf = function(t1, t2){
  hy = entropy(as.data.frame(t2$z), 1)
  hx = entropy(as.data.frame(t1$z), 1)
  t1 = t1 %>% filter(!is.na(t1$z))
  t2 = t2 %>% filter(!is.na(t2$z))
  # t1 indexes and their parents
  print(11)
  t1p = h3_get_direct_parents(t1$h3_ind) %>% as.data.frame()
  t1p = mutate(t1p, h3_chld = rownames(t1p))
  rownames(t1p) = seq(1, length(t1p$h3_chld), 1)
  colnames(t1p) = c('h3_prnt', 'h3_chld')
  print(22)
  # empirical joint distribution
  join_distr = full_join(t1, t1p, by=c("h3_ind" = "h3_chld")) %>%
    full_join(t2, by=c("h3_prnt" = "h3_ind")) %>%
    filter(!is.na(`z.x`)) %>%
    filter(!is.na(`z.y`)) %>%
    select(`z.x`, `z.y`)
  print(33)
  join_distr = join_distr %>%
    group_by(`z.x`, `z.y`) %>%
    summarise(pair_cnt = n())
  s = sum(join_distr$pair_cnt)
  join_distr = join_distr %>%
    mutate(pair_p = pair_cnt / s)
  print(44)
  # theoretical joint distribution
  t1_count = length(t1$z)
  t2_count = length(t2$z)
  px = t1 %>% group_by(z) %>%
    summarise(val_cnt = n()) %>%
    mutate(pij = val_cnt / (t1_count * 1.0))
  py = t2 %>% group_by(z) %>%
    summarise(val_cnt = n()) %>%
    mutate(pij = val_cnt / (t2_count * 1.0))
  vp = full_join(px, py, by='z')
  vp[is.na(vp)] = 0
  print(55)
  a = select(vp, z, pij.x)
  b = select(vp, z, pij.y) %>%
    filter(pij.y > 0)
  c = crossing(a$z, b$z)
  colnames(c) = c('xval', 'yval')
  c = left_join(c, a, by=c('xval' = 'z')) %>%
    left_join(b, by=c('yval' = 'z')) %>%
    mutate(pxy = pij.x * pij.y)
  print(66)
  # table with joints and marginal distribution
  distributions = left_join(join_distr, vp, by=c('z.x' = 'z'))
  for_hxy = filter(distributions, pij.y > 0)
  print(77)
  # empirical conditional entropy H(Y|X)
  hyx_emp = -sum(distributions$pair_p * log2(distributions$pair_p / distributions$pij.x))
  print(88)
  # theoretical conditional entropy H(Y|X)
  hyx_teor = -sum(c$pxy * log2(c$pxy / c$pij.x))
  print(99)
  # I(X,Y) = H(Y) - H(Y|X)
  MIemp = hy - hyx_emp
  MIteor = hy - hyx_teor
  # I(X,Y) = H(X) + H(Y) - H(X,Y)
  MI2 = sum(distributions$pair_p * log2(distributions$pair_p)) + hx + hy
  # I(X,Y) = sum[P(X,Y) * log2(P(X,Y) / P(X)*P(Y))]
  MI3 = sum(for_hxy$pair_p * log2(for_hxy$pair_p / (for_hxy$pij.x * for_hxy$pij.y)))
  print(110)
  mi = data.frame(MI1_emp = c(MIemp),
                  MI1_theor = c(MIteor),
                  MI2 = c(MI2),
                  MI3 = c(MI3))

  return(mi)
}


spectr_dev = function(t1, t2){
  t1 = t1 %>% filter(!is.na(t1$z))
  t2 = t2 %>% filter(!is.na(t2$z))
  # t1 indexes and their parents
  t1p = h3_get_direct_parents(t1$h3_ind) %>% as.data.frame()
  t1p = mutate(t1p, h3_chld = rownames(t1p))
  rownames(t1p) = seq(1, length(t1p$h3_chld), 1)
  colnames(t1p) = c('h3_prnt', 'h3_chld')
  # empirical joint distribution
  join_distr = full_join(t1, t1p, by=c("h3_ind" = "h3_chld")) %>%
    full_join(t2, by=c("h3_prnt" = "h3_ind")) %>%
    filter(!is.na(`z.x`)) %>%
    filter(!is.na(`z.y`)) %>%
    select(`z.x`, `z.y`)
  colnames(join_distr) = c('a', 'b')
  return(sum(abs(join_distr$b - join_distr$a) / join_distr$a) / length(join_distr$a))
}




# КАРУСЕЛЬ
#
# папка с папками lclu
dir = 'D:/1_Аспа/Information loss/hres_crop'
products = list.dirs(dir)

# составляем список уникальных названий городов в папках,
# чтобы потом обращаться именно к городам
citytifs = c()
foldernames = c()
for (product in 2:length(products)){
  p = strsplit(products[product],split='/', fixed=TRUE)
  p = p[[1]][length(p[[1]])] # название продукта
  foldernames = append(foldernames, p, after = length(foldernames))
  # start_symb = nchar(p) + 2
  start_symb = 1
  tifs_infolder = list.files(products[product], pattern = '*.tif$')
  for (tif in tifs_infolder){
    end_symb = nchar(tif) - 4
    cityname = substr(tif, start_symb, end_symb)
    if (cityname %in% citytifs == F){
      citytifs = append(citytifs, cityname, after = length(citytifs))
    }
  }
}

result = data.frame(territory = c("city"),
                    rast_prod = c("territory"),
                    start_res = c(99),
                    start_h3_level = c(99),
                    h3_level_n = c(99),
                    gen_method = c('start'),
                    aggr_method = c('nn'),
                    legend_variety = c(0),
                    entropy = c(0),
                    mi1 = c(0),
                    mi2 = c(0),
                    mi3 = c(0),
                    kl_div = c(0),
                    spectr_dev = c(0))


# ходим по городам
for (city in citytifs){
  print(city)
  current_city_paths = c()
  for (folder in foldernames){
    # path_to_city = paste(dir, '/', folder, '/', folder, '_', city, '.tif', sep='')
    path_to_city = paste(dir, '/', folder, '/', city, '.tif', sep='')
    current_city_paths = append(current_city_paths,
                                path_to_city,
                                after = length(current_city_paths))
  }
  ccp_length = length(current_city_paths)
  for (raster_num in 1:ccp_length){
    print(current_city_paths[raster_num])
    rast1 = read_stars(current_city_paths[raster_num]) # input 1-band raster file


    # information for table and parameters
    initial_resolution = get_rast_cellarea(rast1)
    start_h3_l = choose_h3_level(initial_resolution)
    w = strsplit(current_city_paths[raster_num], split='/')
    territory = substr(w[[1]][length(w[[1]])], 1, nchar(w[[1]][length(w[[1]])]) - 4)
    # сколько уровней будем проходить; 4 - самый высокий уровень, эксп. мнение
    levels_seq = seq(start_h3_l, 5, -1)

    start_tab = h3_raster_to_hex(rast1, start_h3_l, 'bilinear')
    start_tab$z = as.integer(start_tab$z)
    print('transition done')

    p01 = entropy(start_tab, 1)
    p00 = length(unique(start_tab$z))
    this_result = data.frame(territory = c(city),
                             rast_prod = c(territory),
                             start_res = c(initial_resolution),
                             start_h3_level = c(start_h3_l),
                             h3_level_n = c(start_h3_l),
                             gen_method = c('start'),
                             aggr_method = c('bilinear'),
                             legend_variety = c(p00),
                             entropy = c(p01),
                             mi1 = c(0),
                             mi2 = c(0),
                             mi3 = c(0),
                             kl_div = c(0),
                             spectr_dev = c(0)
    )
    result = rbind(result, this_result)

    for (level in levels_seq){
      print(level)
      if (level == start_h3_l){
        coarse_tab = start_tab
        gc(start_tab)
        rm(start_tab)
      } else {
        coarse_tab = finer_tab_inh
        gc(finer_tab_inh)
        rm(finer_tab_inh)
      }
      print('changing level')
      # parent resolution via aggregation

      finer_tab_inh = h3_resample_up("sum",
                                     coarse_tab$h3_ind,
                                     coarse_tab$z) %>%
        as.data.frame()
      finer_tab_inh = mutate(finer_tab_inh,
                             h3_ind = rownames(finer_tab_inh)) %>%
        filter(. > 0)
      rownames(finer_tab_inh) = seq(1, length(finer_tab_inh$h3_ind), 1)
      colnames(finer_tab_inh) = c('z', 'h3_ind')

      print('calculate')
      # 1 - finer_inh, 2 - finer_indep;
      # 0 - variety, 1 - entropy, 2 - mi1, 3 - mi2, 4 - mi3, 5 - KL, 6 - spectr_dev
      p10 = length(unique(finer_tab_inh$z))
      p11 = entropy(finer_tab_inh, 1)
      mimi = mut_inf(coarse_tab, finer_tab_inh)
      p12 = mimi[1, 1]
      p13 = mimi[1, 3]
      p14 = mimi[1, 4]

      p15 = KL_divergence(coarse_tab$z, finer_tab_inh$z)
      p16 = spectr_dev(coarse_tab, finer_tab_inh)

      this_result = data.frame(territory = c(city),
                               rast_prod = c(territory),
                               start_res = c(initial_resolution),
                               start_h3_level = c(start_h3_l),
                               h3_level_n = c(level - 1),
                               gen_method = c('inheritance'),
                               aggr_method = c('sum'),
                               legend_variety = c(p10),
                               entropy = c(p11),
                               mi1 = c(p12),
                               mi2 = c(p13),
                               mi3 = c(p14),
                               kl_div = c(p15),
                               spectr_dev = c(p16)
      )
      result = rbind(result, this_result)
    }
  }
}
