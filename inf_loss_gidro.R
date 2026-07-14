library(dplyr)
library(stars)
library(sf)
library(sp)
library(tidyverse)
library(akima)
library(tidyr)
library(readr)
library(bit64)


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


calc_entropy = function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  # Считаем вероятности (частоты) появления каждого направления
  p = table(x) / length(x)

  # Формула: -sum(p * log2(p))
  # Используем log2, чтобы получить результат в битах
  -sum(p * log(p, base = 2))
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
  colnames(t2) = c('h3_chld', 'z', 'h3_prnt')
  colnames(t1) = c('h3_chld', 'z', 'h3_prnt')

  hy = calc_entropy(t2$z)
  hx = calc_entropy(t1$z)
  t1 = t1 %>% filter(!is.na(t1$z))
  t2 = t2 %>% filter(!is.na(t2$z))

  # empirical joint distribution
  join_distr = full_join(t1, t2, by=c("h3_prnt" = "h3_chld")) %>%
    filter(!is.na(`z.x`)) %>%
    filter(!is.na(`z.y`)) %>%
    select(`z.x`, `z.y`)

  join_distr = join_distr %>%
    group_by(`z.x`, `z.y`) %>%
    summarise(pair_cnt = n())
  s = sum(join_distr$pair_cnt)
  join_distr = join_distr %>%
    mutate(pair_p = pair_cnt / s)

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

  a = select(vp, z, pij.x)
  b = select(vp, z, pij.y) %>%
    filter(pij.y > 0)
  c = crossing(a$z, b$z)
  colnames(c) = c('xval', 'yval')
  c = left_join(c, a, by=c('xval' = 'z')) %>%
    left_join(b, by=c('yval' = 'z')) %>%
    mutate(pxy = pij.x * pij.y)

  # table with joints and marginal distribution
  distributions = left_join(join_distr, vp, by=c('z.x' = 'z'))
  for_hxy = filter(distributions, pij.y > 0)

  # empirical conditional entropy H(Y|X)
  hyx_emp = -sum(distributions$pair_p * log2(distributions$pair_p / distributions$pij.x))

  # theoretical conditional entropy H(Y|X)
  hyx_teor = -sum(c$pxy * log2(c$pxy / c$pij.x))

  # I(X,Y) = H(Y) - H(Y|X)
  MIemp = hy - hyx_emp
  MIteor = hy - hyx_teor
  # I(X,Y) = H(X) + H(Y) - H(X,Y)
  MI2 = sum(distributions$pair_p * log2(distributions$pair_p)) + hx + hy
  # I(X,Y) = sum[P(X,Y) * log2(P(X,Y) / P(X)*P(Y))]
  MI3 = sum(for_hxy$pair_p * log2(for_hxy$pair_p / (for_hxy$pij.x * for_hxy$pij.y)))

  mi = data.frame(MI1_emp = c(MIemp),
                  MI1_theor = c(MIteor),
                  MI2 = c(MI2),
                  MI3 = c(MI3))

  return(mi)
}


mut_inf2 = function(t1, t2){
  # Оставляем только нужные колонки и чистим NA
  t1 = t1[, c(1, 2, 3)] %>% filter(!is.na(.[,2]))
  t2 = t2[, c(1, 2, 3)] %>% filter(!is.na(.[,2]))
  colnames(t1) = c('h3_chld', 'z', 'h3_prnt')
  colnames(t2) = c('h3_chld', 'z', 'h3_prnt')

  # 1. Строим совместное распределение (только пересекающиеся данные)
  join_distr = inner_join(t1, t2, by=c("h3_prnt" = "h3_chld")) %>%
    rename(z_fine = z.x, z_coarse = z.y) %>%
    group_by(z_fine, z_coarse) %>%
    summarise(pair_cnt = n(), .groups = 'drop') %>%
    mutate(p_xy = pair_cnt / sum(pair_cnt))

  # 2. Считаем маржинальные вероятности строго из совместного распределения
  p_x = join_distr %>% group_by(z_fine) %>% summarise(p_x = sum(p_xy), .groups = 'drop')
  p_y = join_distr %>% group_by(z_coarse) %>% summarise(p_y = sum(p_xy), .groups = 'drop')

  # Добавляем маржинальные вероятности в общую таблицу
  distr_full = join_distr %>%
    left_join(p_x, by = "z_fine") %>%
    left_join(p_y, by = "z_coarse")

  # 3. Энтропии
  hx = -sum(p_x$p_x * log2(p_x$p_x))
  hy = -sum(p_y$p_y * log2(p_y$p_y))
  h_xy = -sum(distr_full$p_xy * log2(distr_full$p_xy))

  # 4. Расчет вариантов MI
  # I(X,Y) = H(Y) - H(Y|X)
  # H(Y|X) = H(X,Y) - H(X)
  h_y_cond_x = h_xy - hx
  MI1 = hy - h_y_cond_x

  # I(X,Y) = H(X) + H(Y) - H(X,Y)
  MI2 = hx + hy - h_xy

  # I(X,Y) = sum[P(x,y) * log2(P(x,y) / (P(x)*P(y)))]
  MI3 = sum(distr_full$p_xy * log2(distr_full$p_xy / (distr_full$p_x * distr_full$p_y)))

  return(data.frame(MI_cond = MI1, MI_diff = MI2, MI_kl = MI3, h_yx = h_y_cond_x))
}

calc_kl_divergence = function(df_fine, df_coarse) {
  colnames(df_fine) = c('h3_chld', 'z', 'h3_prnt')
  colnames(df_coarse) = c('h3_chld', 'z', 'h3_prnt')
  # 1. Получаем распределение (гистограмму) направлений для обоих уровней
  # Предполагаем, что z — это колонка с кодом направления (1-6)
  dist_p = df_fine %>%
    count(z) %>%
    mutate(p = n / sum(n))

  dist_q = df_coarse %>%
    count(z) %>%
    mutate(q = n / sum(n))

  # 2. Объединяем, чтобы категории совпадали
  combined = full_join(dist_p, dist_q, by = "z") %>%
    mutate(p = ifelse(is.na(p), 0, p),
           q = ifelse(is.na(q), 1e-10, q)) # Добавляем epsilon для стабильности

  # 3. Считаем KL: сумма P(i) * log2(P(i)/Q(i))
  # Считаем только там, где P(i) > 0
  kl_val = combined %>%
    filter(p > 0) %>%
    summarise(kl = sum(p * log2(p / q))) %>%
    pull(kl)

  return(kl_val)
}

# Пример использования:
# kl_nsa = calc_kl_divergence(data_res8, data_nsa7)
# kl_vvrfra = calc_kl_divergence(data_res8, data_vvrfra7)






# папка с шейпами по населению в городах
dir = 'D:/1_Аспа/Information loss/gidroinfo'



tab0 <- read_csv(paste(dir, '/', 'new_dirs/8_dirs.csv', sep=''), col_types = cols(
  from = col_character(),
  direction = col_integer()
))
tab0 = as.data.frame(tab0)



tab1 <- read_csv(paste(dir, '/', 'new_dirs/vv7_dirs.csv', sep=''), col_types = cols(
  from = col_character(),
  direction = col_integer()
))
tab1 = as.data.frame(tab1)

tab00 = select(tab0, from, direction, parent)

pkl = calc_kl_divergence(tab00, tab1)

pmi = mut_inf2(tab00, tab1)
pmi
pkl





ent = tab0 %>% group_by(tab0$parent) %>%
  summarise(n_cells = n(),                     # Количество ячеек в группе
            ent = calc_entropy(direction)
  ) %>%
  ungroup()

sr = calc_entropy(tab0$direction)
p10 = mean(ent$ent)

sr
p10






# ENTROPY CORRECTION
p = table(tab0$direction) / length(tab0$direction)
if (p[[3]] > p[[1]]){
  kk = c(p[[1]] + abs(p[[1]] - p[[3]])*dolia)
} else {
  kk = c(p[[1]] - abs(p[[1]] - p[[3]])*dolia)
}
if (p[[6]] > p[[2]]){
  kk = append(kk, p[[2]] + abs(p[[2]] - p[[6]])*dolia)
} else {
  kk = append(kk, p[[2]] - abs(p[[2]] - p[[6]])*dolia)
}
if (p[[2]] > p[[3]]){
  kk = append(kk, p[[3]] + abs(p[[3]] - p[[2]])*dolia)
} else {
  kk = append(kk, p[[3]] - abs(p[[3]] - p[[2]])*dolia)
}
if (p[[5]] > p[[4]]){
  kk = append(kk, p[[4]] + abs(p[[4]] - p[[5]])*dolia)
} else {
  kk = append(kk, p[[4]] - abs(p[[4]] - p[[5]])*dolia)
}
if (p[[1]] > p[[5]]){
  kk = append(kk, p[[5]] + abs(p[[5]] - p[[1]])*dolia)
} else {
  kk = append(kk, p[[5]] - abs(p[[5]] - p[[1]])*dolia)
}
if (p[[4]] > p[[6]]){
  kk = append(kk, p[[6]] + abs(p[[6]] - p[[4]])*dolia)
} else {
  kk = append(kk, p[[6]] - abs(p[[6]] - p[[4]])*dolia)
}

-sum(p * log(p, base = 2))
-sum(kk * log(kk, base = 2))


# Kl CORRECTION
tab0 <- read_csv(paste(dir, '/', 'new_dirs/8_dirs.csv', sep=''), col_types = cols(
  from = col_character(),
  direction = col_integer()
))
tab0 = as.data.frame(tab0)



tab1 <- read_csv(paste(dir, '/', 'new_dirs/nsa3_dirs.csv', sep=''), col_types = cols(
  from = col_character(),
  direction = col_integer()
))
tab1 = as.data.frame(tab1)

tab00 = select(tab0, from, direction, parent3)



p = table(tab0$direction) / length(tab0$direction)
colnames(tab00) = c('h3_chld', 'z', 'h3_prnt')
colnames(tab1) = c('h3_chld', 'z', 'h3_prnt')
# 1. Получаем распределение (гистограмму) направлений для обоих уровней
# Предполагаем, что z — это колонка с кодом направления (1-6)
dist_p = tab00 %>%
  count(z) %>%
  mutate(p = n / sum(n))

dist_q = tab1 %>%
  count(z) %>%
  mutate(q = n / sum(n))

# 2. Объединяем, чтобы категории совпадали
combined = full_join(dist_p, dist_q, by = "z") %>%
  mutate(p = ifelse(is.na(p), 0, p),
         q = ifelse(is.na(q), 1e-10, q)) # Добавляем epsilon для стабильности




constanta = 19.10661
dolia = constanta/60


if (p[[3]] > p[[1]]){
  kk = c(combined$n.y[[1]] + abs(combined$n.y[[1]] - combined$n.y[[3]])*dolia)
} else {
  kk = c(combined$n.y[[1]] - abs(combined$n.y[[1]] - combined$n.y[[3]])*dolia)
}
if (combined$n.y[[6]] > combined$n.y[[2]]){
  kk = append(kk, combined$n.y[[2]] + abs(combined$n.y[[2]] - combined$n.y[[6]])*dolia)
} else {
  kk = append(kk, combined$n.y[[2]] - abs(combined$n.y[[2]] - combined$n.y[[6]])*dolia)
}
if (combined$n.y[[2]] > combined$n.y[[3]]){
  kk = append(kk, combined$n.y[[3]] + abs(combined$n.y[[3]] - combined$n.y[[2]])*dolia)
} else {
  kk = append(kk, combined$n.y[[3]] - abs(combined$n.y[[3]] - combined$n.y[[2]])*dolia)
}
if (combined$n.y[[5]] > combined$n.y[[4]]){
  kk = append(kk, combined$n.y[[4]] + abs(combined$n.y[[4]] - combined$n.y[[5]])*dolia)
} else {
  kk = append(kk, combined$n.y[[4]] - abs(combined$n.y[[4]] - combined$n.y[[5]])*dolia)
}
if (combined$n.y[[1]] > combined$n.y[[5]]){
  kk = append(kk, combined$n.y[[5]] + abs(combined$n.y[[5]] - combined$n.y[[1]])*dolia)
} else {
  kk = append(kk, combined$n.y[[5]] - abs(combined$n.y[[5]] - combined$n.y[[1]])*dolia)
}
if (combined$n.y[[4]] > combined$n.y[[6]]){
  kk = append(kk, combined$n.y[[6]] + abs(combined$n.y[[6]] - combined$n.y[[4]])*dolia)
} else {
  kk = append(kk, combined$n.y[[6]] - abs(combined$n.y[[6]] - combined$n.y[[4]])*dolia)
}

combined = cbind(combined, kk) %>% mutate(qq = kk/sum(kk))

kl_val = combined %>%
  filter(p > 0) %>%
  summarise(kl = sum(p * log2(p / q))) %>%
  pull(kl)

kl_val_new = combined %>%
  filter(p > 0) %>%
  summarise(kl = sum(p * log2(p / qq))) %>%
  pull(kl)
kl_val
kl_val_new




# Graphs
angle = c(293, 53, 353, 173, 233, 113)
angle1 = angle - 19
t = table(tab0$direction) %>%
  as.data.frame() %>% mutate(pr = Freq/sum(Freq)*100) %>%
  cbind(angle)
names(t) = c("Xh3", "Yabs", "Y", "X")




ggplot(t, aes(x = X, y = Y, fill = Y)) +
  # geom_bar отрисовывает "лепестки"
  # width подбирается в зависимости от шага азимута (например, 45)
  geom_bar(stat = "identity", width = 12, color = "white") +

  # Превращаем в полярные координаты
  coord_polar(theta = "x", start = 0, direction = 1) +
  # Настройка шкалы X, чтобы 0 был сверху и метки шли по кругу
  scale_x_continuous(
    limits = c(0, 360),
    breaks = seq(0, 330, by = 30),
    labels = c("0°", "30°", "60°", "90°", "120°", "150°", "180°", "210°", "240°", "270°", "300°", "330°")
  ) +
  labs(
    title = "Алгоритм VVRFRA, уровень детализации 2",
    fill = "Частота (%)
    ",    # Название для легенды заливки
  )+
  theme_minimal() +
  theme(
    plot.title = element_text(size = 25, face = "bold"),
    axis.text.y = element_blank(),
    # Подписи осей (названия "Azimuth", "Percent")
    axis.title = element_blank(),

    # Текст меток на осях (числа 0, 45, 90... или "С", "В")
    axis.text = element_text(size = 23, face = "bold"),

    # Текст легенды (если она есть)
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )+
  geom_text(
    data = data.frame(y = c(5, 10, 15, 20, 25)), # Ваши пороги процентов
    aes(x = 77, y = y, label = paste0(y, "%")), # x = 20 — это угол в градусах
    inherit.aes = FALSE,
    size = 6
  )+
  scale_fill_distiller(palette = 'PuBu', limits = c(5, 26), direction = 1)





t |>
  ggplot(mapping = aes(x = X, y = Y, fill = Y)) +
  geom_col() +
  coord_polar()
