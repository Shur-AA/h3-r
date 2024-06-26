#' Converts geographic coordinates into H3 indexes.
#'
#' @param lat
#' @param lon
#' @param res
#'
#' @return string vector of H3 indices
#' @export
h3_points_to_H3 <- function(lat, lon, res) {
  points_to_H3(lat, lon, res)
}


#' Defines H3 indexes and center points geographic coordinates
#' of hexagons which covers defined extent in defined level
#'
#' @param ext_lon
#' @param ext_lat
#' @param res
#'
#' @return map of H3 indexes and corresponding geographic coordinates
#' @export
h3_hex_centers_inbbox <- function(ext_lon, ext_lat, res) {
  hex_centers_inbbox(ext_lon, ext_lat, res)
}


#' Given parent H3 indexes and corresponding values defines
#' list of children in any level and maps parents' values
#'
#' @param level_to
#' @param parent_ind
#' @param parent_vals
#'
#' @return map of H3 indexes and corresponding z-values from parent cells
#' @export
h3_resample_down <- function(level_to, parent_ind, parent_vals) {
  resample_down(level_to, parent_ind, parent_vals)
}


#' Calculates sum of values in equal H3 cells of two rasters
#' - vectors with indexes and values; resamples one raster if
#' they have different levels
#'
#' @param ind1
#' @param z1
#' @param ind2
#' @param z2
#'
#' @return map of H3 indexes and corresponding z-values sum
#' @export
h3_simple_sum <- function(ind1, z1, ind2, z2) {
  simple_sum(ind1, z1, ind2, z2)
}


#' Gets centers of raster cells (stars object) and converts them into the
#' df with their coords and values
#'
#' @param raster_obj
#'
#' @return rast_pnt_grid with columns: ('lon_p', 'lat_p', 'value')
#' @export
h3_raster_to_points <- function(raster_obj) {
  raster_to_points(raster_obj)
}


#' Returns raster's extent coordinates in format
#' which fit to H3 functions
#'
#' @param raster_obj
#'
#' @return df with columns: ('lon_p', 'lat_p')
#' @export
h3_raster_to_bbox <- function(raster_obj) {
  raster_to_bbox(raster_obj)
}


#' Gets df with raster's centers coords and corresponding values
#' and interpolate them into the net of H3 hexagons' centers
#' which covers the raster with linear Akima interpolator
#'
#' @param rast_tab
#' @param hex_tab
#'
#' @return df with hex's centers coords and assigned values (x, y, z)
#' @export
h3_linear_interpolation <- function(rast_tab, hex_tab) {
  linear_interpolation(rast_tab, hex_tab)
}


#' Gets stars raster object, extracts values and interpolates them into
#' the net of H3 hexagons' centers which covers the raster
#' with bilinear interpolator from stars library
#'
#' @param raster_obj
#' @param hex_tab
#'
#' @return df with hex's centers coords and assigned values (z, x, y)
#' @export
h3_bilinear_interpolation <- function(raster_obj, hex_tab) {
  bilinear_interpolation(raster_obj, hex_tab)
}


#' Gets stars raster object and converts it into H3 indexes
#'
#' @param raster_obj
#' @param h3_level
#' @param intrp_method
#'
#' @return df with hex's centers coords, assigned raster values (z, x, y)
#' and H3 indexes
#'
#' @export
h3_raster_to_hex <- function(rast_obj, h3_level, intrp_method = 'bilinear') {
  raster_to_hex(rast_obj, h3_level, intrp_method = 'bilinear')
}


#' Gets stars raster object and returns its min, max or mean value
#'
#' @param ind
#' @param z
#' @param func
#'
#' @return raster min, max or mean value as double
#'
#' @export
h3_global_extremum <- function(ind, z, func) {
  global_extremum(ind, z, func)
}


#' Calculates sum, mean, min, max, majority or minority
#' of values in H3 cells of the
#' second raster which corresponds to H3 cells in the first
#' raster's zones, i.e. cells that share the same value;
#' resamples one raster if they have different levels
#'
#' @param zone_ind
#' @param zone_z
#' @param rast_ind
#' @param rast_z
#' @param stat_type
#'
#'
#' @return map of H3 indexes, corresponding zone codes from first raster
#' and statistics in each zone
#'
#' @export
h3_zonal_statistics <- function(zone_ind, zone_z, rast_ind, rast_z, stat_type, resample_zone) {
  zonal_statistics(zone_ind, zone_z, rast_ind, rast_z, stat_type, resample_zone)
}


#' Returns cell neighbours of defined order and less with the cell itself
#'
#' @param h3s
#' @param radius
#'
#'
#' @return cell neighbors of defined order and less with the cell itself
#'
#' @export
h3_cell_vecinity <- function(h3s, radius) {
  cell_vecinity(h3s, radius)
}


#' Returns hexs after simple focal operations, without matrices
#'
#' @param inds
#' @param z
#' @param stat_type
#' @param vecinity
#'
#' @return map of H3 indexes and corresponding z-values
#'
#' @export
h3_simple_focal <- function(inds, z, stat_type, vecinity) {
  simple_focal(inds, z, stat_type, vecinity)
}


#' Returns geographic coordinates for vector of hex indexes
#'
#' @param inds
#'
#' @return map of H3 indexes and corresponding z-values
#'
#' @export
h3_indexes_to_coords <- function(inds) {
  indexes_to_coords(inds)
}


#' Returns geodesic azimuth for a cell by its hex index
#'
#' @param h3_index string index of the cell
#'
#' @return double value of the azimuth
#'
#' @export
h3_cell_azimuth <- function(h3_index) {
  cell_azimuth(h3_index)
}


#' Returns table of input indexes and corresponding parent indexes
#'
#' @param children_ind vector of indexes
#'
#' @return table of indexes' pairs
#'
#' @export
h3_get_direct_parents <- function(children_ind) {
  get_direct_parents(children_ind)
}


#' Returns table of input indexes and corresponding parent indexes
#'
#' @param func aggregation function (sum, max, avg, majority)
#' @param children_ind vector of indexes
#' @param children_vals vector of values in the indexes
#'
#' @return map of parents' H3 indexes and corresponding z-values
#'
#' @export
h3_resample_up <- function(func, children_ind, children_vals) {
  resample_up(func, children_ind, children_vals)
}


#' Returns table of input indexes and corresponding parent indexes
#' of any coarser level
#'
#' @param func aggregation function (sum, max, avg, majority)
#' @param level_to aim level
#' @param children_ind vector of indexes
#' @param children_vals vector of values in the indexes
#'
#' @return map of parents' H3 indexes and corresponding z-values
#'
#' @export
h3_resample_up_any <- function(func, level_to, children_ind, children_vals) {
  resample_up_any(func, level_to, children_ind, children_vals)
}



#' Calculates water flow directions and fills depressions
#' needs h3 indexes and corresponding heights as well as center cell
#' of the area
#'
#' @param inds vector of H3 indexes
#' @param z corresponding heights
#' @param start_cell center cell of analyzing area
#'
#'
#' @return map of from-to H3 indexes - flow directions
#'
#' @export
h3_flow_dir <- function(inds, z, start_cell) {
  flow_dir(inds, z, start_cell)
}




#' Calculates water flow accumulation from flow direction table
#'
#'
#'
#' @param ifrom vector of H3 indexes
#' @param ito corresponding heights
#'
#'
#' @return map of H3 indexes and cell's accumulation value
#'
#' @export
h3_flow_acc <- function(ifrom, ito) {
  flow_acc(ifrom, ito)
}

