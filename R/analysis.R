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
h3_raster_to_hex <- function(rast_obj, h3_level, intrp_method) {
  raster_to_hex(rast_obj, h3_level, intrp_method)
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


#' Calculates sum of values in H3 cells of the second raster which
#' corresponds to H3 cells in the first raster's zones, i.e. cells
#' that share the same value; resamples one raster if
#' they have different levels
#'
#' @param zone_ind
#' @param zone_z
#' @param rast_ind
#' @param rast_z
#'
#' @return map of H3 indexes, corresponding zone codes from first raster
#' and sums in each zone
#'
#' @export
h3_zonal_sum <- function(zone_ind, zone_z, rast_ind, rast_z) {
  zonal_sum(zone_ind, zone_z, rast_ind, rast_z)
}


#' Calculates sum, mean, min, or max of values in H3 cells of the
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
