//sum.cpp
#include <Rcpp.h>
#include "h3api.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>


using namespace Rcpp;


// internal functions

// converts degrees to radians
double deg_to_rad(const double & degs){
  const double pi = 3.14159265358979323846; // the value of pi
  double radians = degs * pi / 180.0;
  return radians;
}


// checks whether the lon value is in right interval
double check_lon(double longitude) {
  if ((longitude > 180) && (longitude <= 360)) {
    return longitude - 360;
  }
  return longitude;
}


// checks whether the lat value is in right interval
double check_lat(double latitude) {
  if (latitude > 90) {
    return latitude - 180;
  }
  return latitude;
}


// 'exponent it'
double expit(double x, int power){
  double result = 1.0;
  for (int i = 1; i <= power; i++){
    result = result * x;
  }
  return result;
}


// for equal levels only
// sums two vectors z values with corresponding indeces ind
std::map <std::string, double> sum_rasters(const std::vector<std::string> & ind1,
                                           const std::vector<double> & z1,
                                           const std::vector<std::string> & ind2,
                                           const std::vector<double> & z2){
    std::map <std::string, double> resulttab;
    int intersection_num = 0;
    for (int i = 0; i < ind1.size(); i++){
      for (int j = 0; j < ind2.size(); j++){
        if (ind1[i] == ind2[j]){
          intersection_num ++;
          resulttab[ind1[i]] = z1[i] + z2[j];
        }
      }
    }
    if (intersection_num == 0){
      resulttab["result"] = 0;
    }
    return resulttab;
}


double vector_average(std::vector<double> & v){
  int n = v.size();
  double acc = 0.0;
  for (auto i : v){
    if (!std::isnan(i)){
      acc += i;
    }}
  return acc / n;
}


double most_frequent_element(std::vector<double> &arr)
  {
    if (arr.empty())
      return -1;

    std::map<double,double> mydict = {};
    int cnt = 0;
    int itm = 0;

    for (auto&& item : arr) {
      mydict[item] = mydict.emplace(item, 0).first->second + 1;
      if (mydict[item] >= cnt) {
        std::tie(cnt, itm) = std::tie(mydict[item], item);
      }
    }
    return itm;
  }


double least_frequent_element(std::vector<double> &arr)
{
  if (arr.empty())
    return -1;

  int leastCtr = arr.size();
  std::vector<double> temp_arr;
  for(int i = 0; i < arr.size(); i++)
  {temp_arr.push_back(arr[i]);}
  double leastElement = -99;
  double currentCtr =  1;

  std::sort(temp_arr.begin(), temp_arr.end());
  for (int i = 0; i < (temp_arr.size()-1); i++){
    if (temp_arr[i] == temp_arr[i + 1])
    {currentCtr += 1;}
    else{
      if (currentCtr < leastCtr){
        leastCtr = currentCtr;
        leastElement = temp_arr[i];
        currentCtr = 1;}}
  }

  if (currentCtr < leastCtr){
    leastCtr = currentCtr;
    leastElement = temp_arr[temp_arr.size() - 1];
  }
  return leastElement ;
}

// [[Rcpp::export]]
std::string H3_to_parent(std::string h3s,
                         int res) {
  std::string z;
  H3Index h3 = stringToH3(h3s.std::string::c_str());
  H3Index h3Parent = h3ToParent(h3, res);
  char h3ParentStr[17];
  h3ToString(h3Parent, h3ParentStr, sizeof(h3ParentStr));
  z = h3ParentStr;
  return z;
}


// END internal functions







// converts geographic coordinates into H3 indexes

// [[Rcpp::export]]
std::vector<std::string> points_to_H3(const std::vector<double> & lon,
                                      const std::vector<double> & lat,
                                      const int & res) {
  int coord_size = lat.size();
  std::vector<std::string> indexes(coord_size);
  for (int i = 0; i < coord_size; ++i){
    GeoCoord location;
    location.lat = deg_to_rad(check_lat(lat[i]));
    location.lon = deg_to_rad(check_lon(lon[i]));
    H3Index h3 = geoToH3(&location, res);
    char h3s[17];
    std::cout << h3 << std::endl;
    h3ToString(h3, h3s, sizeof(h3s));
    indexes[i] = h3s;
  }
  return indexes;
}


// defines H3 indexes and center points geographic coordinates
// of hexagons which covers defined extent in defined level
// returns H3 indexes and corresponding geo coords

// [[Rcpp::export]]
std::map <std::string, std::vector<double>> hex_centers_inbbox(const std::vector<double> & ext_lon,
                                                               const std::vector<double> & ext_lat,
                                                               const int & res){
  // number of vertices in the extent
  int vert_num = ext_lon.size();
  // converting extent vertices coords in
  // special polygon structure which H3 understands

  // make geofence from geocoord
  Geofence geofence;
  geofence.numVerts= vert_num;
  GeoCoord* sfVerts = new GeoCoord[vert_num];
  for (int i = 0; i < vert_num; i++) {
    GeoCoord coord;
    coord.lat = deg_to_rad(check_lat(ext_lat[i]));
    coord.lon = deg_to_rad(check_lon(ext_lon[i]));
    sfVerts[i] = coord;
  }
  geofence.verts = sfVerts;

  // Make polygon
  GeoPolygon extent;
  extent.geofence = geofence;
  extent.numHoles = 0;

  // h3 indexes in extent
  int numHexagons = maxPolyfillSize(&extent, res);

  H3Index* hexagons = new H3Index[numHexagons]();
  polyfill(&extent, res, hexagons);

  std::map <std::string, std::vector<double>> geotab;

  for (int i = 0; i < numHexagons; i++) {
    H3Index hexagon = hexagons[i];
    if (hexagon != 0) {
      char h3Str[17];
      h3ToString(hexagons[i], h3Str, sizeof(h3Str));
      GeoCoord hex_center;
      std::vector<double> centroid_coords;
      h3ToGeo(hexagons[i], &hex_center);
      centroid_coords.push_back(radsToDegs(hex_center.lon));
      centroid_coords.push_back(radsToDegs(hex_center.lat));
      geotab[h3Str] = centroid_coords;
    }
  }
  delete[] hexagons;
  return geotab;
}


// given parent H3 indexes and corresponding values defines
// list of children in any level and maps parents' values

// [[Rcpp::export]]
std::map <std::string, double> resample_down(const int & level_to,
                                             const std::vector<std::string> & parent_ind,
                                             const std::vector<double> & parent_vals){
  int level_from = h3GetResolution(stringToH3(parent_ind[0].std::string::c_str()));

  try{
    if (level_from == level_to){
      throw 0; // nothing to do exception
    }}
  catch(int x){
    std::cout<<"nothing to do exception - equal levels" << std::endl;
  }

  try{
    if (level_from > level_to){
      throw 1; // wrong direction of down resampling exception
    }}
  catch(int x){
    std::cout<<"wrong direction of down resampling exception" << std::endl;
  }

  try{
    if (parent_ind.size() != parent_vals.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal lengths exception - unpredictable result" << std::endl;
  }


  std::map <std::string, double> childtab;

  // loop through parent string indexes
  for (int i = 0; i < parent_ind.size(); i++) {
    H3Index h3 = stringToH3(parent_ind[i].std::string::c_str());
    int n = maxH3ToChildrenSize(h3, level_to);  // define maximum number of the hex's children
    H3Index* h3Children = new H3Index[n];  // vector for children
    h3ToChildren(h3, level_to, h3Children);
    for (int j = 0; j < n; ++j) {
      char h3Str[17];
      h3ToString(h3Children[j], h3Str, sizeof(h3Str));
      childtab[h3Str] = parent_vals[i];
    }
    delete[] h3Children;
  }
  return childtab;
}


// calculates sum of values in equal H3 cells of two rasters
// (vectors with indexes and values); resamples one raster if
// they have different levels

// [[Rcpp::export]]
std::map <std::string, double> simple_sum(std::vector<std::string> & ind1,
                                          std::vector<double> & z1,
                                          std::vector<std::string> & ind2,
                                          std::vector<double> & z2){
  int h3_level1 = h3GetResolution(stringToH3(ind1[0].std::string::c_str()));
  int h3_level2 = h3GetResolution(stringToH3(ind2[0].std::string::c_str()));

  try{
  if ((ind1.size() != z1.size()) || (ind2.size() != z2.size())){
    throw 2; // not equal lengths exception
  }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  // проверить, как себя ведёт с NA и NaN
  // написать обработчики ошибок (почти)

  if (h3_level1 > h3_level2){
    std::map <std::string, double> new_ind2 = resample_down(h3_level1, ind2, z2);
    ind2.resize(new_ind2.size());
    z2.resize(new_ind2.size());
    int k = 0;
    for (auto const & pair : new_ind2){
      ind2[k] = pair.first;
      z2[k] = pair.second;
      k++;
    }
  }
  if (h3_level1 < h3_level2){
    std::map <std::string, double> new_ind1 = resample_down(h3_level2, ind1, z1);
    ind1.resize(new_ind1.size());
    z1.resize(new_ind1.size());
    int k = 0;
    for (auto const & pair : new_ind1){
      ind1[k] = pair.first;
      z1[k] = pair.second;
      k++;
    }
  }

  std::map <std::string, double> sumtab = sum_rasters(ind1, z1, ind2, z2);

  return sumtab;

  // продумать разные типы данных для z (перегрузка, шаблоны?)
}



// example of a global function: returns max, min or mean value from raster

// [[Rcpp::export]]
double global_extremum(std::vector<std::string> & ind,
                       std::vector<double> & z,
                       const std::string func){
  try{
    if (ind.size() != z.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  double result;
  if (func == "min"){
    result = *min_element(std::begin(z), std::end(z));
  } else if (func == "max"){
    result = *max_element(std::begin(z), std::end(z));
  } else if (func == "avg"){
    result = vector_average(z);
  } else {
    result = -0.0;
  }
  return result;
}


// Calculates sum, mean, min, or max of values in H3 cells of the
// second raster which corresponds to H3 cells in the first
// raster's zones, i.e. cells that share the same value;
// resamples one raster if they have different levels

// [[Rcpp::export]]
std::map <std::string, std::vector<double>> zonal_statistics(
                                          std::vector<std::string> & zone_ind,
                                          std::vector<double> & zone_z,
                                          std::vector<std::string> & rast_ind,
                                          std::vector<double> & rast_z,
                                          const std::string stat_type,
                                          const bool resample_zone){
  int h3_level_zone = h3GetResolution(stringToH3(zone_ind[0].std::string::c_str()));
  int h3_level_rast = h3GetResolution(stringToH3(rast_ind[0].std::string::c_str()));

  try{
    if ((zone_ind.size() != zone_z.size()) || (rast_ind.size() != rast_z.size())){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  // проверить, как себя ведёт с NA и NaN
  // написать обработчики ошибок (почти)

  // level relationships: 0 - equal; 1 - raster's coarser; 2 - zones're coarser
  int resample_status = 0;

  // if zones have finer resolution - resample raster
  if (h3_level_zone > h3_level_rast){
    resample_status = 1;
    std::map <std::string, double> new_rast_ind = resample_down(h3_level_zone,
                                                                rast_ind,
                                                                rast_z);
    rast_ind.resize(new_rast_ind.size());
    rast_z.resize(new_rast_ind.size());
    int k = 0;
    for (auto const & pair : new_rast_ind){
      rast_ind[k] = pair.first;
      rast_z[k] = pair.second;
      k++;
    }
  }

  // otherwise - resample zones
  if (h3_level_zone < h3_level_rast){
    resample_status = 2;
    std::map <std::string, double> new_zone_ind = resample_down(h3_level_rast,
                                                                zone_ind,
                                                                zone_z);
    zone_ind.resize(new_zone_ind.size());
    zone_z.resize(new_zone_ind.size());
    int k = 0;
    for (auto const & pair : new_zone_ind){
      zone_ind[k] = pair.first;
      zone_z[k] = pair.second;
      k++;
    }
  }

  std::vector<double> result_zone_vals(zone_ind.size());

  // unique zones' codes (zone values)
  std::set<double> zone_codes(zone_z.begin(),
                              zone_z.end());
  // iterate through zones
  for (double zone_code : zone_codes){
    std::vector<double> vals_in_zone;
    // iterate through all zone layers hex indexes
    // to find those in the zone
    for (int i = 0; i < zone_ind.size(); i++){
      // if the code is of needed zone
      if (zone_z[i] == zone_code){
        // find the same index in raster layer
        for (int j = 0; j < rast_ind.size(); j++){
          if (rast_ind[j] == zone_ind[i]){
            vals_in_zone.push_back(rast_z[j]);
            break;
          }
        }
      }
    }
    // result of a statistical operation on the vector
    double zone_result = 0.0;
    if (stat_type == "min"){
      zone_result = *min_element(std::begin(vals_in_zone),
                                 std::end(vals_in_zone));
    } else if (stat_type == "max"){
      zone_result = *max_element(std::begin(vals_in_zone),
                                 std::end(vals_in_zone));
    } else if (stat_type == "avg"){
      zone_result = vector_average(vals_in_zone);
    } else if (stat_type == "sum"){
      for (auto val : vals_in_zone){
        if (!std::isnan(val)){
          zone_result += val;
        }}
    } else if (stat_type == "majority"){
      zone_result = most_frequent_element(vals_in_zone);
    } else if (stat_type == "minority"){
      zone_result = least_frequent_element(vals_in_zone);
    } else {
      zone_result = -0.0;
    }
    // write the result into the zone's hexs
    for (int i = 0; i < zone_ind.size(); i++){
      if (zone_z[i] == zone_code){
        result_zone_vals[i] =  zone_result;
      }
    }
  }
  // final dict with hex indexes, zones codes and statistic values
  std::map <std::string, std::vector<double>> zonalstats;
  for (int i = 0; i < zone_ind.size(); i++){
    std::vector<double> zone_result;
    zone_result.push_back(zone_z[i]);
    zone_result.push_back(result_zone_vals[i]);
    zonalstats[zone_ind[i]] = zone_result;
  }

  if ((resample_status == 2 && resample_zone) ||
      (resample_status == 1) ||
      (resample_status == 0)){
    return zonalstats;
  } else if (resample_status == 2 && !resample_zone){
    std::map <std::string, std::vector<double>> old_zone;
    for (auto const & pair : zonalstats){
      std::string zone_cell = H3_to_parent(pair.first, h3_level_zone);
      old_zone[zone_cell] = pair.second;
  }
    return old_zone;
  }
}


// Returns cell neighbors of defined order and less with the cell itself

// [[Rcpp::export]]
std::vector<std::string> cell_vecinity(std::string h3s, int radius) {
    H3Index h3 = stringToH3(h3s.std::string::c_str());
    int n = maxKringSize(radius);
    H3Index* out = new H3Index[n]();
    kRing(h3, radius, out);
    int counter = 0;
    for (int i = 0; i < n; ++i) {
      if (out[i] != 0) {
        ++counter;
      }
    }

    std::vector<std::string> v(counter);
    for(int i = 0; i < counter; ++i) {
      char h3s[17];
      h3ToString(out[i], h3s, sizeof(h3s));
      v[i] = h3s;
    }

    // free(out);
    delete[] out;
    return v;
  }


// Calculates base statistics in focal window

// [[Rcpp::export]]
std::map <std::string, double> simple_focal(std::vector<std::string> & inds,
                                            std::vector<double> & z,
                                            const std::string stat_type){
  int h3_level = h3GetResolution(stringToH3(inds[0].std::string::c_str()));
  try{
    if (inds.size() != z.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  int n = inds.size();
  std::map <std::string, double> result;
  for (int i = 0; i < n; i++){
    std::string this_ind = inds[i]; // current cell
    std::vector<std::string> this_vecinity = cell_vecinity(this_ind, h3_level); // current cell neighbors
    // get values in focal window
    std::vector<double> initial_vals;
    for (int j = 0; j < n; j++){
      for (auto const & vec_ind : this_vecinity){
        if (inds[j] == vec_ind){
          initial_vals.push_back(z[j]);
        }
      }
    }
    // calculate statistic in focal window

    double w_result = 0.0;
    if (stat_type == "min"){
      w_result = *min_element(std::begin(initial_vals),
                                 std::end(initial_vals));
    } else if (stat_type == "max"){
      w_result = *max_element(std::begin(initial_vals),
                                 std::end(initial_vals));
    } else if (stat_type == "avg"){
      w_result = vector_average(initial_vals);
    } else if (stat_type == "sum"){
      for (auto val : initial_vals){
        if (!std::isnan(val)){
          w_result += val;
        }}
    } else if (stat_type == "majority"){
      w_result = most_frequent_element(initial_vals);
    } else if (stat_type == "minority"){
      w_result = least_frequent_element(initial_vals);
    } else {
      w_result = -0.0;
    }
    // write the result into hexs
    result[inds[i]] = w_result;
  }
  return result;
  }




// Calculates hex's centers coordinates

// [[Rcpp::export]]
std::map <std::string, std::vector<double>> indexes_to_coords(std::vector<std::string> & inds){
  std::map <std::string, std::vector<double>> geotab;
  int numHexagons = inds.size();
  for (int i = 0; i < numHexagons; i++) {
      GeoCoord hex_center;
      std::vector<double> centroid_coords;
      H3Index h3 = stringToH3(inds[i].std::string::c_str());
      h3ToGeo(h3, &hex_center);
      centroid_coords.push_back(radsToDegs(hex_center.lon));
      centroid_coords.push_back(radsToDegs(hex_center.lat));
      geotab[inds[i]] = centroid_coords;
  }
  return geotab;
}













