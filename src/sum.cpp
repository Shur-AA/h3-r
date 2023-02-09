//sum.cpp
#include <Rcpp.h>
#include "h3api.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <math.h>


using namespace Rcpp;


// internal functions


// Dictionary of base cell's azimuth

const std::map<std::string, int> bc_azimuth = {{"0", 218.290385106721},
                                               {"1", 342.583876444099},
                                               {"2", 339.565502292778},
                                               {"3", 84.3385298378604},
                                               {"4", 309.370561217946},
                                               {"5", 276.5900234123},
                                               {"6", 30.8781696149839},
                                               {"7", 49.4270946789126},
                                               {"8", 302.412537716117},
                                               {"9", 10.309674146324},
                                               {"10", 313.713538667284},
                                               {"11", 351.262408773729},
                                               {"12", 38.6755104962345},
                                               {"13", 16.203220736166},
                                               {"14", 51.8304080271509},
                                               {"15", 351.418122951382},
                                               {"16", 324.654296125929},
                                               {"17", 183.619278114153},
                                               {"18", 347.400874108958},
                                               {"19", 31.5197544157351},
                                               {"20", 194.292400567454},
                                               {"21", 53.6359872457783},
                                               {"22", 305.92918709196},
                                               {"23", 156.106574124094},
                                               {"24", 12.1123540059164},
                                               {"25", 179.562675767721},
                                               {"26", 26.6518863318642},
                                               {"27", 11.1114107906054},
                                               {"28", 342.054997414918},
                                               {"29", 11.171294409814},
                                               {"30", 332.477144661756},
                                               {"31", 352.642886995835},
                                               {"32", 156.342216403357},
                                               {"33", 316.16716182498},
                                               {"34", 210.123976043561},
                                               {"35", 3.49775618632267},
                                               {"36", 201.943635115208},
                                               {"37", 339.144602932168},
                                               {"38", 64.9026064632485},
                                               {"39", 160.943040770228},
                                               {"40", 192.008848045182},
                                               {"41", 2.03374355122804},
                                               {"42", 32.5983362216963},
                                               {"43", 21.8281266741111},
                                               {"44", 344.122232295351},
                                               {"45", 172.289232907796},
                                               {"46", 8.69772604604191},
                                               {"47", 200.987258088485},
                                               {"48", 152.053073103989},
                                               {"49", 16.0766806573519},
                                               {"50", 154.286266938449},
                                               {"51", 26.2878483356397},
                                               {"52", 339.664766874379},
                                               {"53", 352.791355525612},
                                               {"54", 28.1169255496121},
                                               {"55", 200.509893355518},
                                               {"56", 2.54477171830757},
                                               {"57", 340.264280022697},
                                               {"58", 43.8339582496727},
                                               {"59", 342.67812404629},
                                               {"60", 11.7546917007767},
                                               {"61", 168.74765842802},
                                               {"62", 197.866305499552},
                                               {"63", 208.932797684863},
                                               {"64", 200.701445646274},
                                               {"65", 177.167110342362},
                                               {"66", 339.725522359102},
                                               {"67", 152.83178440706},
                                               {"68", 188.364439173619},
                                               {"69", 201.275860528836},
                                               {"70", 154.992950055343},
                                               {"71", 25.7049920546076},
                                               {"72", 237.314462121582},
                                               {"73", 27.8526255167888},
                                               {"74", 338.916723724987},
                                               {"75", 172.831789016884},
                                               {"76", 7.52112264225278},
                                               {"77", 195.273860935452},
                                               {"78", 160.001250135935},
                                               {"79", 148.967469115864},
                                               {"80", 179.511950449523},
                                               {"81", 347.493721316135},
                                               {"82", 18.5181979360588},
                                               {"83", 189.546051224453},
                                               {"84", 203.489398519783},
                                               {"85", 337.120605701631},
                                               {"86", 177.462866751602},
                                               {"87", 328.943684095856},
                                               {"88", 226.006616665999},
                                               {"89", 22.7337106549417},
                                               {"90", 190.184780622089},
                                               {"91", 210.807506107484},
                                               {"92", 171.242030503069},
                                               {"93", 199.934580706627},
                                               {"94", 172.4823059236},
                                               {"95", 157.203247377055},
                                               {"96", 8.33793896205384},
                                               {"97", 244.235133490634},
                                               {"98", 22.2006219025585},
                                               {"99", 236.745771728884},
                                               {"100", 129.376345184138},
                                               {"101", 344.089509437182},
                                               {"102", 154.191315929153},
                                               {"103", 196.938994225748},
                                               {"104", 354.47578585963},
                                               {"105", 221.087675552151},
                                               {"106", 194.939837836892},
                                               {"107", 207.994313007577},
                                               {"108", 169.135403579493},
                                               {"109", 146.515544789238},
                                               {"110", 197.672098858928},
                                               {"111", 231.5110548796},
                                               {"112", 177.628842629001},
                                               {"113", 242.195406460574},
                                               {"114", 138.903570883011},
                                               {"115", 161.433226551426},
                                               {"116", 260.50009696699},
                                               {"117", 87.0886164367714},
                                               {"118", 41.5910731408811},
                                               {"119", 221.70146541449},
                                               {"120", 173.054981967017},
                                               {"121", 314.358924253677}};

// Dictionary of base cell's hemisphere (where is 'up' for them)
// 0 - northern hemisphere
// 1 - southern hemisphere
// 2 - pentagon

const std::map<std::string, int> bc_hms = {{"0", 0},
                                          {"1", 0},
                                          {"2", 0},
                                          {"3", 0},
                                          {"4", 2},
                                          {"5", 0},
                                          {"6", 0},
                                          {"7", 0},
                                          {"8", 0},
                                          {"9", 0},
                                          {"10", 0},
                                          {"11", 0},
                                          {"12", 0},
                                          {"13", 0},
                                          {"14", 2},
                                          {"15", 0},
                                          {"16", 0},
                                          {"17", 1},
                                          {"18", 0},
                                          {"19", 0},
                                          {"20", 1},
                                          {"21", 0},
                                          {"22", 0},
                                          {"23", 1},
                                          {"24", 2},
                                          {"25", 1},
                                          {"26", 0},
                                          {"27", 0},
                                          {"28", 0},
                                          {"29", 0},
                                          {"30", 0},
                                          {"31", 0},
                                          {"32", 1},
                                          {"33", 0},
                                          {"34", 1},
                                          {"35", 0},
                                          {"36", 1},
                                          {"37", 0},
                                          {"38", 2},
                                          {"39", 1},
                                          {"40", 1},
                                          {"41", 0},
                                          {"42", 0},
                                          {"43", 0},
                                          {"44", 0},
                                          {"45", 1},
                                          {"46", 0},
                                          {"47", 1},
                                          {"48", 1},
                                          {"49", 2},
                                          {"50", 1},
                                          {"51", 0},
                                          {"52", 0},
                                          {"53", 0},
                                          {"54", 0},
                                          {"55", 1},
                                          {"56", 0},
                                          {"57", 0},
                                          {"58", 2},
                                          {"59", 0},
                                          {"60", 0},
                                          {"61", 1},
                                          {"62", 1},
                                          {"63", 2},
                                          {"64", 1},
                                          {"65", 1},
                                          {"66", 0},
                                          {"67", 1},
                                          {"68", 1},
                                          {"69", 1},
                                          {"70", 1},
                                          {"71", 0},
                                          {"72", 2},
                                          {"73", 0},
                                          {"74", 0},
                                          {"75", 1},
                                          {"76", 0},
                                          {"77", 1},
                                          {"78", 1},
                                          {"79", 1},
                                          {"80", 1},
                                          {"81", 0},
                                          {"82", 0},
                                          {"83", 2},
                                          {"84", 1},
                                          {"85", 0},
                                          {"86", 1},
                                          {"87", 0},
                                          {"88", 1},
                                          {"89", 0},
                                          {"90", 1},
                                          {"91", 1},
                                          {"92", 1},
                                          {"93", 1},
                                          {"94", 1},
                                          {"95", 1},
                                          {"96", 0},
                                          {"97", 2},
                                          {"98", 0},
                                          {"99", 1},
                                          {"100", 1},
                                          {"101", 0},
                                          {"102", 1},
                                          {"103", 1},
                                          {"104", 0},
                                          {"105", 1},
                                          {"106", 1},
                                          {"107", 2},
                                          {"108", 1},
                                          {"109", 1},
                                          {"110", 1},
                                          {"111", 1},
                                          {"112", 1},
                                          {"113", 1},
                                          {"114", 1},
                                          {"115", 1},
                                          {"116", 1},
                                          {"117", 2},
                                          {"118", 1},
                                          {"119", 1},
                                          {"120", 1},
                                          {"121", 1}};


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

// добавить присвоение центрального, если все разные
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


// which part of an array has the same value as given
double how_many_with_thesame_val(std::vector<double> &arr,
                                 double that_val)
{
  int n = arr.size();
  int dublicate_counter = 0;
  for (int i = 0; i < n; ++i){
    if (arr[i] == that_val){
      ++dublicate_counter;
    }
  if (dublicate_counter > 0){
    --dublicate_counter; // as kring vicinity includes center cell
    --n;
  }
  }
  return dublicate_counter / n;
}



// which part of an array has less values than given
double how_many_with_less_val(std::vector<double> &arr,
                              double that_val)
{
  int n = arr.size();
  int counter = 0;
  for (int i = 0; i < n; ++i){
    if (arr[i] < that_val){
      ++counter;
    }
    if (counter > 0){
      --n;
    }
  }
  return (counter * 1.0) / (n * 1.0);
}


// which part of an array has greater values than given
double how_many_with_greater_val(std::vector<double> &arr,
                              double that_val)
{
  int n = arr.size();
  int counter = 0;
  for (int i = 0; i < n; ++i){
    if (arr[i] > that_val){
      ++counter;
    }
    if (counter > 0){
      --n;
    }
  }
  return (counter * 1.0) / (n * 1.0);
}


// how many distinct values are in the array
double unique_elements_num(std::vector<double> &arr)
{
  std::set<double> s(arr.begin(), arr.end());
  return s.size();
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


// calculates 3x3 matrices determinant
double get_33_determinant(double a11, double a12, double a13,
                          double a21, double a22, double a23,
                          double a31, double a32, double a33){
  return a11 * a22 * a33 - a11 * a23 * a32 - a12 * a21 * a33 +
          a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31;
}


// calculates 2x2 matrices determinant
double get_22_determinant(double a11, double a12,
                          double a21, double a22){
  return a11 * a22 - a12 * a21;
}



// finds coefficients of the plane equation which approximates several points (N > 3)
// https://pikabu.ru/story/postroenie_ploskosti_oblaku_tochek_metodom_naimenshikh_kvadratov_8125372

std::vector<double> get_plane_coefs(std::vector<double> & xx,
                                    std::vector<double> & yy,
                                    std::vector<double> & zz){
  try{
    if (xx.size() != yy.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  try{
    if (xx.size() != zz.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }

  int n = xx.size();

  double sumX = 0;
  double sumY = 0;
  double sumZ = 0;
  double sumX2 = 0;
  double sumY2 = 0;
  double sumZ2 = 0;
  double sumXY = 0;
  double sumYZ = 0;
  double sumZX = 0;
  for (int i; i < n; ++i){
    sumX += xx[i];
    sumY += yy[i];
    sumZ += zz[i];
    sumX2 += (xx[i] * xx[i]);
    sumY2 += (yy[i] * yy[i]);
    sumZ2 += (zz[i] * zz[i]);
    sumXY += (xx[i] * yy[i]);
    sumYZ += (yy[i] * zz[i]);
    sumZX += (zz[i] * xx[i]);
  }
  double x_x = sumX2 - (sumX * sumX) / n;
  double y_y = sumY2 - (sumY * sumY) / n;
  double z_z = sumZ2 - (sumZ * sumZ) / n;
  double x_y = sumXY - (sumX * sumY) / n;
  double y_z = sumYZ - (sumY * sumZ) / n;
  double z_x = sumZX - (sumZ * sumX) / n;

  double detA1 = get_33_determinant(1, x_y, z_x, 1, y_y, y_z, 1, y_z, z_z);
  double detB1 = get_33_determinant(x_x, 1, z_x, x_y, 1, y_z, z_x, 1, z_z);
  double detC1 = get_33_determinant(x_x, x_y, 1, x_y, y_y, 1, z_x, y_z, 1);

  double l1 = sqrt(detA1 * detA1 + detB1 * detB1 + detC1 * detC1);

  double detA2 = get_33_determinant(-1, x_y, z_x, 1, y_y, y_z, 1, y_z, z_z);
  double detB2 = get_33_determinant(x_x, -1, z_x, x_y, 1, y_z, z_x, 1, z_z);
  double detC2 = get_33_determinant(x_x, x_y, -1, x_y, y_y, 1, z_x, y_z, 1);

  double l2 = sqrt(detA2 * detA2 + detB2 * detB2 + detC2 * detC2);

  double detA3 = get_33_determinant(1, x_y, z_x, 1, y_y, y_z, -1, y_z, z_z);
  double detB3 = get_33_determinant(x_x, 1, z_x, x_y, 1, y_z, z_x, -1, z_z);
  double detC3 = get_33_determinant(x_x, x_y, 1, x_y, y_y, 1, z_x, y_z, -1);

  double l3 = sqrt(detA3 * detA3 + detB3 * detB3 + detC3 * detC3);

  double a = 0;
  double b = 0;
  double c = 0;
  if (l1 >= l2 && l1 >= l3){
    a = detA1 / l1;
    b = detB1 / l1;
    c = detC1 / l1;
  } else if (l2 >= l1 && l2 >= l3){
    a = detA2 / l2;
    b = detB2 / l2;
    c = detC2 / l2;
  } else if (l3 >= l1 && l3 >= l2){
    a = detA3 / l3;
    b = detB3 / l3;
    c = detC3 / l3;
  }
  double d = -(a * sumX + b * sumY + c * sumZ) / n;

  std::vector<double> eq_coefs{a, b, c, d};

  return eq_coefs;
}



// finds coefficients of the plane equation which approximates several points (N > 3)
// https://medium.com/swlh/how-to-find-the-least-squares-plane-from-a-cloud-of-point-using-excel-numbers-etc-92b66f852522

std::vector<double> get_plane_coefs2(std::vector<double> & xx,
                                    std::vector<double> & yy,
                                    std::vector<double> & zz){
  try{
    if (xx.size() != yy.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  try{
    if (xx.size() != zz.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }

  int n = xx.size();

  double sumX = 0;
  double sumY = 0;
  double sumZ = 0;
  double sumX2 = 0;
  double sumY2 = 0;
  double sumXY = 0;
  double sumYZ = 0;
  double sumZX = 0;
  for (int i; i < n; ++i){
    sumX += xx[i];
    sumY += yy[i];
    sumZ += zz[i];
    sumX2 += (xx[i] * xx[i]);
    sumY2 += (yy[i] * yy[i]);
    sumXY += (xx[i] * yy[i]);
    sumYZ += (yy[i] * zz[i]);
    sumZX += (zz[i] * xx[i]);
  }

  double detA = get_33_determinant(sumX2, sumXY, sumX,
                                   sumXY, sumY2, sumY,
                                    sumX,  sumY,  n);
  if (detA != 0){
    detA = 1 / detA;
    // finding minor matrix
    double m11 = get_22_determinant(sumY2, sumY, sumY, n);
    double m12 = get_22_determinant(sumXY, sumY, sumX, n);
    double m13 = get_22_determinant(sumXY, sumY2, sumX, sumY);
    double m21 = get_22_determinant(sumXY, sumX, sumY, n);
    double m22 = get_22_determinant(sumX2, sumX, sumX, n);
    double m23 = get_22_determinant(sumX2, sumXY, sumX, sumY);
    double m31 = get_22_determinant(sumXY, sumX, sumY2, sumY);
    double m32 = get_22_determinant(sumX2, sumX, sumXY, sumY);
    double m33 = get_22_determinant(sumX2, sumXY, sumXY, sumY2);
    // finding adjugate matrix
    m12 = -m12;
    m21 = -m21;
    m23 = -m23;
    m32 = -m32;
    // finding transposed adjugate matrix
    double temp_station;
    temp_station = m12;
    m12 = m21;
    m21 = temp_station;
    temp_station = m13;
    m13 = m31;
    m31 = temp_station;
    temp_station = m23;
    m23 = m32;
    m32 = temp_station;
    // final inverted matrix
    m11 = m11 * detA;
    m12 = m12 * detA;
    m13 = m13 * detA;
    m21 = m21 * detA;
    m22 = m22 * detA;
    m23 = m23 * detA;
    m31 = m31 * detA;
    m32 = m32 * detA;
    m33 = m33 * detA;
    // coefs
    double a = m11 * sumZX + m12 * sumYZ + m13 * sumZ;
    double b = m21 * sumZX + m22 * sumYZ + m23 * sumZ;
    double c = m31 * sumZX + m32 * sumYZ + m33 * sumZ;
    double d = -(a * sumX + b * sumY + c * sumZ) / n;

    std::vector<double> eq_coefs{a, b, c, d};
    return eq_coefs;

  } else {
    std::vector<double> eq_coefs{0, 0, 0, 0};
    return eq_coefs;
  }
}


// angle between plane and horizontal plane
double get_slope(double aa, double bb, double cc){
  double cosA = abs(cc) / sqrt(aa * aa + bb * bb + cc * cc);
  double angle = std::acos(cosA) * 180.0 / 3.14159265358979323846;
  return angle;
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


// given children H3 indexes returns
// list of parents in the closest level

// [[Rcpp::export]]
std::map <std::string, std::string> get_direct_parents(const std::vector<std::string> & children_ind){
  int level_from = h3GetResolution(stringToH3(children_ind[0].std::string::c_str()));
  // parents of every child
  std::map <std::string, std::string> chp;
  if (level_from > 0){
    // loop through children string indexes
    for (int i = 0; i < children_ind.size(); i++) {
      H3Index h3 = stringToH3(children_ind[i].std::string::c_str());
      H3Index h3Parent = h3ToParent(h3, level_from - 1);
      char h3ParentStr[17];
      h3ToString(h3Parent, h3ParentStr, sizeof(h3ParentStr));
      chp[children_ind[i]] = h3ParentStr;
    }
  }
  return chp;
}


// given children H3 indexes and corresponding values defines
// list of parents in the closest level and aggregates values

// [[Rcpp::export]]
std::map <std::string, double> resample_up(const std::string func,
                                           const std::vector<std::string> & children_ind,
                                           const std::vector<double> & children_vals){
  int level_from = h3GetResolution(stringToH3(children_ind[0].std::string::c_str()));

  try{
    if (children_ind.size() != children_vals.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal lengths exception - unpredictable result" << std::endl;
  }
  std::map <std::string, double> parentab;
  if (level_from > 0){
    // parents of every child
    std::vector<std::string> chp;
    // loop through children string indexes
    for (int i = 0; i < children_ind.size(); i++) {
      H3Index h3 = stringToH3(children_ind[i].std::string::c_str());
      H3Index h3Parent = h3ToParent(h3, level_from - 1);
      char h3ParentStr[17];
      h3ToString(h3Parent, h3ParentStr, sizeof(h3ParentStr));
      chp.push_back(h3ParentStr);
      }

    std::set<std::string> prnts(chp.begin(),
                                chp.end());
    // loop through unique parent string indexes
    for (std::string pind : prnts) {
      H3Index h3 = stringToH3(pind.std::string::c_str());
      int n = maxH3ToChildrenSize(h3, level_from);  // define maximum number of the hex's children
      H3Index* h3Children = new H3Index[n];  // vector for children
      h3ToChildren(h3, level_from, h3Children); // children list in H3 format
      // list of children of every parent
      std::vector<std::string> children_list;
      for (int j = 0; j < n; ++j) {
        char h3Str[17];
        h3ToString(h3Children[j], h3Str, sizeof(h3Str));
        children_list.push_back(h3Str);
      }
      delete[] h3Children;
      // list of children;s values of every parent
      std::vector<double> children_list_vals;
      for (std::string chind : children_list) {
        for (int k = 0; k < children_ind.size(); ++k) {
          if (chind == children_ind[k]){
            children_list_vals.push_back(children_vals[k]);
          }
        }
      }
      // result value of aggregation
      double w_result = 0.0;
      if (func == "max"){
        w_result = *max_element(std::begin(children_list_vals),
                                std::end(children_list_vals));
      } else if (func == "avg"){
        w_result = vector_average(children_list_vals);
      } else if (func == "sum"){
        for (auto val : children_list_vals){
          if (!std::isnan(val)){
            w_result += val;
          }}
      } else if (func == "majority"){
        w_result = most_frequent_element(children_list_vals);
      } else {
        w_result = -0.0;
      }
      // write the result into parent hexs
      parentab[pind] = w_result;
    }
  }
  return parentab;
}


// given children H3 indexes and corresponding values defines
// list of parents in any coarser level and aggregates values
// this function is overloading

// [[Rcpp::export]]
std::map <std::string, double> resample_up_any(const std::string func,
                                           const int level_to,
                                           const std::vector<std::string> & children_ind,
                                           const std::vector<double> & children_vals){
  int level_from = h3GetResolution(stringToH3(children_ind[0].std::string::c_str()));

  try{
    if (children_ind.size() != children_vals.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal lengths exception - unpredictable result" << std::endl;
  }
  std::map <std::string, double> parentab;
  if (level_from > 0 && level_from > level_to){
    // parents of every child
    std::vector<std::string> chp;
    // loop through children string indexes
    for (int i = 0; i < children_ind.size(); i++) {
      H3Index h3 = stringToH3(children_ind[i].std::string::c_str());
      H3Index h3Parent = h3ToParent(h3, level_to);
      char h3ParentStr[17];
      h3ToString(h3Parent, h3ParentStr, sizeof(h3ParentStr));
      chp.push_back(h3ParentStr);
    }

    std::set<std::string> prnts(chp.begin(),
                                chp.end());
    // loop through unique parent string indexes
    for (std::string pind : prnts) {
      H3Index h3 = stringToH3(pind.std::string::c_str());
      int n = maxH3ToChildrenSize(h3, level_from);  // define maximum number of the hex's children
      H3Index* h3Children = new H3Index[n];  // vector for children
      h3ToChildren(h3, level_from, h3Children); // children list in H3 format
      // list of children of every parent
      std::vector<std::string> children_list;
      for (int j = 0; j < n; ++j) {
        char h3Str[17];
        h3ToString(h3Children[j], h3Str, sizeof(h3Str));
        children_list.push_back(h3Str);
      }
      delete[] h3Children;
      // list of children's values of every parent
      std::vector<double> children_list_vals;
      for (std::string chind : children_list) {
        for (int k = 0; k < children_ind.size(); ++k) {
          if (chind == children_ind[k]){
            children_list_vals.push_back(children_vals[k]);
          }
        }
      }
      // result value of aggregation
      double w_result = 0.0;
      if (func == "max"){
        w_result = *max_element(std::begin(children_list_vals),
                                std::end(children_list_vals));
      } else if (func == "avg"){
        w_result = vector_average(children_list_vals);
      } else if (func == "sum"){
        for (auto val : children_list_vals){
          if (!std::isnan(val)){
            w_result += val;
          }}
      } else if (func == "majority"){
        w_result = most_frequent_element(children_list_vals);
      } else {
        w_result = -0.0;
      }
      // write the result into parent hexs
      parentab[pind] = w_result;
    }
  }
  return parentab;
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


// Calculates sum, mean, min, max, majority or minority
// of values in H3 cells of the second raster
// which corresponds to H3 cells in the first
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
                                            const std::string stat_type,
                                            int vecinity){
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
    std::vector<std::string> this_vecinity = cell_vecinity(this_ind, vecinity); // current cell neighbors
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
    } else if (stat_type == "procent"){
      double this_val = z[i]; // current central cell value
      w_result = how_many_with_thesame_val(initial_vals, this_val);
    } else if (stat_type == "procentile_less"){
      double this_val = z[i];
      w_result = how_many_with_less_val(initial_vals, this_val);
    } else if (stat_type == "procentile_more"){
      double this_val = z[i];
      w_result = how_many_with_greater_val(initial_vals, this_val);
    } else if (stat_type == "variety"){
      w_result = unique_elements_num(initial_vals);
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



// Calculates azimuth of a cell

// [[Rcpp::export]]
double cell_azimuth(std::string & h3_index){
  int h3_level = h3GetResolution(stringToH3(h3_index.std::string::c_str()));
  H3Index h3 = stringToH3(h3_index.std::string::c_str());
  int bc = h3GetBaseCell(h3);
  double azimuth = bc_azimuth.at(std::to_string(bc));
  if (h3_level % 2 != 0){
    azimuth = azimuth - 19.1066053508691;
  }
  return azimuth;
}


// Calculates gradient and aspect in focal window

// [[Rcpp::export]]
std::map <std::string, std::vector<double>> gradient_aspect(std::vector<std::string> & inds,
                                                           std::vector<double> & z,
                                                           const std::string stat_type){
  std::map <std::string, std::vector<double>> geotab;
  int h3_level = h3GetResolution(stringToH3(inds[0].std::string::c_str()));
  try{
    if (inds.size() != z.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  int n = inds.size();

  for (int i = 0; i < n; i++){
    std::string this_ind = inds[i]; // current cell
    std::vector<std::string> this_vecinity = cell_vecinity(this_ind, 1); // current cell immediate neighbors
    // get values in focal window
    // assume that order in this_vecinity = in initial_vals
    std::vector<double> initial_vals;
    // coords of central cell vertices
    std::vector<double> px1;
    std::vector<double> py1;
    std::vector<double> px2;
    std::vector<double> py2;
    for (auto const & vec_ind : this_vecinity){
      for (int j = 0; j < n; j++){
        if (inds[j] == vec_ind){
          initial_vals.push_back(z[j]);
        }
      }
      if (vec_ind != this_ind){
        // find the edge between center and current neighbor cell
        H3Index h3Origin = H3_EXPORT(stringToH3)(this_ind.std::string::c_str());
        H3Index h3Destination = H3_EXPORT(stringToH3)(vec_ind.std::string::c_str());
        H3Index h3Edge = H3_EXPORT(getH3UnidirectionalEdge)(h3Origin, h3Destination);
        // get coords of the edge vertices
        GeoBoundary geoBoundary;
        H3_EXPORT(getH3UnidirectionalEdgeBoundary)(h3Edge, &geoBoundary);
        px1.push_back(radsToDegs(geoBoundary.verts[0].lon));
        py1.push_back(radsToDegs(geoBoundary.verts[0].lat));
        px2.push_back(radsToDegs(geoBoundary.verts[1].lon));
        py2.push_back(radsToDegs(geoBoundary.verts[1].lat));
      } else{
        // get coords of the current central cell
        uint64_t h3 = stringToH3(this_ind.std::string::c_str());
        GeoCoord geoCoord;
        h3ToGeo(h3, &geoCoord);
        px1.push_back(radsToDegs(geoCoord.lon));
        py1.push_back(radsToDegs(geoCoord.lat));
        px2.push_back(-1);
        py2.push_back(-1);
      }
    }
    // calculating average values for central cell vertices on 3 cells
    int vnum = this_vecinity.size();
    std::vector<double> vert_avg; // array for average values plus central value
    for (int k = 0; k < vnum; ++k){
      if (this_vecinity[k] == this_ind){
        vert_avg.push_back(z[i]);
      } else{
        for (int g = 0; g < vnum; ++g){
          if ((px1[k] == px2[g]) && (py1[k] == py2[g])){
            std::vector<double> triple_vals{z[i], initial_vals[k], initial_vals[g]};
            vert_avg.push_back(vector_average(triple_vals));
          }
        }
      }
    }

    std::vector<double> plane_eq_coefs = get_plane_coefs2(px1, py1, vert_avg);
    plane_eq_coefs.push_back(get_slope(plane_eq_coefs[0], plane_eq_coefs[1], plane_eq_coefs[2]));
  // тестовая выгрузка!
    geotab[this_ind] = plane_eq_coefs;

  }
  return geotab;
}



// Calculates drainage in focal window

// [[Rcpp::export]]
std::map <std::string, std::string> drainage(std::vector<std::string> & inds,
                                             std::vector<double> & z){
  std::map <std::string, std::string> geotab;
  int h3_level = h3GetResolution(stringToH3(inds[0].std::string::c_str()));
  try{
    if (inds.size() != z.size()){
      throw 2; // not equal lengths exception
    }}
  catch(int x){
    std::cout<<"not equal vector lengths exception - unpredictable result" << std::endl;
  }
  int n = inds.size();

  for (int i = 0; i < n; i++){
    std::string this_ind = inds[i]; // current cell
    std::vector<std::string> this_vecinity = cell_vecinity(this_ind, 1); // current cell immediate neighbors
    // get values in focal window
    // assume that order in this_vecinity = in initial_vals
    std::string desc_slope_hex = this_ind;  // h3 index of the cell with descent slope
    double desc_slope_val = 0;  // h difference with the cell with descent slope
    for (auto const & vec_ind : this_vecinity){
      for (int j = 0; j < n; j++){
        if (inds[j] == vec_ind && vec_ind != this_ind && ((z[i] - z[j]) > desc_slope_val)){
          desc_slope_val = z[i] - z[j];
          desc_slope_hex = vec_ind;
        }
      }
    }
    geotab[this_ind] = desc_slope_hex;

  }
  return geotab;
}







