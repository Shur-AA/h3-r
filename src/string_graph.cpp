#include <iostream>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <tuple>
#include <fstream>
#include <set>
#include <unordered_set>
#include <type_traits>
#include <queue>
#include <bitset>
#include <iostream>
#include <sstream>
#include "h3api.h"
#include <map>


#include "string_graph.hpp"


std::vector<std::string> cell_vecinity1(std::string h3s, int radius) {
  H3Index h3 = stringToH3(h3s.std::string::c_str());
  int n = maxKringSize(radius) - 1; // exclude center cell
  H3Index* out = new H3Index[n]();
  hexRing(h3, radius, out);
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




void stringgraph::Graph::convert_h3(std::unordered_map<std::string, std::string> h3_data){

    auto dpoint_ring = [](const std::string cell){
        std::string mod_string = cell;
        mod_string.erase(std::remove(mod_string.begin(), mod_string.end(), '"'), mod_string.end());
        unsigned long long result = strtoull(mod_string.data(), NULL, 16);
        GeoBoundary boundary;
        H3Index indexed = result;
        hexahonal_vec hvector;
        h3ToGeoBoundary(indexed, &boundary);
        for (int v = 0; v < boundary.numVerts; v++){
            double latit = radsToDegs(boundary.verts[v].lat);
            double lngit = radsToDegs(boundary.verts[v].lon);
            dpoint lonlat = std::make_tuple(lngit, latit);
            hvector.push_back(lonlat);
        }
        return hvector;
    };


    for (const auto& pair : h3_data) {
      const auto& cell = pair.first;
      const auto& group_id = pair.second;
      hexahonal_vec cell_vector = dpoint_ring(cell);
      umap_view.emplace(cell, cell_vector);
    }
}


void stringgraph::Graph::make_graph_data(std::unordered_map<std::string, std::string> h3_data){

    std::string groupid;
    std::vector<std::string> zd_group;
    for (const auto& pair : h3_data) {
      const auto& cell = pair.first;
      const auto& group_id = pair.second;
      groupid = pair.second;
      all_indexes.insert(cell);
      all_indexes.insert(group_id);
    }

    for (const auto& pair : umap_view) {
      const auto& cellname = pair.first;
      const auto& hexvec = pair.second;


      std::unordered_map<std::string, int> inner_weight;

      for (const auto& pair2 : umap_view) {
        const auto& cellname2 = pair2.first;
        const auto& hexvec2 = pair2.second;

        if ((cellname != cellname2) && (h3_data.at(cellname) == h3_data.at(cellname2))) {
          if (std::find_first_of(hexvec.begin(), hexvec.end(), hexvec2.begin(), hexvec2.end()) != hexvec.end()) {
            inner_weight.emplace(cellname2, 1);
          }
        }
      }
      weighted_graph_data.emplace(cellname, inner_weight);
      zd_group.push_back(cellname);
    }


    // defining neighbors of the aim cell

      std::vector<std::string> this_vicinity = cell_vecinity1(groupid, 1);
      std::unordered_map<std::string, int> inner_weight;
      for (auto const & vic_ind : this_vicinity){
        if (std::find(zd_group.begin(), zd_group.end(), vic_ind) != zd_group.end()){
          inner_weight.emplace(vic_ind, 1);
        }
      }
      weighted_graph_data.emplace(groupid, inner_weight);





}


void stringgraph::Graph::make_graph(std::unordered_map<std::string, std::string> h3_data){
    convert_h3(h3_data);
    make_graph_data(h3_data);
}


std::unordered_map<std::string, std::string> stringgraph::Graph::find_path(std::string cell){

    std::priority_queue<std::pair<int, std::string>, std::vector<std::pair<int, std::string>>> frontier;
    std::pair<int, std::string> start_pair = std::make_pair(0, cell);
    frontier.push(start_pair);

    std::unordered_map<std::string, std::string> came_from;
    std::unordered_map<std::string, int> cost_of_travel;

    if (all_indexes.find(cell) == all_indexes.end()) {
      std::cout<<"Cell is not in all_indexes"<<std::endl;
      return came_from;
    }


    came_from.emplace(cell, std::string("-1"));
    cost_of_travel.emplace(cell, 0);

    if (weighted_graph_data.find(cell) != weighted_graph_data.end()) {
      int total_dist = 0;

      while (!frontier.empty()) {
        std::pair<int, std::string> current_n = frontier.top();
        frontier.pop();

        if (weighted_graph_data.find(current_n.second) != weighted_graph_data.end() &&
            cost_of_travel.find(current_n.second) != cost_of_travel.end()) {

          std::unordered_map<std::string, int> neighbors = weighted_graph_data.at(current_n.second);

          for (auto& pair : neighbors) {
            const auto& neighb_index = pair.first;
            int weight_min = pair.second;

            int new_cost = cost_of_travel.at(current_n.second) + weight_min;
            total_dist = new_cost;

            if (cost_of_travel.find(neighb_index) == cost_of_travel.end()) {
              cost_of_travel.emplace(neighb_index, new_cost);
              std::pair<int, std::string> newpriority = std::make_pair(new_cost, neighb_index);
              frontier.push(newpriority);
              came_from.emplace(neighb_index, current_n.second);
            } else {
                        if (new_cost < cost_of_travel.at(neighb_index)){

                            cost_of_travel.at(neighb_index) = new_cost;
                            std::pair<int, std::string> newpriority = std::make_pair(new_cost, neighb_index);
                            frontier.push(newpriority);
                            came_from.at(neighb_index) = current_n.second;

                        }
                    }
                }
            }
        }
    }

    std::unordered_map<std::string, std::string> reverse_back;

    std::vector<std::string> fromnodes;
    std::transform(came_from.begin(), came_from.end(), std::back_inserter(fromnodes),
                   [](const std::pair<std::string, std::string>& kv) { return kv.first; });

    for (std::string fnode:fromnodes){
        std::string current_node = fnode;
        while (current_node != cell){
            std::string prev_node = came_from.at(current_node);
            reverse_back.emplace(current_node, prev_node);
            current_node = prev_node;
        }
    }

    return reverse_back;
}
