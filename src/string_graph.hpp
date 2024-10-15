#ifndef STRING_GRAPH_HPP
#define STRING_GRAPH_HPP



#include <unordered_map>
#include <string>
#include <algorithm>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <tuple>
#include <set>
#include <unordered_set>
#include "h3api.h"
#include <map>


namespace stringgraph {

    typedef std::tuple<double, double> dpoint;
    typedef std::vector<dpoint> hexahonal_vec;

    class Graph {
        public:
            void make_graph(std::unordered_map<std::string, std::string> h3_data);
            std::unordered_map<std::string, std::string> find_path(std::string cell);
        private:
            std::unordered_set<std::string> all_indexes;
            std::unordered_map<std::string, hexahonal_vec> umap_view;
            std::unordered_map<std::string, std::unordered_map<std::string, int>> weighted_graph_data;
            
            void convert_h3(std::unordered_map<std::string, std::string> h3_data);
            void make_graph_data(std::unordered_map<std::string, std::string> h3_data);
    };
}

#endif