#include "cpu_nodynmem.hpp"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <optional>
#include <regex>
#include <string>
#include <algorithm>
#include <unordered_map>

#include <cmath>
#include <queue>

#include "types.hpp"
int hexStringToInt(const std::string &hexString)
{
    int result;
    std::stringstream ss;
    ss << std::hex << hexString;
    ss >> result;
    return result;
}

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        // show error msg - like usage program.exe filename
        std::cout << "ERROR: missing relative filepath" << std::endl;
        return 1;
    }
    // Path to the input text file
    std::string filename = argv[1];

    // Open the file
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;
    }

    // Define your regex pattern (for example: match lines with "error")
    std::regex pattern("0\\ ([A-Fa-f0-9]{3})/([A-Fa-f0-9]{3})", std::regex_constants::icase); // case-insensitive
    std::string line;
    int line_num = 0;

    // Read and process each line

    std::vector<Point> event_points;
    while (std::getline(file, line))
    {
        ++line_num;
        std::smatch match;
        if (std::regex_match(line, match, pattern))
        {
            // std::cout << "Line " << line_num << ": " << hexStringToInt(match[1]) <<" "<< hexStringToInt(match[2]) << std::endl;
            Point p = {(double)hexStringToInt(match[2]), (double)hexStringToInt(match[1])};
            event_points.push_back(p);
        }
    }

    auto max_id = dbscan_algo_fixed_mem(event_points.data(), 1, 1, event_points.size());
    std::cout << max_id << std::endl;

    // std::vector<cluster> clusters;
    // clusters.resize(max_id);
    // cluster current_cluster;
    // current_cluster.cluster_id = 0;

    // for (auto point : event_points)
    // {
    //     clusters[point.cluster_id].members.push_back(point);
    //     clusters[point.cluster_id].cluster_id = point.cluster_id;
        
    //     auto& cluster = clusters[point.cluster_id];
    //     if (cluster.bbox_max_x == -1 || point.x > cluster.bbox_max_x)
    //         cluster.bbox_max_x = point.x;
    //     if (cluster.bbox_max_y == -1 || point.y > cluster.bbox_max_y)
    //         cluster.bbox_max_y = point.y;
    //     if (cluster.bbox_min_x == -1 || point.x < cluster.bbox_min_x)
    //         cluster.bbox_min_x = point.x;
    //     if (cluster.bbox_min_y == -1 || point.y < cluster.bbox_min_y)
    //         cluster.bbox_min_y = point.y;
    // }

    // for (auto cluster : clusters)
    // {
    //     std::cout << cluster.bbox_min_y << "," << cluster.bbox_min_x << "," << cluster.bbox_max_y << "," << cluster.bbox_max_x << std::endl;
    // }

    file.close();
    return 0;
}