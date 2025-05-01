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


double distance(const Point &a, const Point &b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) +
                     (a.y - b.y) * (a.y - b.y));
}

struct cluster_member
{
    int zig_zag_index = -1;
    uint16_t row;
    uint16_t col;
};



struct point
{
    int x = -1;
    int y = -1;
};

int hexStringToInt(const std::string &hexString)
{
    int result;
    std::stringstream ss;
    ss << std::hex << hexString;
    ss >> result;
    return result;
}

std::vector<int> region_query(const std::vector<Point> &points, int idx, double eps)
{
    std::vector<int> neighbors;
    for (int i = 0; i < points.size(); ++i)
    {
        if (distance(points[idx], points[i]) <= eps)
        {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

void expand_cluster(std::vector<Point> &points, int idx, int clusterID, double eps, int minPts, std::vector<cluster> &clusters)
{
    std::queue<int> q;

    // get the neighbors around the point we are indexed at currently
    auto neighbors = region_query(points, idx, eps);

    if (neighbors.size() < minPts)
    {
        points[idx].label = PointLabel::NOISE;
        return;
    }

    // if the point has neighbors set its cluster id and set its status as clustered
    points[idx].cluster_id = clusterID;
    points[idx].label = PointLabel::CLUSTERED;

    // loop through all of the neighbors of the specific point being clustered and 
    // see if they have been categorized yet. if they havent, add them to the cluster
    // and then see if it has neighbors as well
    for (int n : neighbors)
    { // add all the neighbors to the cluster id
        if (points[n].label == PointLabel::UNVISITED)
        {
            points[n].label = PointLabel::CLUSTERED;
            points[n].cluster_id = clusterID;

            // get the neighbors of these neighbors as well
            auto subNeighbors = region_query(points, n, eps);
            if (subNeighbors.size() >= minPts)
            {
                for (int sn : subNeighbors)
                {
                    q.push(sn);
                }
            }
        }
        else if (points[n].label == PointLabel::NOISE)
        {
            points[n].label = PointLabel::CLUSTERED;
            points[n].cluster_id = clusterID;
        }
    }

    while (!q.empty())
    {
        int current = q.front();
        q.pop();

        auto currentNeighbors = region_query(points, current, eps);
        if (currentNeighbors.size() >= minPts)
        {
            for (int n : currentNeighbors)
            {
                if (points[n].label == PointLabel::UNVISITED)
                {
                    points[n].label = PointLabel::CLUSTERED;
                    points[n].cluster_id = clusterID;
                    q.push(n);
                }
            }
        }
    }
}

int dbscan_algo(std::vector<Point> &points, double eps, double min_points)
{
    int cluster_id = 0;
    std::vector<cluster> clusters;
    for (int i = 0; i < points.size(); i++)
    {
        if (points[i].label == PointLabel::UNVISITED)
        {
            expand_cluster(points, i, cluster_id, eps, min_points, clusters);
            if (points[i].cluster_id == cluster_id)
            {
                ++cluster_id;
            }
        }
    }
    return cluster_id;
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

    auto max_id = dbscan_algo(event_points, 1, 1);
    std::cout << max_id << std::endl;

    std::vector<cluster> clusters;
    clusters.resize(max_id);
    cluster current_cluster;
    current_cluster.cluster_id = 0;

    for (auto point : event_points)
    {
        clusters[point.cluster_id].members.push_back(point);
        clusters[point.cluster_id].cluster_id = point.cluster_id;
        
        auto& cluster = clusters[point.cluster_id];
        if (cluster.bbox_max_x == -1 || point.x > cluster.bbox_max_x)
            cluster.bbox_max_x = point.x;
        if (cluster.bbox_max_y == -1 || point.y > cluster.bbox_max_y)
            cluster.bbox_max_y = point.y;
        if (cluster.bbox_min_x == -1 || point.x < cluster.bbox_min_x)
            cluster.bbox_min_x = point.x;
        if (cluster.bbox_min_y == -1 || point.y < cluster.bbox_min_y)
            cluster.bbox_min_y = point.y;
    }

    for (auto cluster : clusters)
    {
        std::cout << cluster.bbox_min_y << "," << cluster.bbox_min_x << "," << cluster.bbox_max_y << "," << cluster.bbox_max_x << std::endl;
    }

    file.close();
    return 0;
}