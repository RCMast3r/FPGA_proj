#include "types.hpp"
#include <cmath>
#include <iostream>
double distance(const Point &a, const Point &b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) +
                     (a.y - b.y) * (a.y - b.y));
}
// return the number of neighbors that were detected
size_t region_query_no_dyn(Point points[], int idx, double eps, size_t num_points_to_cluster, size_t neighbors[])
{
    size_t neighbors_size= 0;
    for(int i=0; i<num_points_to_cluster; i++)
    {
        if(distance(points[idx], points[i]) <= eps)
        {
            neighbors[neighbors_size]=i;
            neighbors_size++;
        }
    }
    return neighbors_size;
}

void expand_cluster_no_dyn(Point points[], int idx, int clusterID, double eps, int minPts, cluster clusters[], size_t num_points_to_cluster)
{
    StaticQueue<size_t, MAX_SEARCH_QUEUE_SIZE> q;
    
    // get the neighbors around the point we are indexed at currently
    size_t neighbors[MAX_NUM_NEIGHBORS];
    // for(size_t i =0; i<MAX_NUM_NEIGHBORS; i++){
    //     neighbors[i]=0;
    // }
    auto num_neighbors = region_query_no_dyn(points, idx, eps, num_points_to_cluster, neighbors);

    if (num_neighbors < minPts)
    {
        points[idx].label = PointLabel::NOISE;
        return;
    }

    points[idx].cluster_id = clusterID;
    points[idx].label = PointLabel::CLUSTERED;

    // for (int n : neighbors)
    for(size_t i=0; i<num_neighbors; i++)
    { // add all the neighbors to the cluster id
        auto n = neighbors[i];
        if (points[n].label == PointLabel::UNVISITED)
        {
            points[n].label = PointLabel::CLUSTERED;
            points[n].cluster_id = clusterID;

            // get the neighbors of these neighbors as well
            size_t sub_neighbors[MAX_NUM_NEIGHBORS];
            auto sub_neighbors_size = region_query_no_dyn(points, n, eps, num_points_to_cluster, sub_neighbors);

            if (sub_neighbors_size >= minPts)
            {
                for(size_t j=0; j < sub_neighbors_size; j++)
                {
                    auto sn = sub_neighbors[j];
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

        size_t current_neighbors[MAX_NUM_NEIGHBORS];
        auto currentNeighbors_size = region_query_no_dyn(points, current, eps, num_points_to_cluster, current_neighbors);
        if (currentNeighbors_size >= minPts)
        {
            for(int k = 0; k < currentNeighbors_size; k++)
            {
                auto n = current_neighbors[k];
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

int dbscan_algo_fixed_mem(Point points[], double eps, double min_points, size_t num_points_to_cluster)
{
    int cluster_id = 0;
    cluster clusters[MAX_CLUSTERS];

    for(int i=0; i<num_points_to_cluster; i++)
    {
        if(points[i].label == PointLabel::UNVISITED)
        {
            expand_cluster_no_dyn(points, i, cluster_id, eps, min_points, clusters, num_points_to_cluster);
            if (points[i].cluster_id == cluster_id)
            {
                ++cluster_id;
            }
        }
    }
    return cluster_id;
}
