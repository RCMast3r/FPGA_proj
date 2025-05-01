#ifndef __CPU_NODYNMEM_H__
#define __CPU_NODYNMEM_H__

#include "types.hpp"

double distance(const Point &a, const Point &b);

size_t region_query_no_dyn(Point points[], int idx, double eps, size_t num_points_to_cluster, size_t neighbors[]);
void expand_cluster_no_dyn(Point points[], int idx, int clusterID, double eps, int minPts, cluster clusters[], size_t num_points_to_cluster);
int dbscan_algo_fixed_mem(Point points[], double eps, double min_points, size_t num_points_to_cluster);
#endif // __CPU_NODYNMEM_H__