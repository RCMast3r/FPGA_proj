#ifndef COL_CLUST_HEADER
#define COL_CLUST_HEADER

#include <ap_int.h>
#include <cstdint>
//#include <cstdlib>
//#include <fstream>
//#include <hls_math.h>
//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
#include <algorithm>

typedef ap_uint<10> col_idx_t; // Column Index: 1024 Columns -> log2(1024) = 10
typedef ap_uint<9>  row_idx_t; // Row Index: 512 Rows -> log2(512) = 9
typedef ap_uint<19> ID_t; // Event ID Number: ID seems to go up to 7FFFF -> ceil(log2(0x7FFFF)) = 19
typedef ap_uint<1>  bit_t; // Single Bit: I read that Vitis treats booleans as 8 bits in a struct, so this may pack better.

/**
 * @brief stores column and row indices with minimum bits necessary
 */
struct coordinates {
    col_idx_t col; // Column Index
    row_idx_t row; // Row Index
};
/**
 * @brief stores column and row indices with minimum bits necessary
 * 
 * wrapper of "coordinates" type with stream information
 */
struct fired_pixel {
    bit_t is_end; // marks the end of streamed data
    /**
     * @brief marks a new event
     * 
     * also indicates the type stored in the anonymous union
     * 
     * 0: coordinates
     * 1: ID_t 
     */
    bit_t is_new_event;
    union {
        coordinates coords; // coordinates of the fired pixel in the event grid
        ID_t ID; // ID of the new event
    };

    // we need a default constructor b/c compiler can't know which union member we want to initialize
    // data is not an event nor end marker by default, so initialize primary data type
    fired_pixel() : is_end{}, is_new_event{}, coords{} {}

    // copy assignment operator
    fired_pixel& operator=(const fired_pixel& other) {
        if (this != &other) {  // protect against self-assignment
            is_end = other.is_end;
            is_new_event = other.is_new_event;

            if (is_new_event == 0) {
                coords = other.coords;
            } else {
                ID = other.ID;
            }
        }
        return *this;
    }

    // equality test
    bool operator==(const fired_pixel& other) const {
        if (is_end != other.is_end || is_new_event != other.is_new_event) {
            return false;
        }
        if (is_new_event == 0) {
            return (
                coords.col == other.coords.col &&
                coords.row == other.coords.row
            );
        }
        else { return ID == other.ID; }
    }

    // non-equality test
    bool operator!=(const fired_pixel& other) const {
        return !(*this == other);
    }
};

/**
 * @brief stores the row and indices of the bounds of a (sub)cluster
 */
struct box_bounds {
    col_idx_t L;
    col_idx_t R;
    row_idx_t T;
    row_idx_t B;
};
/**
 * @brief stores the row and indices of the bounds of a (sub)cluster
 * 
 * wrapper of "box_bounds" type with stream information
 */
struct cluster_bounds {
    bit_t is_end; // marks the end of streamed data
    /**
     * @brief marks a new event
     * 
     * also indicates the type stored in the anonymous union
     * 
     * 0: coordinates
     * 1: ID_t 
     */
    bit_t is_new_event;
    union {
        box_bounds bounds; // indices of the (sub)cluster's bounding box
        ID_t ID; // ID of the new event
    };

    // we need a default constructor b/c compiler can't know which union member we want to initialize
    // data is not an event nor end marker by default, so initialize primary data type
    cluster_bounds() : is_end{}, is_new_event{}, bounds{} {}

    // copy assignment operator
    cluster_bounds& operator=(const cluster_bounds& other) {
        if (this != &other) {  // protect against self-assignment
            is_end = other.is_end;
            is_new_event = other.is_new_event;

            if (is_new_event == 0) {
                bounds = other.bounds;
            } else {
                ID = other.ID;
            }
        }
        return *this;
    }

    // equality test
    bool operator==(const cluster_bounds& other) const {
        if (is_end != other.is_end || is_new_event != other.is_new_event) {
            return false;
        }
        if (is_new_event == 0) {
            return (
                bounds.L == other.bounds.L &&
                bounds.R == other.bounds.R &&
                bounds.T == other.bounds.T &&
                bounds.B == other.bounds.B
            );
        } else { return ID == other.ID; }
    }

    // non-equality test
    bool operator!=(const cluster_bounds& other) const {
        return !(*this == other);
    }
};

// Topic 3 says max. fired pixels in a cluster is 16
// log2(16) = 4
#define fired_pixels_per_cluster_bits 4
const unsigned int max_num_fired_pixels_per_cluster = (1 << fired_pixels_per_cluster_bits);

// Largest bounding box for a # of fired pixels is if the pixels are in a diagonal line
// 16 x 16 box -> 256 pixels in the largest bounding box
const unsigned int max_num_pixels_per_cluster = max_num_fired_pixels_per_cluster * max_num_fired_pixels_per_cluster;

// minimum bits needed to store the maximum number of fired pixels 
// add one to get the actual number of bits (saves bitwidth)
// used to know how far to index into the key array
typedef ap_uint<fired_pixels_per_cluster_bits> num_fired_t;

//typedef ap_uint<(max_num_pixels_per_cluster)> key_t; // stores which pixels in the cluster were fired (book ordering)

// the default key as described in the PowerPoint is actually inefficient b/c only 16 of the 256 bits will ever be set
// instead, have an array of 16 indices (where the index is the book ordering of the bounding box)
// index should be 8 bits to represent the 256 possible positions
// Nice thing about using the index instead of coords is that it can adapt to varying bounding box dimensions
typedef ap_uint<(2 * fired_pixels_per_cluster_bits)> box_bounds_idx_t;

/**
 * @brief stores all information about a pixel cluster
 */
struct cluster {
    ID_t ID; // Event ID Number
    box_bounds bounds; // bounding box of the cluster
    num_fired_t num_fired; // number of fired pixels in the cluster
    box_bounds_idx_t key[max_num_fired_pixels_per_cluster]; // locations of the fired pixels in the cluster
    // tbd? add center of mass AKA centroid
};

/**
 * @brief This function finds pixel clusters given a file of events.
 * 
 * @param input_file_lines DRAM array of file lines stored in a fired_pixel struct.
 * @param num_lines The number of lines contained in the file.
 * @param clusters DRAM array which stores any found clusters.
 */
void HLS_kernel_columnar_cluster(fired_pixel input_file_lines[], unsigned int num_lines, cluster clusters[]);

#endif