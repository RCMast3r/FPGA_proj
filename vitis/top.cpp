#include <header.h>
#include "hls_stream.h"

/**
 * @brief This functions reads in lines from the input file in DRAM over AXI into a stream and removes empty events
 * 
 * @param input_file_lines array of input file lines each stored in the "fired_pixel" struct
 * @param num_lines the number of lines in the input file to read
 * @param fired_pixel_stream_out the stream of fired pixels, non-empty events, and end markers
 */
void read_input_lines(fired_pixel input_file_lines[], unsigned int num_lines, hls::stream<fired_pixel>& fired_pixel_stream_out) {
    fired_pixel prev_line = input_file_lines[0];

    for (int i = 1; i < num_lines; i++) {
        fired_pixel curr_line = input_file_lines[i];

        bool prev_is_pixel = !(prev_line.is_new_event);
        bool prev_not_empty_event = !(prev_line.is_new_event && curr_line.is_new_event);

        // see if we can send prior_line
        if (prev_is_pixel || prev_not_empty_event) {
            fired_pixel_stream_out << prev_line;
        }

        prev_line = curr_line;
    }

    fired_pixel stream_end_marker;
    stream_end_marker.is_end = 1;
    fired_pixel_stream_out << stream_end_marker;

    return;
}

void HLS_kernel_columnar_cluster(fired_pixel input_file_lines[], unsigned int num_lines, cluster clusters[]) {
    //#pragma hls interface …

    hls::stream<fired_pixel> fired_pixel_stream_A, fired_pixel_stream_B, fired_pixel_stream_C;
    hls::stream<cluster_bounds> subcluster_stream, cluster_bounds_stream;
    hls::stream<cluster> cluster_stream;

    #pragma hls dataflow
    // Read Stage
    // Stage 1 – Column Pair Clustering
    // Stage 2 – Cluster Stitching
    // Stage 3 – Find Cluster Keys
    // Write Stage

    read_input_lines(input_file_lines, num_lines, fired_pixel_stream_A); // suppresses empty events
    //add_pixel_to_subcluster(fired_pixel_stream_A, subcluster_stream, fired_pixel_stream_B);
    //stitch_subclusters(subcluster_stream, fired_pixel_stream_B, cluster_bounds_stream, fired_pixel_stream_C);
    //analyze_clusters(cluster_bounds_stream, fired_pixel_stream_C, cluster_stream);
    //write_clusters(cluster_stream, clusters);
}
