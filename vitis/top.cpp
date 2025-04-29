#include <header.h>
#include "hls_stream.h"

#define DEBUG_TOP 1

/**
 * @brief This functions reads in lines from the input file in DRAM over AXI into a stream and removes empty events
 * 
 * @param input_file_lines array of input file lines each stored in the "fired_pixel" struct
 * @param num_lines the number of lines in the input file to read
 * @param fired_pixel_stream_out the stream of fired pixels, non-empty events, and end markers
 */
void read_input_lines(
    fired_pixel input_file_lines[],
    unsigned int num_lines,
    hls::stream<fired_pixel>& fired_pixel_stream_out)
{
    fired_pixel prev_line = input_file_lines[0]; // first entry

    for (int i = 1; i < num_lines; i++) // "i < num_lines" is the ONLY exit condition b/c dataflow loops can only have one exit
    {
        fired_pixel curr_line = input_file_lines[i];

        bool prev_is_pixel = !(prev_line.is_new_event); // no need to check for end marker, since lines are not in a stream
        bool prev_not_empty_event = !(prev_line.is_new_event && curr_line.is_new_event); // 2 subsequent events means the first was empty

        // send along fired pixels and NON-empty events 
        if (prev_is_pixel || prev_not_empty_event)
        {
            fired_pixel_stream_out << prev_line;
        }

        prev_line = curr_line;
    }

    fired_pixel stream_end_marker;
    stream_end_marker.is_end = 1;
    fired_pixel_stream_out << stream_end_marker; // let next stage know that the stream is over

    return;
}

/**
 * @brief This functions reads fired pixels and finds clusters within column pairs
 * 
 * @param fired_pixel_stream_in stream of fired pixels (and event/end markers)
 * @param subcluster_stream stream to output subclusters to
 * @param fired_pixel_stream_out passes along the fired pixels to the next stage
 */
void add_pixel_to_subcluster(
    hls::stream<fired_pixel>& fired_pixel_stream_in,
    hls::stream<cluster_bounds>& subcluster_stream,
    hls::stream<fired_pixel>& fired_pixel_stream_out)
{
    fired_pixel fp;
    fired_pixel_stream_in >> fp;

    cluster_bounds sc; // subcluster

    cluster_bounds sc_new_event;
    sc_new_event.is_new_event = 1;
    sc_new_event.ID = fp.ID; // we know first stream entry is the new event info
    subcluster_stream      << sc_new_event;

    fired_pixel_stream_out << fp; // first fired pixel entry will just be a new event marker, so pass it along

    col_idx_t prev_C; row_idx_t prev_R;

    bool new_sc_event = true;

    while (true) // loop MUST only have ONE exit condition ("break" in this case) for dataflow 
    {

        fired_pixel_stream_in >> fp;

        if (fp.is_end)
        {
            subcluster_stream << sc; // b/c no more pixels will be added to subcluster

            cluster_bounds subcluster_stream_end_marker;
            subcluster_stream_end_marker.is_end = 1;
            subcluster_stream << subcluster_stream_end_marker; // let outputs know about end

            fired_pixel_stream_out << fp; // always pass along the fired_pixel data (no need to alter it for next stage)
            break;
        }
        else if (fp.is_new_event)
        {
            subcluster_stream << sc; // b/c no more pixels will be added to subcluster

            new_sc_event == true; // if a new event starts, the next pixel will be the first added to the subcluster

            sc_new_event.ID = fp.ID;
            subcluster_stream << sc_new_event; // let outputs know about new events

            fired_pixel_stream_out << fp; // always pass along the fired_pixel data (no need to alter it for next stage)
        }
        else // is a fired pixel
        {
            col_idx_t C = fp.coords.col;
            row_idx_t R = fp.coords.row;

            bool new_col_pair = (C / 2) != (prev_C / 2); // is this pixel in a different column-pair than the prev pixel?
            bool not_adjacent = std::max(C - prev_C, (ap_int<11>)(R - prev_R)) > 1; // must cast right argument to match left arg

            if (new_sc_event || not_adjacent || new_col_pair) // this pixel is NOT part of prev subcluster
            {
                subcluster_stream << sc; // the prev subcluster is complete

                sc.bounds.L = C; // init new subcluster based on first fired pixel
                sc.bounds.R = C;
                sc.bounds.T = R;
                sc.bounds.B = R;
            }
            else // fired pixel is part of current subcluster
            {
                sc.bounds.L = std::min(sc.bounds.L, C);
                sc.bounds.R = std::max(sc.bounds.R, C);
                // T doesn’t change after first fired pixel
                sc.bounds.B = std::max(sc.bounds.B, R);
            }

            prev_C = C; // we need to keep track of subsequent pixel coords for adjacency tests
            prev_R = R;

            fired_pixel_stream_out << fp; // always pass along the fired_pixel data (no need to alter it for next stage)
        }
    }

    // no need to send along a subcluster to the output stream here
    // b/c the end marker condition already did so
    return;
}

/**
 * @brief Find final cluster boundaries given a streams of subclusters and pixels
 * 
 * @param subcluster_stream input stream of subclusters from column-pairs
 * @param fired_pixel_stream_in stream of fired pixels (and event/end markers)
 * @param cluster_bounds_stream output stream of final cluster boundaries
 * @param fired_pixel_stream_out passes along the fired pixels to the next stage
 */
void stitch_subclusters(
    hls::stream<cluster_bounds>& subcluster_stream,
    hls::stream<fired_pixel>& fired_pixel_stream_in,
    hls::stream<cluster_bounds>& cluster_bounds_stream,
    hls::stream<fired_pixel>& fired_pixel_stream_out)
{
    ap_uint<2> R_edge = 0, L_edge = 1, next_R_edge = 2; // avoids copying between arrays
    static bit_t edges[3][512]; // partition on the [3] dimension // is bit_t actually efficient or should we pack a uint? : static for zero init.
    static col_idx_t edge_valid_range[3]; // 0=none valid : avoids zeroing BRAM : ?also part. on the [3] dimension? : 10 bits needed, hence col type
    col_idx_t R_edge_col_idx, L_edge_col_idx, next_R_edge_col_idx;
    fired_pixel fp;
    cluster_bounds sc;
    fired_pixel_stream_in >> fp; // first stream entries will be new event markers
    subcluster_stream >> sc;
    ID_t curr_event = fp.ID;
    static cluster_bounds accum_subclusters[256];
    bool sc_end = false, fp_end = false, sc_wait = false, fp_wait = false, reinit = false, first_col_pair = true;

    while (true) // only ONE exit condition - when both streams have ended
    {
        if (!fp_end && !fp_wait)
        {
            fired_pixel_stream_in >> fp;
            if (fp.is_end)
            {
                fp_end = true;
            }
            else if (fp.is_new_event)
            {
                if (sc_wait) // if sc was waiting on fp to catch up, it no longer has to
                {
                    sc_wait = false;
                    curr_event = fp.ID;
                    reinit = true;
                }
                else // if not, then fp will have to wait on sc to catch up
                {
                    fp_wait = true;
                }
            }
            else // is a fired pixel
            {
                col_idx_t C = fp.coords.col; row_idx_t R = fp.coords.row;
                bool is_left_edge = (C % 2) == 0;
                ap_uint<2> curr_edge = (is_left_edge ? L_edge : next_R_edge);
                edges[curr_edge][R] = 1;
                edge_valid_range[curr_edge] = 1 + R; // # is the *exclusive* right extent of validity, hence + 1
            }
            fired_pixel_stream_out << fp; // pass along fired pixels
        }
        if (!sc_end && !sc_wait) {
            subcluster_stream >> sc;
            if (sc.is_end)
            {
                sc_end = true;
            }
            else if (sc.is_new_event)
            {
                if (fp_wait) // if fp was waiting on sc to catch up, it no longer has to
                {
                    fp_wait = false;
                    curr_event = sc.ID;
                    reinit = true;
                }
                else // if not, then sc will have to wait on fp to catch up
                {
                    sc_wait = true;
                }
            }
            else // is a subcluster
            {
                // TBD
            }
        }

        if (reinit)
        {
            // init local variables
            // send off remaining clusters b/c we are in a new event
            // send cluster output stream the new event marker

            // TBD
        }

        if (fp_end && sc_end)
        {
            // TBD

            // send off remaining subclusters b/c there are no more subcluster to stitch

            //cluster_bounds_stream << ; // let next stage know there are no more clusters
            break;
        }
    }

    return;
};

fired_pixel fp_array[num_lines];

void copy_fp_to_bram(hls::stream<fired_pixel>&fired_pixel_stream_C)
{

int idx=0;

while(!fired_pixel_stream_C.empty())
{
    fired_pixel_stream_C >> fp_array[idx];
    idx++;
}

}

void analyze_clusters(
    hls::stream<cluster_bounds>& cluster_bounds_stream,
    hls::stream<fired_pixel>&fired_pixel_stream_C, 
    hls::stream<cluster>& cluster_stream)
    {
        ID_t present_cluster_id;
        cluster_bounds final_cluster;
        bit_t have_to_process=0;

        bit_t key[256];

        typedef ap_uint<4> relative_pos;

        ap_uint<12> sum_r=0;
        ap_uint<10> sum_c=0;

        relative_pos rr, cc;



        for(int i=0;i<256;i++)
        {
            key[i]=0;
        }

        while(true)
        {
            cluster_bounds_stream >> final_cluster;

            if(final_cluster.is_end)
            {
                cluster cluster_end;
                cluster_end.is_end =1;
                cluster_stream << cluster_end;

                break;

            }

            else if(final_cluster.is_new_event)
            {
                present_cluster_id = final_cluster.ID;
                sum_r=0;
                sum_c=0;
                for(int i=0;i<256;i++)
                    {
                        key[i]=0;
                    }
                continue;

            }

            else
            {
                cluster output_cluster;
                output_cluster.num_fired =0;


                output_cluster.num_columns = final_cluster.bounds.R - final_cluster.bounds.L ;
                output_cluster.num_rows = final_cluster.bounds.B - final_cluster.bounds.T ;
                

                //fired_pixel fp;

                for(idx=0;idx<num_lines;idx++)
                {
                
                if(fp_array[idx].is_end)
                {
                    break;
                }



                if(fp_array[idx].is_new_event)
                {
                    if(fp_array[idx].ID == present_cluster_id)
                    {

                        have_to_process =1;
                        continue;
                    }
                    else
                    {
                        have_to_process =0;
                    }
                }

                else if(have_to_process)
                {
                    if(fp_array[idx].coords.col <= final_cluster.bounds.R && fp_array[idx].coords.col>= final_cluster.bounds.L && fp_array[idx].coords.row <= final_cluster.bounds.B && fp_array[idx].coords.row >=  final_cluster.bounds.T )
                    {
                        output_cluster.num_fired++;

                        rr = fp_array[idx].coords.row - final_cluster.bounds.T;
                        cc = fp_array[idx].coords.col - final_cluster.bounds.L;

                        key[rr*output_cluster.num_columns +cc] = 1;

                        sum_r += fp_array[idx].coords.row;
                        sum_c += fp_array[idx].coords.col;

                    }


                }

                }

                output_cluster.key = key;
                output_cluster.centre_of_mass_x_cord = sum_r/output_cluster.num_fired;
                output_cluster.centre_of_mass_y_cord = sum_c/output_cluster.num_fired;

                cluster_stream << output_cluster;

                



            }
        }
    }


void write_clusters(
    hls::stream<cluster> &cluster_stream,
    cluster              *out_buf,       // pointer to a pre‐allocated DRAM region
    int                    max_clusters, // size of out_buf[]
) {


    int idx = 0;
    cluster c;
    
    // keep reading until we see the end marker
    while (true) {
    #pragma HLS PIPELINE II=1
        if (!cluster_stream.empty()) {
            cluster_stream >> c;
            if (c.is_end) {
                break;
            }
            // guard against overflowing the DRAM buffer
            if (idx < max_clusters) {
                out_buf[idx] = c;
                idx++;
            }
        }
    }
    

}

void HLS_kernel_columnar_cluster(fired_pixel input_file_lines[], unsigned int num_lines, cluster clusters[], Max_cluster_count)
{
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
    add_pixel_to_subcluster(fired_pixel_stream_A, subcluster_stream, fired_pixel_stream_B);
    //stitch_subclusters(subcluster_stream, fired_pixel_stream_B, cluster_bounds_stream, fired_pixel_stream_C);
    //analyze_clusters(cluster_bounds_stream, fired_pixel_stream_C, cluster_stream);
    //write_clusters(cluster_stream, clusters, Max_cluster_count);
}
