#include "header.h"
#include "hls_stream.h"

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

#if DEBUG==1
            if (prev_is_pixel)
            {
                std::cout << "pixel (C: " <<
                    std::hex <<
                    (unsigned int)(prev_line.coords.col) <<
                    ", R: " <<
                    (unsigned int)(prev_line.coords.row) <<
                    ")" <<
                    std::endl;
            }
            else
            {
                std::cout << "(non-empty) event ID: " <<
                std::hex <<
                (unsigned int)prev_line.ID <<
                std::endl;
            }
#endif
        }

        prev_line = curr_line;
    }

    // check if the last line was a fired pixel, if so, send it off too
    if (!(prev_line.is_new_event))
    {
        fired_pixel_stream_out << prev_line;

#if DEBUG==1
        std::cout << "Final Case: sent pixel (C: " <<
            std::hex <<
            (unsigned int)(prev_line.coords.col) <<
            ", R: " <<
            (unsigned int)(prev_line.coords.row) <<
            ")" <<
            std::endl;
#endif
    }

    fired_pixel stream_end_marker;
    stream_end_marker.is_end = 1;
    fired_pixel_stream_out << stream_end_marker; // let next stage know that the stream is over

    return;
}

#if DEBUG>=2
void log_sent_sc(cluster_bounds sc)
{
    if (sc.is_end)
    {
        std::cout << "sent end marker" << std::endl;
    }
    else if (sc.is_new_event)
    {
        std::cout << "sent event ID: " <<
            std::hex <<
            (unsigned int)(sc.ID) <<
            std::endl;
    }
    else
    {
        std::cout << "sent subcluster (L: " <<
            std::hex <<
            (unsigned int)(sc.bounds.L) <<
            ", R: " <<
            (unsigned int)(sc.bounds.R) <<
            ", T: " <<
            (unsigned int)(sc.bounds.T) <<
            ", B: " <<
            (unsigned int)(sc.bounds.B) <<
            ")" <<
            std::endl;
    }
}
#endif

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
    subcluster_stream << sc_new_event;
#if DEBUG==2
    log_sent_sc(sc_new_event);
#endif

    fired_pixel_stream_out << fp; // first fired pixel entry will just be a new event marker, so pass it along

    col_idx_t prev_C; row_idx_t prev_R;

    bool new_sc_event = true;

    while (true) // loop MUST only have ONE exit condition ("break" in this case) for dataflow 
    {
        fired_pixel_stream_in >> fp;

        if (fp.is_end)
        {
            subcluster_stream << sc; // b/c no more pixels will be added to subcluster
#if DEBUG==2
            log_sent_sc(sc);
#endif

            cluster_bounds subcluster_stream_end_marker;
            subcluster_stream_end_marker.is_end = 1;
            subcluster_stream << subcluster_stream_end_marker; // let outputs know about end
#if DEBUG==2
            log_sent_sc(subcluster_stream_end_marker);
#endif

            fired_pixel_stream_out << fp; // always pass along the fired_pixel data (no need to alter it for next stage)
            break;
        }
        else if (fp.is_new_event)
        {
            subcluster_stream << sc; // b/c no more pixels will be added to subcluster
#if DEBUG==2
            log_sent_sc(sc);
#endif

            new_sc_event = true; // if a new event starts, the next pixel will be the first added to the subcluster

            sc_new_event.ID = fp.ID;
            subcluster_stream << sc_new_event; // let outputs know about new events
#if DEBUG==2
            log_sent_sc(sc_new_event);
#endif

            fired_pixel_stream_out << fp; // always pass along the fired_pixel data (no need to alter it for next stage)
        }
        else // is a fired pixel
        {
            col_idx_t C = fp.coords.col;
            row_idx_t R = fp.coords.row;

            bool new_col_pair = (C / 2) != (prev_C / 2); // is this pixel in a different column-pair than the prev pixel?
            bool not_adjacent = std::max(C - prev_C, (ap_int<11>)(R - prev_R)) > 1; // must cast right argument to match left arg

            if (new_sc_event) // this pixel is NOT part of prev subcluster
            {
                // no need to write the last sc b/c the is_new_event check already did that 

                new_sc_event = false;

                sc.bounds.L = C; // init new subcluster based on first fired pixel
                sc.bounds.R = C;
                sc.bounds.T = R;
                sc.bounds.B = R;
            }
            else if (not_adjacent || new_col_pair) // this pixel is NOT part of prev subcluster
            {
                subcluster_stream << sc; // the prev subcluster is complete
#if DEBUG==2
                log_sent_sc(sc);
#endif

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
    // we have a growing accumulated (acc) region
    // which grows as we append adjacent (adj) column pair regions to the right

    // we don't need to know the left edge of the accumulated region
    col_idx_t acc_right_edge_C;
    bit_t acc_right_edge[512];

    col_idx_t adj_left_edge_C; // no need to store the right col idx, as it is this plus one
    bit_t adj_left_edge[512];
    bit_t adj_right_edge[512]; // not needed for current stitching, but will be used as acc right edge next time

    cluster_bounds curr_acc_subclusters[256];
    cluster_bounds next_acc_subclusters[256];

    fired_pixel first_fp_of_next_col_pair;
    cluster_bounds first_sc_of_next_col_pair;

    col_idx_t fp_prev_C;
    row_idx_t fp_prev_R;
    row_idx_t fp_curr_col_pair_idx;
    row_idx_t fp_prev_col_pair_idx;


    fired_pixel fp;
    cluster_bounds sc;

    // -------------------

    // initially read out the first event markers & pass along
    fired_pixel_stream_in >> fp;
    fired_pixel_stream_out << fp;
    subcluster_stream >> sc;
    cluster_bounds_stream << sc;

    //bool fp_is_end = false; // no need to track since sc stream will be sync with fp stream
    bool sc_is_end = false;
    //bool is_end = false;

    bool is_first_fp_of_event = true;

    bool have_next_fp_buffered = false;
    bool have_next_sc_buffered = false;

    while(!sc_is_end)
    {
        #pragma hls pipeline

        // void reinit_local_variables()
        // for each in array -> set to zero
        //      acc arrays should use the is_end bit to mark empty scs 
        //      don't reset arrays that should carry over between outer loop iterations

        bool fp_in_same_col_pair = true;
        bool sc_in_current_col_pair = true;

        bool sc_is_new_event = false;



        // void read_next_col_pair()
        while (fp_in_same_col_pair)
        {
            // if we had a carry over between column pairs, use that one before 
            if (have_next_fp_buffered) 
            {
                fp = first_fp_of_next_col_pair;
                have_next_fp_buffered = false;
#if DEBUG==3
                std::cout << "loaded the buffered fp" <<
                    std::endl;
#endif
            }
            else
            {
                fired_pixel_stream_in >> fp;
#if DEBUG==3
                std::cout << "read next fp entry from stream" <<
                    std::endl;
#endif
            }

            // should move one once the column pair is done
            // - an fp from a new col pair in the same event arrives
            // - a new event is seen
            // - the stream end arrived
            if (fp.is_end)
            {
                // a stream end in the first column means we
                //fp_is_end = true;
                fp_in_same_col_pair = false;
                fired_pixel_stream_out << fp;
#if DEBUG==3
                std::cout << "fp read/sent end marker" <<
                    std::endl;
#endif
            }
            else if (fp.is_new_event)
            {
                // no need to set a flag to track this
                // because the sc stream will be in sync with this one and encounter the event too
                fp_in_same_col_pair = false;
                have_next_fp_buffered = false;
                fired_pixel_stream_out << fp;
                is_first_fp_of_event = true;
#if DEBUG==3
                std::cout << "fp read/sent new event id: " <<
                    std::hex <<
                    (unsigned int)(fp.ID) <<
                    std::endl;
#endif
            }
            else // is a fired pixel
            {
                // get the coordinates of the fired pixel
                col_idx_t C = fp.coords.col;
                row_idx_t R = fp.coords.row;

#if DEBUG==3
                std::cout << "fp: (C: " <<
                    std::hex <<
                    (unsigned int)(C) <<
                    ", R: " <<
                    (unsigned int)(R) <<
                    ")" <<
                    std::endl;
#endif

                // is this pixel in a different column-pair than the prev pixel?
                fp_curr_col_pair_idx = (C / 2); 
                fp_in_same_col_pair = (fp_curr_col_pair_idx == fp_prev_col_pair_idx);

                // pixel is in the next column pair, so save it for later
                if (!fp_in_same_col_pair && !is_first_fp_of_event) 
                {
#if DEBUG==3
                    std::cout << "buffered the fp, as its from next col-pair" <<
                        std::endl;
#endif
                    first_fp_of_next_col_pair = fp;
                    fp_in_same_col_pair = false;
                }
                else // pixel is in the same column pair, so add it to the edge array
                {
#if DEBUG==3
                    std::cout << "fp is in curr col-pair" <<
                        std::endl;
#endif
                    
                    bool is_left_edge ((C % 2) == 0);

                    if (is_left_edge)
                    {
                        adj_left_edge[R] = 1;
#if DEBUG==3
                        std::cout << "fp is in left col of pair" <<
                            std::endl;
#endif
                    }
                    else
                    {
                        adj_right_edge[R] = 1;
#if DEBUG==3
                            std::cout << "fp is in right col of pair" <<
                            std::endl;
#endif
                    }

                    fp_prev_C = C;
                    //fp_prev_R = R; // is this needed?
                    fp_prev_col_pair_idx = fp_curr_col_pair_idx;

                    is_first_fp_of_event = false;

                    fired_pixel_stream_out << fp;
                }
            }
            // would normally put the fired_pixel_stream_out << fp; here,
            // but we shouldn't double-write the buffered fp to the out stream
        }
        
        // void stitch_next_subclusters()
        while (sc_in_current_col_pair)
        {
            // if we had a carry over between column pairs, 
            if (have_next_sc_buffered) 
            {
                sc = first_sc_of_next_col_pair;
            }
            else
            {
                subcluster_stream >> sc;
            }

            // should move one once the column pair is done
            // - an sc from a new col pair in the same event arrives
            // - a new event is seen
            // - the stream end arrived
            if (sc.is_end)
            {
                // a stream end in the first column means we
                sc_is_end = true;
                sc_in_current_col_pair = false;

                // send out reamining sc

                // send along end marker

                // do NOT break, as new event would also do that
                // loops in dataflows can't have two exit conditions
            }
            else if (sc.is_new_event)
            {
                sc_is_new_event = true;
                sc_in_current_col_pair = false;
                have_next_sc_buffered = false;

                // send out remaining sc

                // send along event marker

                // do NOT break, as new event would also do that
                // loops in dataflows can't have two exit conditions
            }
            else // is a subcluster
            {                
                sc_in_current_col_pair; // TBD

                if (!sc_in_current_col_pair) // pixel is in the next column pair, so save it for later
                {
                    first_sc_of_next_col_pair = sc;
                    fp_in_same_col_pair = false;

                    //send out remaining sc in acc
                }
                else // sc is in the same column pair
                {
                    // is it within the adj left edge?
                    bool in_left_edge; //TBD
                    if (!in_left_edge) {
                        // subclusters in the right edge can't get stitched
                        // so don't bother trying to stitch it

                        // add it to the next set of acc subclusters
                    }
                    else // it is possible the sc could be stitched
                    {
                        // go through all curr sc from the acc region
                        bool acc_sc_is_end = false;
                        unsigned int acc_sc_idx = 0;
                        do
                        {
                            cluster_bounds acc_sc = curr_acc_subclusters[acc_sc_idx];
                            acc_sc_is_end = acc_sc.is_end;
                            if (!acc_sc_is_end) // we have another subcluster from the acc. region
                            {
                                
                            }
                            acc_sc_idx++;
                            if (acc_sc_idx >= 256)
                            {
                                acc_sc_is_end = true;
                            }
                        }
                        while(!acc_sc_is_end);
                    }
                }
            }
        }

        // void handle_end_conditions()
        if (sc_is_end)
        {
            // if the streams are over, no need to reinit or manage arrays
            // stitch_sc has already written out the remaining sc
        }
        else
        {
            if (sc_is_new_event) // new event means we may need to reset some local variables
            {
                // do we need to do anything
                // stitch_sc has already written out the remaining sc
            }
            else // the next column pair is in the same event, so manage variables as normal
            {
                // the adj_right_edge array will become the acc_left_edge, so do a copy
                // go through the sc arrays, send out sc's that won't be in the acc right edge
                // b/c they can't be stitched next time
                // also 
            }
        }
    }
}

/**
 * @brief OLD Version Find final cluster boundaries given a streams of subclusters and pixels
 * 
 * @param subcluster_stream input stream of subclusters from column-pairs
 * @param fired_pixel_stream_in stream of fired pixels (and event/end markers)
 * @param cluster_bounds_stream output stream of final cluster boundaries
 * @param fired_pixel_stream_out passes along the fired pixels to the next stage
 */
void stitch_subclusters_OLD(
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

    col_idx_t prev_C;

    // Stream synchronization booleans
    bool sc_end = false, fp_end = false, sc_wait_for_fp_to_get_event = false, fp_wait_for_sc_to_get_event = false;
    
    bool sc_waiting_for_next_col_pair = false;
    bool sc_waiting_for_edge_validity = false;
    bool fp_has_resolved_edge_validity = false;

    bool fp_waiting_for_stiching_to_finish = false;

    // Other state booleans
    bool fp_first_pixel_of_new_event = true;

    bool reinit = false, first_col_pair = true;

    fired_pixel buffer_fp;
    cluster_bounds buffer_sc;

    while (true) // only ONE exit condition - when both streams have ended
    {
        bool fp_is_waiting = fp_wait_for_sc_to_get_event;
        if (!fp_end && !fp_is_waiting)
        {
            fired_pixel_stream_in >> fp;
            if (fp.is_end)
            {
                fp_end = true;
            }
            else if (fp.is_new_event)
            {
                if (sc_wait_for_fp_to_get_event) // if sc was waiting on fp to catch up, it no longer has to
                {
                    sc_wait_for_fp_to_get_event = false;
                    curr_event = fp.ID;
                    reinit = true;
                }
                else // if not, then fp will have to wait on sc to catch up
                {
                    fp_wait_for_sc_to_get_event = true;
                }
            }
            else // is a fired pixel
            {
                // TBD: If the fired pixel stream gets to a new column-pair ahead of the stitching
                //      Save the fired_pixel b/c we won't have an edge array to add it to yet
                //      Signal a need to change edges and pause the fired pixel ingest

                col_idx_t C = fp.coords.col; row_idx_t R = fp.coords.row;

                bool new_col_pair = (C / 2) != (prev_C / 2); // is this pixel in a different column-pair than the prev pixel?

                if (fp_first_pixel_of_new_event)
                {
                    // skip the new_col_pair_check
                }

                bool is_left_edge = ((C % 2) == 0);
                ap_uint<2> curr_edge = (is_left_edge ? L_edge : next_R_edge);
                edges[curr_edge][R] = 1;
                edge_valid_range[curr_edge] = 1 + R; // # is the *exclusive* right extent of validity, hence + 1

                prev_C = C;
            }
            fired_pixel_stream_out << fp; // pass along fired pixels
        }
        
        bool sc_is_waiting = sc_wait_for_fp_to_get_event || sc_waiting_for_edge_validity || sc_waiting_for_next_col_pair;
        if (!sc_end && !sc_is_waiting) {
            subcluster_stream >> sc;
            if (sc.is_end)
            {
                sc_end = true;
            }
            else if (sc.is_new_event)
            {
                if (fp_wait_for_sc_to_get_event) // if fp was waiting on sc to catch up, it no longer has to
                {
                    fp_wait_for_sc_to_get_event = false;
                    curr_event = sc.ID;
                    reinit = true;
                }
                else // if not, then sc will have to wait on fp to catch up
                {
                    sc_wait_for_fp_to_get_event = true;
                }
            }
            else // is a subcluster
            {
                // TBD
            }
        }

        if (reinit) // occrs when both streams have hit a new event
        {
            // init local variables
            fp_first_pixel_of_new_event = true;

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

void HLS_kernel_columnar_cluster(fired_pixel input_file_lines[], unsigned int num_lines, cluster clusters[])
{
    //#pragma hls interface …

    hls::stream<fired_pixel> fired_pixel_stream_A, fired_pixel_stream_B, fired_pixel_stream_C;
    hls::stream<cluster_bounds> subcluster_stream, cluster_bounds_stream;
    hls::stream<cluster> cluster_stream;

    #pragma hls dataflow
    // Stage 1 - Read to Stream
    // Stage 2 – Column Pair Clustering
    // Stage 3 – Cluster Stitching
    // Stage 4 – Find Cluster Keys
    // Stage 5 - Write from Stream

    read_input_lines(input_file_lines, num_lines, fired_pixel_stream_A); // suppresses empty events
    add_pixel_to_subcluster(fired_pixel_stream_A, subcluster_stream, fired_pixel_stream_B);
    stitch_subclusters(subcluster_stream, fired_pixel_stream_B, cluster_bounds_stream, fired_pixel_stream_C);
    //analyze_clusters(cluster_bounds_stream, fired_pixel_stream_C, cluster_stream);
    //write_clusters(cluster_stream, clusters);
}

// Stage Debugging for C-sim

#ifdef DEBUG
void debug_stage(fired_pixel input_file_lines[], unsigned int num_lines)
{
    //#pragma hls interface …
#if DEBUG >= 1
    hls::stream<fired_pixel> fired_pixel_stream_A, fired_pixel_stream_B, fired_pixel_stream_C;
#endif

#if DEBUG >= 2
    hls::stream<cluster_bounds> subcluster_stream, cluster_bounds_stream;
#endif

#if DEBUG >= 4
    hls::stream<cluster> cluster_stream;
#endif

    #pragma hls dataflow
    // Stage 1 - Read to Stream
    // Stage 2 – Column Pair Clustering
    // Stage 3 – Cluster Stitching
    // Stage 4 – Find Cluster Keys
    // Stage 5 - Write from Stream
#if DEBUG >= 1
    read_input_lines(input_file_lines, num_lines, fired_pixel_stream_A); // suppresses empty events
#endif

#if DEBUG >= 2
    add_pixel_to_subcluster(fired_pixel_stream_A, subcluster_stream, fired_pixel_stream_B);
#endif

#if DEBUG >= 3
    stitch_subclusters(subcluster_stream, fired_pixel_stream_B, cluster_bounds_stream, fired_pixel_stream_C);
#endif

#if DEBUG >= 4
    analyze_clusters(cluster_bounds_stream, fired_pixel_stream_C, cluster_stream);
#endif

#if DEBUG >= 5
    write_clusters(cluster_stream, clusters);
#endif
}

#endif
