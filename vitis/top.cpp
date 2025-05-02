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
 * @brief stitches one boundary into another
 * 
 * @param source this struct has its bounds altered
 * @param addition bounds get stitched into the source
 */
void stitch_bounds(cluster_bounds& source, cluster_bounds& addition)
{
    col_idx_t L = source.bounds.L;
    col_idx_t R = source.bounds.R;
    row_idx_t T = source.bounds.T;
    row_idx_t B = source.bounds.B;

    col_idx_t aL = addition.bounds.L;
    col_idx_t aR = addition.bounds.R;
    row_idx_t aT = addition.bounds.T;
    row_idx_t aB = addition.bounds.B;

    source.bounds.L = std::min(L, aL);
    source.bounds.R = std::max(R, aR);
    source.bounds.T = std::min(T, aT); // potentially can be skipped given our data ordering
    source.bounds.B = std::max(B, aB);
}

/**
 * @brief writes unstitched subclustes left in array to the cluster_bounds stream
 * 
 * @param curr_acc_subclusters subcluster array to be flushed
 * @param cluster_bounds_stream stream to write final cluster bounds to
 */
void flush_subclusters(cluster_bounds subclusters_arr[], hls::stream<cluster_bounds>& cluster_bounds_stream)
{
    for (unsigned int i = 0; ((i < 256) && !(subclusters_arr[i].is_end)); i++)
    {
        if (!(subclusters_arr[i].is_new_event)) // is_new_event is used in these arrays as a marker for stitched scs 
        {
            cluster_bounds_stream << subclusters_arr[i];
#if DEBUG==3
            col_idx_t L = subclusters_arr[i].bounds.L;
            col_idx_t R = subclusters_arr[i].bounds.R;
            row_idx_t T = subclusters_arr[i].bounds.T;
            row_idx_t B = subclusters_arr[i].bounds.B;

            std::cout << "sent final cb: (L: " <<
                (unsigned int)L <<
                ", R: " <<
                (unsigned int)R <<
                ", T: " <<
                (unsigned int)T <<
                ", B: " <<
                (unsigned int)B <<
                ")" <<
                std::endl;
#endif
        }
    }
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
    ap_uint<8> place_idx_next_acc;

    fired_pixel first_fp_of_next_col_pair;
    cluster_bounds first_sc_of_next_col_pair;

    //col_idx_t fp_prev_C;
    //row_idx_t fp_prev_R;
    row_idx_t fp_curr_col_pair_idx;
    row_idx_t fp_prev_col_pair_idx;

    row_idx_t sc_curr_col_pair_idx;
    row_idx_t sc_prev_col_pair_idx;

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
    bool is_first_sc_of_event = true;
    bool is_first_col_pair_of_event = true;

    bool have_next_fp_buffered = false;
    bool have_next_sc_buffered = false;

    while(!sc_is_end)
    {
        //#pragma hls pipeline // don't do, else it tries to unroll
        // potential future fix: convert the three stages to functions and prevent them from inlining
#if DEBUG==3
        std::cout << "------------------------------" <<
            std::endl <<
            "begin next upper while loop, reinit" <<
            std::endl;
#endif

        // void reinit_local_variables()
        // for each in array -> set to zero
        //      acc arrays should use the is_end bit to mark empty scs 
        //      don't reset arrays that should carry over between outer loop iterations

        bool fp_in_same_col_pair = true;
        bool sc_in_same_col_pair = true;

        bool sc_is_new_event = false;

        // swap edges
        acc_right_edge_C = adj_left_edge_C + 1;
        for (unsigned int i = 0; i < sizeof(acc_right_edge); i++)
        {
            // shift edges
            bit_t zero_bit = 0; // gets around a ternary type constraint
            acc_right_edge[i] = (is_first_sc_of_event ? zero_bit : adj_right_edge[i]);
            adj_left_edge[i] = zero_bit; // assume pixels are not fired
            adj_right_edge[i] = zero_bit;
        }

        cluster_bounds empty_marked_subcluster;
        empty_marked_subcluster.is_end = 1;

        // swap in new subclusters
        for (unsigned int i = 0; i < 256; i++)
        {
            // prev iter's next clusters are now the curr clusters
            curr_acc_subclusters[i] = (is_first_sc_of_event ? empty_marked_subcluster : next_acc_subclusters[i]);
            next_acc_subclusters[i] = empty_marked_subcluster; // assume pixels are not fired
        }
        place_idx_next_acc = 0;

        cluster_bounds prior_sc_from_next_acc = empty_marked_subcluster;
        bool sc_accum_into_prior_stitch = false;

#if DEBUG==3
        std::cout << "------------------------------" <<
            std::endl <<
            "begin next fp while loop" <<
            std::endl;
#endif
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
                is_first_fp_of_event = true;
                fired_pixel_stream_out << fp;
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
                fp_prev_col_pair_idx = fp_curr_col_pair_idx;

                // pixel is in the next column pair, so save it for later
                if (!fp_in_same_col_pair && !is_first_fp_of_event) 
                {
                    first_fp_of_next_col_pair = fp;
                    fp_in_same_col_pair = false;
                    have_next_fp_buffered = true;
#if DEBUG==3
                    std::cout << "buffered the fp, as its from next col-pair" <<
                        std::endl;
#endif
                }
                else // pixel is in the same column pair, so add it to the edge array
                {
#if DEBUG==3
                    std::cout << "fp is in curr col-pair" <<
                        std::endl;
#endif
                    
                    bool is_left_edge = ((C % 2) == 0);

                    if (is_left_edge)
                    {
                        adj_left_edge[R] = 1;
                        adj_left_edge_C = C;
#if DEBUG==3
                        std::cout << "fp is in left col of pair" <<
                            std::endl;
#endif
                    }
                    else
                    {
                        adj_right_edge[R] = 1;
                        adj_left_edge_C = (C - 1);
#if DEBUG==3
                            std::cout << "fp is in right col of pair" <<
                            std::endl;
#endif
                    }

                    //fp_prev_C = C;
                    //fp_prev_R = R; // is this needed?

                    is_first_fp_of_event = false;

                    fired_pixel_stream_out << fp;
                }
            }
            // would normally put the fired_pixel_stream_out << fp; here,
            // but we shouldn't double-write the buffered fp to the out stream
        }

#if DEBUG==3
        unsigned int num_acc_pixels = 0;
        unsigned int num_left_pixels = 0;
        unsigned int num_right_pixels = 0;
        for (int i = 0; i < sizeof(adj_left_edge); i++)
        {
            num_acc_pixels   += acc_right_edge[i];
            num_left_pixels  += adj_left_edge[i];
            num_right_pixels += adj_right_edge[i];
        }

        std::cout << "acc region right edge col idx: " <<
            (unsigned int)(acc_right_edge_C) <<
            std::endl <<
            "adj region left edge col idx:  " <<
            (unsigned int)(adj_left_edge_C) <<
            std::endl;

        std::cout << "acc region right edge has " <<
            num_acc_pixels <<
            " pixels" <<
            std::endl <<
            "adj region left edge has " <<
            num_left_pixels <<
            " pixels" <<
            std::endl <<
            "adj region right edge has " <<
            num_right_pixels <<
            " pixels" <<
            std::endl;
#endif
        
#if DEBUG==3
        std::cout << "------------------------------" <<
            std::endl <<
            "begin next sc while loop" <<
            std::endl;
#endif
        // void stitch_next_subclusters()
        while (sc_in_same_col_pair)
        {
            // if we had a carry over between column pairs, 
            if (have_next_sc_buffered) 
            {
                sc = first_sc_of_next_col_pair;
                have_next_sc_buffered = false;
#if DEBUG==3
                std::cout << "loaded the buffered sc" <<
                    std::endl;
#endif
            }
            else
            {
                subcluster_stream >> sc;
#if DEBUG==3
                std::cout << "read next sc entry from stream" <<
                    std::endl;
#endif
            }

            // should move one once the column pair is done
            // - an sc from a new col pair in the same event arrives
            // - a new event is seen
            // - the stream end arrived
            if (sc.is_end)
            {
                sc_is_end = true;
                sc_in_same_col_pair = false;

                // send out reamining sc
                flush_subclusters(curr_acc_subclusters, cluster_bounds_stream);
                flush_subclusters(next_acc_subclusters, cluster_bounds_stream); // these wont carry over to next event

                cluster_bounds_stream << sc; // send along end marker

#if DEBUG==3
                std::cout << "sc read/sent end marker" <<
                    std::endl;
#endif

                // do NOT break, as new event would also do that
                // loops in dataflows can't have two exit conditions
            }
            else if (sc.is_new_event)
            {
                sc_is_new_event = true;
                sc_in_same_col_pair = false;
                have_next_sc_buffered = false;
                is_first_sc_of_event = true;
                is_first_col_pair_of_event = true;

                // send out remaining sc
                flush_subclusters(curr_acc_subclusters, cluster_bounds_stream);
                flush_subclusters(next_acc_subclusters, cluster_bounds_stream); // these wont carry over to next event

                cluster_bounds_stream << sc; // send along event marker
#if DEBUG==3
                std::cout << "sc read/sent new event id: " <<
                    std::hex <<
                    (unsigned int)(sc.ID) <<
                    std::endl;
#endif
                // do NOT break, as new event would also do that
                // loops in dataflows can't have two exit conditions

            }
            else // is a subcluster
            {
                // get the bounds of the subcluster
                col_idx_t L = sc.bounds.L;
                col_idx_t R = sc.bounds.R;
                row_idx_t T = sc.bounds.T;
                row_idx_t B = sc.bounds.B;

#if DEBUG==3
                std::cout << "sc: (L: " <<
                    std::hex <<
                    (unsigned int)(L) <<
                    ", R: " <<
                    (unsigned int)(R) <<
                    ", T: " <<
                    (unsigned int)(T) <<
                    ", B: " <<
                    (unsigned int)(B) <<
                    ")" <<
                    std::endl;
#endif        
                
                sc_curr_col_pair_idx = (R / 2);
                sc_in_same_col_pair = (sc_curr_col_pair_idx == sc_prev_col_pair_idx);
                sc_prev_col_pair_idx = sc_curr_col_pair_idx;

                if (!sc_in_same_col_pair && !is_first_sc_of_event) // pixel is in the next column pair, so save it for later
                {
                    first_sc_of_next_col_pair = sc;
                    sc_in_same_col_pair = false;
                    is_first_col_pair_of_event = false;
                    have_next_sc_buffered = true;

#if DEBUG==3
                    std::cout << "buffered the sc, as it's from next col-pair" <<
                        std::endl;
#endif

                    // send out remaining sc in acc
                    flush_subclusters(curr_acc_subclusters, cluster_bounds_stream);
                    // don't flush next_acc b/c it continues
                }
                else // sc is in the same column pair
                {
#if DEBUG==3
                    std::cout << "sc is in curr col-pair" <<
                        std::endl;
#endif


                    bool adj_touches_acc_region = ((adj_left_edge_C - acc_right_edge_C) == 1);

                    if (!adj_touches_acc_region || is_first_col_pair_of_event) 
                    {

                        // subclusters wont stitch if adj is separated from acc region
                        // or if they are in the first col pair of an event
                        // b/c there is no prior acc region to stitch to

                        // add it to the next set of acc subclusters
                        next_acc_subclusters[place_idx_next_acc] = sc;
                        place_idx_next_acc += 1;
#if DEBUG==3
                        std::cout << "sc skips stitch b/c " <<
                            (is_first_col_pair_of_event ?
                                "in first col pair of event" :
                                "separate adj & acc reg"
                            ) <<
                            std::endl;
#endif
                    }
                    else // it is possible the sc could be stitched
                    {
#if DEBUG==3
                        std::cout << "sc could be stitched" <<
                            std::endl;
#endif
                        // 1. Check if sc overlaps with prior saved next_acc sc (buffer this?)

                        if (place_idx_next_acc > 0) // don't check if there is no prior stitch
                        {
                            // get the bounds of the subcluster
                            col_idx_t pL = prior_sc_from_next_acc.bounds.L;
                            col_idx_t pR = prior_sc_from_next_acc.bounds.R;
                            row_idx_t pT = prior_sc_from_next_acc.bounds.T;
                            row_idx_t pB = prior_sc_from_next_acc.bounds.B;

                            // compare bounds to check for overlap on each axis
                            bool is_LR_overlap = (pR >= L);
                            bool is_TB_overlap = (pB >= T);
                            
                            if (is_TB_overlap && is_LR_overlap) // if they overlap
                            {
                                // stitch bounds
                                stitch_bounds(prior_sc_from_next_acc, sc);
                                next_acc_subclusters[place_idx_next_acc - 1] = prior_sc_from_next_acc;

                                sc_accum_into_prior_stitch = true;
#if DEBUG==3
                                std::cout << "sc overlapped with prior stitch & accumulated into prior" <<
                                    std::endl;
#endif
                            }
                            else // sc didn't overlap with prior switch
                            {
                                sc_accum_into_prior_stitch = false;
#if DEBUG==3
                                std::cout << "sc did not overlap with prior stitch" <<
                                    std::endl;
#endif
                            }
                        }

                        // 2. Then, if sc is not in left edge, skip checking for sc with adjacent bounds in curr_acc

                        // is it within the adj left edge?
                        bool in_adj_left_edge = ((L % 2) == 0);

                        if (!in_adj_left_edge) // sc can't stitch with anything currently
                        {
#if DEBUG==3
                            std::cout << "sc is not in the left edge" <<
                                std::endl;
#endif

                            // subclusters in the right edge can't get stitched with sc in acc region
                            // so don't bother trying to stitch it

                            if (!sc_accum_into_prior_stitch)
                            {
                                // add it to the next set of acc subclusters
                                next_acc_subclusters[place_idx_next_acc] = sc;
                                place_idx_next_acc += 1;
                                prior_sc_from_next_acc = sc;
                                
#if DEBUG==3
                                std::cout << "sc not acc. with prior nor in left, so added to next" <<
                                    std::endl;
#endif

                            }
                            // if sc was accumulated into prior stitch, then don't add it to the next acc sc
                        }
                        else // sc could stitch with acc sc
                        {
#if DEBUG==3
                            std::cout << "sc in left edge, check curr acc sc" <<
                                std::endl;
#endif
                            // 3. Check sc in acc for adjacent bounds (early skip chance when next acc sc T > sc.B)

                            bool early_exit = false;

                            cluster_bounds sc_stitch_accum = sc;

                            // go through all curr sc from the acc region
                            for (unsigned int i = 0; ((i < 256) && !(curr_acc_subclusters[i].is_end) && !early_exit); i++) // go through curr acc sc with early exit
                            {
                                col_idx_t aL = curr_acc_subclusters[i].bounds.L;
                                col_idx_t aR = curr_acc_subclusters[i].bounds.R;
                                row_idx_t aT = curr_acc_subclusters[i].bounds.T;
                                row_idx_t aB = curr_acc_subclusters[i].bounds.B;

#if DEBUG==3
                                std::cout << "checking against acc sc: (aL: " <<
                                    (unsigned int)aL <<
                                    ", aR: " <<
                                    (unsigned int)aR <<
                                    ", aT: " <<
                                    (unsigned int)aT <<
                                    ", aB: " <<
                                    (unsigned int)aB <<
                                    ")" <<
                                    std::endl;
#endif                                
                                // check if acc sc is too far below (assumes acc sc are ordered, which should be true)
                                if ((aT - B) > 1)
                                {
                                    early_exit = true;
#if DEBUG==3
                                    std::cout << "acc sc have gone past sc range, so stop searching" <<
                                        std::endl;
#endif
                                }
                                else if (aR != acc_right_edge_C) // check if acc sc is in the right edge
                                {
#if DEBUG==3
                                    std::cout << "acc sc is not in the right edge, can't be stitched" <<
                                        std::endl;
#endif
                                }
                                else
                                {
#if DEBUG==3
                                    std::cout << "acc sc is in the right edge" <<
                                        std::endl;
#endif

                                    // check for overlap of their vertical ranges
                                    bool bounds_overlap = (std::max(T, aT) <= std::min(B, aB));

                                    // check for 
                                    bool sc_bound_top_corner_adj    = ((T - aB) == 1);
                                    bool sc_bound_bottom_corner_adj = ((aT - B) == 1);

                                    bool has_adj_pixel = false;

                                    if (bounds_overlap)
                                    {
                                        // 4. If adjacent bounds, check the overlapping range for adj pixels in (min R, max R)
                                        //    Add bound if an adj pixel is found, and mark acc sc as stitched (is_new_event)
                                        
                                        // get the overlapping range
                                        row_idx_t oT = std::max(T, aT);
                                        row_idx_t oB = std::min(B, aB);

#if DEBUG==3
                                        std::cout << "acc sc has range overlap with sc: (oT: " <<
                                            (unsigned int)oT <<
                                            ", oB: " <<
                                            (unsigned int)oB <<
                                            ")" <<
                                            std::endl;
#endif

                                        // check edges
                                        //(dont check up/down if outside of overlap range or if outside of array range)

                                        // for each pixel in the overlap range of the adj sc
                                        for (unsigned int j = oT; ((j <= oB) && !has_adj_pixel); j++)
                                        {
                                            if (adj_left_edge[j]) // if adj has a fired pixel there
                                            {
                                                // check each pixel at the obove, equal, and below position in acc sc (if possible)
                                                bool is_acc_equal_fired = acc_right_edge[j]; // center is always in bounds
                                                bool is_acc_above_fired = false;
                                                if (((j - 1) >= 0) && ((j - 1) >= oT)) // safety and overlap range check
                                                {
                                                    is_acc_above_fired = acc_right_edge[j - 1];
                                                }
                                                bool is_acc_below_fired = false;
                                                if (((j + 1) < 256) && ((j + 1) <= oB)) // safety and overlap range check
                                                {
                                                    is_acc_below_fired = acc_right_edge[j + 1];
                                                }
                                                has_adj_pixel = is_acc_equal_fired || is_acc_above_fired || is_acc_below_fired;
#if DEBUG==3
                                                if (has_adj_pixel)
                                                {
                                                    std::cout << "adj_pixels: adj_sc C: " <<
                                                        (unsigned int)j <<
                                                        " w/ (above,equal,below:" <<
                                                        (unsigned int)is_acc_above_fired <<
                                                        "," <<
                                                        (unsigned int)is_acc_equal_fired <<
                                                        "," <<
                                                        (unsigned int)is_acc_below_fired <<
                                                        ")" <<
                                                        std::endl;
                                                }
#endif
                                            }
                                        }
                                    }
                                    else if (sc_bound_top_corner_adj)
                                    {
#if DEBUG==3
                                        std::cout << "sc bound top corner is adj to acc sc's" <<
                                            std::endl;
#endif
                                        bool sc_top_corner_pixel_adj = ((adj_left_edge[T] == 1) && (acc_right_edge[aB] == 1));

                                        if (sc_top_corner_pixel_adj)
                                        {
                                            has_adj_pixel = true;
#if DEBUG==3
                                            std::cout << "sc top corner pixel is adj to acc sc's" <<
                                                std::endl;
#endif
                                        }
                                    }
                                    else if (sc_bound_bottom_corner_adj)
                                    {
#if DEBUG==3
                                        std::cout << "sc bound bottom corner is adj to acc sc's" <<
                                            std::endl;
#endif
                                        bool sc_bottom_corner_pixel_adj = ((adj_left_edge[B] == 1) && (acc_right_edge[aT] == 1));

                                        if (sc_bottom_corner_pixel_adj)
                                        {
                                            has_adj_pixel = true;
#if DEBUG==3
                                            std::cout << "sc bottom corner pixel is adj to acc sc's" <<
                                                std::endl;
#endif
                                        }
                                    }

                                    if (has_adj_pixel)
                                    {
                                        // stitch bounds and add to previous stitch
                                        stitch_bounds(sc_stitch_accum, curr_acc_subclusters[i]);

                                        // mark acc sc as stitched so it doesn't get written to cluster stream
                                        curr_acc_subclusters[i].is_new_event = 1;

#if DEBUG==3
                                        std::cout << "interm. stitch accum to: (sL: " <<
                                            (unsigned int)(sc_stitch_accum.bounds.L) <<
                                            ", sR: " <<
                                            (unsigned int)(sc_stitch_accum.bounds.R) <<
                                            ", sT: " <<
                                            (unsigned int)(sc_stitch_accum.bounds.T) <<
                                            ", sB: " <<
                                            (unsigned int)(sc_stitch_accum.bounds.B) <<
                                            ")" <<
                                            std::endl <<
                                            "marked acc sc as stitched" <<
                                            std::endl;
#endif
                                    }
                                }
                            }

                            bool stitch_accum_into_prior_stitch = false;

                            // check if stitch accmum overlaps the previous stitch
                            if (place_idx_next_acc > 0)
                            {
                                // get the bounds of the prior subcluster
                                col_idx_t pL = prior_sc_from_next_acc.bounds.L;
                                col_idx_t pR = prior_sc_from_next_acc.bounds.R;
                                row_idx_t pT = prior_sc_from_next_acc.bounds.T;
                                row_idx_t pB = prior_sc_from_next_acc.bounds.B;

                                // compare bounds to check for overlap on each axis
                                bool is_LR_overlap = (pR >= sc_stitch_accum.bounds.L);
                                bool is_TB_overlap = (pB >= sc_stitch_accum.bounds.T);
                                
                                if (is_TB_overlap && is_LR_overlap) // if they overlap
                                {
                                    stitch_accum_into_prior_stitch = true;
#if DEBUG==3
                                    std::cout << "stitch accum overlapped with prior stitch & accumulated into prior" <<
                                        std::endl;
#endif
                                }
                            }

                            // check whether to add to prior stitch or start a new one
                            if (sc_accum_into_prior_stitch || stitch_accum_into_prior_stitch)
                            {
                                // add to previous stitch
                                stitch_bounds(prior_sc_from_next_acc, sc_stitch_accum);
                                next_acc_subclusters[place_idx_next_acc - 1] = prior_sc_from_next_acc;

#if DEBUG==3
                                std::cout << "saved stitch accum to prior sc in next_acc_sc: (sL: " <<
                                    (unsigned int)(prior_sc_from_next_acc.bounds.L) <<
                                    ", sR: " <<
                                    (unsigned int)(prior_sc_from_next_acc.bounds.R) <<
                                    ", sT: " <<
                                    (unsigned int)(prior_sc_from_next_acc.bounds.T) <<
                                    ", sB: " <<
                                    (unsigned int)(prior_sc_from_next_acc.bounds.B) <<
                                    ")" <<
                                    std::endl;
#endif
                            }
                            else
                            {
                                // stitch it in
                                next_acc_subclusters[place_idx_next_acc] = sc_stitch_accum;
                                place_idx_next_acc += 1;
                                prior_sc_from_next_acc = sc_stitch_accum;

#if DEBUG==3
                                std::cout << "saved stitch accum to next_acc_sc: (sL: " <<
                                    (unsigned int)(sc_stitch_accum.bounds.L) <<
                                    ", sR: " <<
                                    (unsigned int)(sc_stitch_accum.bounds.R) <<
                                    ", sT: " <<
                                    (unsigned int)(sc_stitch_accum.bounds.T) <<
                                    ", sB: " <<
                                    (unsigned int)(sc_stitch_accum.bounds.B) <<
                                    ")" <<
                                    std::endl;
#endif
                            }
                        }
                    }

                    is_first_sc_of_event = false;
                }
            }
        }
    }
#if DEBUG==3
    std::cout << "End of Stage 3" <<
        std::endl;
#endif
}

#if STAGE_1_3_TEST==1
void consume_remaining_streams(hls::stream<fired_pixel>& fired_pixel_stream_in, hls::stream<cluster_bounds>& cluster_bounds_stream)
{
    
}
#endif

void HLS_kernel_columnar_cluster(fired_pixel input_file_lines[], unsigned int num_lines, cluster clusters[])
{
    #pragma HLS interface m_axi port = input_file_lines offset = slave bundle = mem1
    //#pragma HLS interface m_axi port = num_lines offset = slave bundle = mem1 // not needed, as not a mem location
    #pragma HLS interface m_axi port = clusters offset = slave bundle = mem2
    #pragma HLS interface s_axilite port = return

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
#if STAGE_1_3_TEST==1
    consume_remaining_streams(fired_pixel_stream_C, cluster_bounds_stream);
#endif
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
