#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <optional>
#include <regex>
#include <string>
#include <algorithm>
#include <unordered_map>

#define COLUMNS 1024
// #define SAMPLE_MAX 512

struct cluster_member
{
    int zig_zag_index = -1;
    uint16_t row;
    uint16_t col;
};

struct cluster
{
    std::vector<cluster_member> members;
};



struct point {
    int x = -1;
    int y = -1;
};
struct bbox {
    point top_left;
    point bottom_right;
};
struct sample
{
    int row = -1;
    int col = -1;
};
struct subcluster
{
    bbox bound_box;
    std::vector<sample> members;
    int cluster_id = -1;
    int event_id = -1;
};


int hexStringToInt(const std::string& hexString) {
    int result;
    std::stringstream ss;
    ss << std::hex << hexString;
    ss >> result;
    return result;
}

void process_event(const std::vector<sample>& samples, 
                   std::vector<std::vector<sample>>& col_pair_sets,
                   int current_cluster_index)
{
    
    int prev_col_pair_index = (samples.front().col / 2);
    auto prev_zig_zag_index = 2*samples.front().row + ((samples.front().row + samples.front().row)%2);

    int ci = current_cluster_index;
    
    std::vector<sample> col_pair_set;
    // 
    for(const auto & sample : samples)
    {
        auto col_pair_index = sample.col / 2;
        auto zig_zag_index = (2*sample.row) + ((sample.row + sample.row)%2);

        if(prev_col_pair_index != col_pair_index)
        {
            col_pair_sets.push_back(col_pair_set);
            col_pair_set.clear();
        } else {
            col_pair_set.push_back(sample);
        }
        prev_col_pair_index = col_pair_index;
    }

    // stitching the column pair sets

    std::vector<std::vector<sample>> clusters;

    for(size_t ind = 0; ind < (col_pair_sets.size()-1); ind+=2)
    {
        auto col_r = col_pair_sets[ind+1];
        auto col_l = col_pair_sets[ind];
    }
    
}



// this is done on each events samples
std::vector<subcluster> cluster_column(const std::vector<sample>& samples, int & current_cluster_index, int event_id)
{
    std::vector<subcluster> cluster_res;

    std::optional<sample> prev_sample = std::nullopt;

    subcluster current_subcluster;
    for (const auto sample : samples)
    {
        auto column = sample.col;
        auto row = sample.row;
        if(!prev_sample) // first sample, must be new sub-cluster
        {   
            current_cluster_index++;
            current_subcluster.event_id = event_id;
            current_subcluster.cluster_id = current_cluster_index;
            current_subcluster.bound_box.bottom_right.x = sample.col;
            current_subcluster.bound_box.bottom_right.y = sample.row;
            current_subcluster.bound_box.top_left = current_subcluster.bound_box.bottom_right;
            current_subcluster.members.push_back(sample);
            prev_sample = sample;
        } else {
            bool new_col_pair = ( (column/2) != ((prev_sample->col)/2));
            bool new_subcluster = std::max(column - (prev_sample->col), row-(prev_sample->row)) > 1;
            if(new_subcluster || new_col_pair)
            {
                cluster_res.push_back(current_subcluster);
                current_subcluster = {};
                current_cluster_index++;
                current_subcluster.event_id = event_id;
                current_subcluster.cluster_id = current_cluster_index;
                current_subcluster.bound_box.bottom_right.x = sample.col;
                current_subcluster.bound_box.bottom_right.y = sample.row;
                current_subcluster.bound_box.top_left = current_subcluster.bound_box.bottom_right;
                current_subcluster.members.push_back(sample);

            } else {
                current_subcluster.members.push_back(sample);
                current_subcluster.bound_box.bottom_right.x = std::max(current_subcluster.bound_box.bottom_right.x, sample.col);
                current_subcluster.bound_box.top_left.x = std::min(current_subcluster.bound_box.top_left.x, sample.col);
                current_subcluster.bound_box.bottom_right.y = std::max(current_subcluster.bound_box.bottom_right.y, sample.row);
                current_subcluster.bound_box.top_left.y = std::min(current_subcluster.bound_box.top_left.y, sample.row);
            }
        }
    }
    cluster_res.push_back(current_subcluster);
    return cluster_res;
}

void stich_subclusters(std::vector<subcluster> subclusters)
{
    // first we go through all of the subclusters to see if any of the bboxes are connected

    // for(const auto sc : subclusters)
    // {
    //     sc.bound_box.bottom_right
    // }

}
int main(int argc, char *argv[]) {

    if (argc != 2)
    {
        // show error msg - like usage program.exe filename
        std::cout << "ERROR: missing relative filepath" <<std::endl;
        return 1;
    }
    // Path to the input text file
    std::string filename = argv[1];

    // Open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;
    }

    // Define your regex pattern (for example: match lines with "error")
    std::regex pattern("0\\ ([A-Fa-f0-9]{3})/([A-Fa-f0-9]{3})", std::regex_constants::icase); // case-insensitive
    std::string line;
    int line_num = 0;

    // Read and process each line

    while (std::getline(file, line)) {
        ++line_num;
        std::vector<sample> events_sample;
        std::smatch match;
        if (std::regex_match(line,match, pattern)) {
            std::cout << "Line " << line_num << ": " << hexStringToInt(match[1]) <<" "<< hexStringToInt(match[2]) << std::endl;
        } else {
            std::cout << "miss line" << line_num << ": " << line << std::endl;
        }
    }

    file.close();
    return 0;
}