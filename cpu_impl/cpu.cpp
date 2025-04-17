#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
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
    std::vector<cluster_member>
};
struct sample
{
    int row = -1;
    int col = -1;
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

std::vector<sample> cluster_column(const std::vector<sample>& samples_l, const std::vector<sample>& samples_r)
{
    std::vector<sample> cluster_res;
    
    return cluster_res;
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