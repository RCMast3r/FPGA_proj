#include <cstdint>
#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <algorithm>
#include <unordered_map>

#define COLUMNS 1024
// #define SAMPLE_MAX 512
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
                   std::unordered_map<int, std::vector<sample>>& current_clusters, 
                   int current_cluster_index)
{
    int starting_col = samples.front().col;

    int ci = current_cluster_index;
    for(const auto & sample : samples)
    {

    }
}

std::vector<sample> cluster_column(const std::vector<sample>& samples_l, const std::vector<sample>& samples_r)
{
    std::vector<sample> cluster_res;
    for()
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