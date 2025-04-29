#include <header.h>



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "fired_pixel.h"    // your fired_pixel definition
#include "cluster.h"        // your cluster definition

using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>\n";
        return 1;
    }

    // 1) Open input file
    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Error: cannot open " << argv[1] << "\n";
        return 1;
    }

    // 2) Read & parse lines into a vector
    std::vector<fired_pixel> inputs;
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        char flag;           // '0' or '1'
        std::string payload; // e.g. "68A98" or "09F/0EB"
        if (!(iss >> flag >> payload)) {
            std::cerr << "Skipping malformed line: " << line << "\n";
            continue;
        }

        fired_pixel fp;
        fp.is_end = 0;        // always zero for file
        if (flag == '1') {
            // new event → ID
            fp.is_new_event = 1;
            fp.ID = (ID_t)std::stoul(payload, nullptr, 16);
        } else {
            // coordinate line → col/row
            fp.is_new_event = 0;
            auto slash = payload.find('/');
            std::string col_s = payload.substr(0, slash);
            std::string row_s = payload.substr(slash + 1);
            fp.coords.col = (col_idx_t)std::stoul(col_s, nullptr, 16);
            fp.coords.row = (row_idx_t)std::stoul(row_s, nullptr, 16);
        }
        inputs.push_back(fp);
    }

    // 3) Number of lines
    unsigned int num_lines = inputs.size();
    std::cout << "Read " << num_lines << " fired_pixel entries\n";
    int Max_cluster_count = 100000000;

    // 4) Allocate output buffer (size <= num_lines)
    std::vector<cluster> clusters(Max_cluster_count);

    // 5) Call your HLS kernel
    HLS_kernel_columnar_cluster(
        inputs.data(),
        num_lines,
        clusters.data(),
        Max_cluster_count
    );



    return 0;
}
