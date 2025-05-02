#include "header.h"
using namespace std;

int main(int argc, char *argv[]) {
    // if the cosim harness gives you a filename, use it;
    // otherwise fall back to a hard-coded test vector in the cosim working directory:
    const char* fname = (argc > 1)
        ? argv[1]
        : "input_test_file.txt";   // put your vectors here

    // 1) Open input file
    ifstream infile(fname);
    if (!infile) {
        cerr << "Error: cannot open " << fname << "\n";
        return 1;
    }

    // 2) Read & parse lines into a vector
    vector<fired_pixel> inputs;
    string line;
    while (getline(infile, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        char flag;
        string payload;
        if (!(iss >> flag >> payload)) {
            cerr << "Skipping malformed line: " << line << "\n";
            continue;
        }
        fired_pixel fp;
        fp.is_end = 0;
        if (flag == '1') {
            fp.is_new_event = 1;
            fp.ID = (ID_t)stoul(payload, nullptr, 16);
        } else {
            fp.is_new_event = 0;
            auto slash = payload.find('/');
            fp.coords.col = (col_idx_t)stoul(payload.substr(0, slash), nullptr, 16);
            fp.coords.row = (row_idx_t)stoul(payload.substr(slash+1), nullptr, 16);
        }
        inputs.push_back(fp);
    }

    // 3) Report how many we got
    unsigned num_lines = inputs.size();
    cout << "Read " << num_lines << " fired_pixel entries\n";

    // 4) Prepare output buffer
    constexpr int max_clusters = 100;
    vector<cluster> clusters(max_clusters);

    // 5) Call your kernel
    HLS_kernel_columnar_cluster(
      inputs.data(),
      num_lines,
      clusters.data(),
      max_clusters
    );

    return 0;
}

