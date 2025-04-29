#include <header.h>

const bool IS_INPUT_LOCAL = false;

using namespace std;

void read_input_file(const char *filename, fired_pixel file_line_array[])
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        perror("Failed to open file");
        exit(1);
    }
}

int main () {
    // conditional addition of "../data/" to start of filename
    std::string redirect = "../data/";

    return 0;
}