#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <iomanip>

// Print up to N×N of a square matrix stored as:
// [ int64_t rows ][ int64_t cols ][ rows*cols (double) ]
// row-major, contiguous.
void print_first_N(const std::string& filename, std::size_t N) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return;
    }

    // read header
    int64_t rows, cols;
    in.read(reinterpret_cast<char*>(&rows), sizeof rows);
    in.read(reinterpret_cast<char*>(&cols), sizeof cols);
    if (!in) {
        std::cerr << "ERROR: failed reading header\n";
        return;
    }

    // reduce N to file dimensions
    std::size_t R = static_cast<std::size_t>(rows);
    std::size_t C = static_cast<std::size_t>(cols);
    std::size_t Nr = std::min(N, R);
    std::size_t Nc = std::min(N, C);

    // buffer to hold a partial row
    std::vector<double> buffer(Nc);

    // advance to the first double
    in.seekg(2 * sizeof(int64_t), std::ios::beg);

    // for each of the first Nr rows:
    for (std::size_t i = 0; i < Nr; ++i) {
        // read Nc doubles for this row
        in.read(reinterpret_cast<char*>(buffer.data()),
                Nc * sizeof(double));
        if (!in) {
            std::cerr << "ERROR: read failed at row " << i << "\n";
            return;
        }

        // skip the rest of this row (if C > Nc)
        std::size_t to_skip = (C - Nc) * sizeof(double);
        in.seekg(to_skip, std::ios::cur);

        // print them
        for (std::size_t j = 0; j < Nc; ++j) {
            std::cout << std::fixed << std::setprecision(6)
                      << buffer[j]
                      << (j + 1 < Nc ? ", " : "\n");
        }
    }
}

int main(int argc, char* argv[]) {
    const std::string prog = argv[0];
    if (argc < 2) {
        std::cerr << "Usage: " << prog << " [-h|--help] [-n N] <matrix_file>\n"
                  << "  -h, --help    Show this help message and exit\n"
                  << "  -n N          Number of rows/cols to display (default 10)\n"
                  << "  matrix_file   Binary file containing the matrix\n";
        return 1;
    }

    std::string filename;
    std::size_t N = 10;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << prog << " [-h|--help] [-n N] <matrix_file>\n"
                      << "  -h, --help    Show this help message and exit\n"
                      << "  -n N          Number of rows/cols to display (default 10)\n"
                      << "  matrix_file   Binary file containing the matrix\n";
            return 0;
        }
        else if (arg == "-n") {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: -n requires a numeric argument\n";
                return 1;
            }
            try {
                N = std::stoul(argv[++i]);
            } catch (...) {
                std::cerr << "ERROR: invalid number for -n: '" << argv[i] << "'\n";
                return 1;
            }
        }
        else if (filename.empty()) {
            filename = arg;
        }
        else {
            std::cerr << "ERROR: unexpected argument '" << arg << "'\n";
            return 1;
        }
    }

    if (filename.empty()) {
        std::cerr << "ERROR: missing matrix file\n";
        return 1;
    }

    print_first_N(filename, N);
    return 0;
}
