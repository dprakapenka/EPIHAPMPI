#include "test_read_binmat_shared.cpp"

const std::string SMALL_MATRIX_FILE = "test/small_matrix.bin";

int main() {
    std::vector<std::pair<MKL_INT, MKL_INT>> test_sizes = {
        {2, 2}, {3, 3}, {5, 5}, {10, 10}, {4, 6}, {7, 8}
    };

    for (const auto& size : test_sizes) {
        MKL_INT rows = size.first, cols = size.second;
        std::cout << "\n Running SMALL matrix test for size: (" << rows << ", " << cols << ")...\n";

        std::vector<double> original_matrix(rows * cols);
        generate_test_matrix_file(SMALL_MATRIX_FILE, rows, cols, original_matrix);
        test_read_binmat_size(SMALL_MATRIX_FILE, rows, cols);
        validate_matrix(SMALL_MATRIX_FILE, rows, cols, original_matrix);
    }

    std::cout << " All SMALL matrix tests passed!\n";
    MPI_Finalize();
    return 0;
}
