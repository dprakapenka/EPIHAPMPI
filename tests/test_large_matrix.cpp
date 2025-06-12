#include "test_read_binmat_shared.cpp"

const std::string LARGE_MATRIX_FILE = "test/large_matrix.bin";

int main() {
    MKL_INT rows = 46341, cols = 46341;
    std::cout << "\n Running LARGE matrix test for size: (" << rows << ", " << cols << ")...\n";

    std::vector<double> original_matrix(rows * cols);
    generate_test_matrix_file(LARGE_MATRIX_FILE, rows, cols, original_matrix);
    test_read_binmat_size(LARGE_MATRIX_FILE, rows, cols);
    validate_matrix(LARGE_MATRIX_FILE, rows, cols, original_matrix);

    std::cout << " Large matrix test completed successfully!\n";
    MPI_Finalize();
    return 0;
}
