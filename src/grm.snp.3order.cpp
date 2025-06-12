#include <iostream>
#include <cstdint>
#include <vector> // change the vector of tuple to array and remove?
#include <array> // change the vector of tuple to array and remove?
#include <tuple>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h> //fstat
#include <chrono>  // timing
#include <string>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
#include <mkl.h>

// remove if removing available memory function that is not used
#ifdef __linux__
#include <sys/sysinfo.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

// default buffer size
static constexpr std::size_t DBUFF_ = 32;            // MiB

// compute the largest io size
static size_t get_max_io_size() {
	static size_t cap = [](){
		long page = sysconf(_SC_PAGESIZE);
		return (size_t)(INT_MAX & ~(page - 1));
	}();
	return cap;
}


// Function to get available memory in bytes
size_t get_available_memory() {
	size_t available_memory = 0;

#ifdef __linux__
	struct sysinfo info;
	if (sysinfo(&info) == 0) {
		available_memory = info.freeram;
	}
#elif defined(__APPLE__)
	uint64_t memsize;
	size_t len = sizeof(memsize);
	if (sysctlbyname("hw.memsize", &memsize, &len, NULL, 0) == 0) {
		available_memory = static_cast<size_t>(memsize);
	}
#elif defined(_WIN32)
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	if (GlobalMemoryStatusEx(&status)) {
		available_memory = status.ullAvailPhys;
	}
#endif

	return available_memory;
}


// Ensure buffer size is a multiple of sizeof(double)
size_t align_to_double(size_t size) {
	const size_t alignment = sizeof(double);  // 8 bytes
	return (size / alignment) * alignment;  // Faster than `size - (size % alignment)`
}

// Function to parse buffer size input (supports MB, GB, M, G and ensures memory alignment)
size_t parse_buffer_size(const std::string& buffer_arg) {
	if (buffer_arg.empty()) {
		return static_cast<size_t>(DBUFF_) * 1024 * 1024; // Default
	}

	size_t multiplier = 1;
	std::string num_part;

	// Convert string to lowercase & extract numeric part
	for (char c : buffer_arg) {
		if (std::isdigit(c)) {
			num_part += c;
		} else {
			c = std::tolower(c);
			if (c == 'g') {
				multiplier = 1024 * 1024 * 1024;  // Convert to bytes
				break;
			} else if (c == 'm') {
				multiplier = 1024 * 1024;
				break;
			}
		}
	}

	// Handle invalid input
	if (num_part.empty()) {
		std::cerr << "Error: Invalid buffer size input '" << buffer_arg << "'. Using default." << std::endl;
		return static_cast<size_t>(8) * 1024 * 1024; // Default: 8MB
	}

	size_t buffer_size = std::stoull(num_part) * multiplier;

	// Set limits on buffer size for performance
	size_t max_buffer_size = 80L * 1024 * 1024 * 1024; // 80GB max buffer
	buffer_size = std::min(buffer_size, max_buffer_size);

	size_t min_buffer_size = 1 * 1024 * 1024; // 1MB minimum for efficiency
	buffer_size = std::max(buffer_size, min_buffer_size);

	// Ensure it's aligned
	return align_to_double(buffer_size);
}


// Configure OpenMP and MKL threads based on system capabilities
void configure_threads(int user_threads) {
	int max_threads = omp_get_max_threads();
	int num_threads = (user_threads > 0) ? user_threads : max_threads;

	std::cout << "setting openmp threads to: " << num_threads << "/" << max_threads << std::endl;

	// Set OpenMP thread count
	omp_set_num_threads(num_threads);

}

int check_file(const std::string &path,
		std::array<uint64_t,2> &dims,
		size_t                 &file_size)
{
	// open the file
	int fd = open(path.c_str(), O_RDONLY);
	if (fd < 0) {
		std::cerr << "ERROR: cannot open \"" << path << "\": "
			<< strerror(errno) << "\n";
		return -1;
	}

	// read the dimensions
	uint64_t hdr[2];
	ssize_t r = pread(fd, hdr, sizeof(hdr), 0);
	if (r != (ssize_t)sizeof(hdr)) {
		std::cerr << "ERROR: reading header of \"" << path << "\": "
			<< (r < 0 ? strerror(errno) : "short read")
			<< "\n";
		close(fd);
		return -1;
	}

	// verify dimensions
	if (dims[0] == 0 && dims[1] == 0) {
		dims[0] = hdr[0];
		dims[1] = hdr[1];
	}
	else if (dims[0] != hdr[0] || dims[1] != hdr[1]) {
		std::cerr << "ERROR: dimension mismatch for \"" << path << "\":\n"
			<< "  expected " << dims[0] << "×" << dims[1]
			<< ", got "      << hdr[0] << "×" << hdr[1] << "\n";
		close(fd);
		return -1;
	}

	// expected file size
	uint64_t num_elements = dims[0] * dims[1];
	size_t   expect = 2 * sizeof(uint64_t)
		+ num_elements * sizeof(double);

	// init or verify file_size
	if (file_size == 0) {
		file_size = expect;
	}
	else if (file_size != expect) {
		std::cerr << "ERROR: file_size mismatch for \"" << path << "\":\n"
			<< "  expected " << expect
			<< " bytes, caller has " << file_size << "\n";
		close(fd);
		return -1;
	}

	// 5) stat the fd to confirm on-disk size
	struct stat st;
	if (fstat(fd, &st) < 0) {
		std::cerr << "ERROR: fstat(\"" << path << "\"): "
			<< strerror(errno) << "\n";
		close(fd);
		return -1;
	}
	if (size_t(st.st_size) != expect) {
		std::cerr << "ERROR: actual size for \"" << path << "\":\n"
			<< "  expected " << expect
			<< " bytes, got "      << st.st_size << "\n";
		close(fd);
		return -1;
	}

	return fd;
}


void hadamard_product(const std::string& file1, const std::string& file2, const std::string& output_file, size_t buffer_size)
{
	// use check_file to open & verify inputs, and get dims + file_size
	std::array<uint64_t,2> dims = {0,0};
	size_t file_size = 0;
	int fd1 = check_file(file1, dims, file_size);
	if (fd1 < 0) return;
	int fd2 = check_file(file2, dims, file_size);
	if (fd2 < 0) { close(fd1); return; }

	// open and size the output file, write header
	int fd_out = open(output_file.c_str(),
			O_RDWR | O_CREAT | O_TRUNC,
			0666);
	if (fd_out < 0) {
		std::cerr << "Error opening output “" << output_file << "”: "
			<< strerror(errno) << "\n";
		close(fd1); close(fd2);
		return;
	}
	if (ftruncate(fd_out, file_size) < 0) {
		std::cerr << "Error truncating output: " << strerror(errno) << "\n";
		close(fd1); close(fd2); close(fd_out);
		return;
	}
	if (pwrite(fd_out, dims.data(), sizeof(dims), 0)
			!= (ssize_t)sizeof(dims)) {
		std::cerr << "Error writing output header\n";
		close(fd1); close(fd2); close(fd_out);
		return;
	}

	// compute total elements and header size
	uint64_t num_elements = dims[0] * dims[1];
	size_t header_bytes   = sizeof(dims);
	size_t data_bytes     = file_size - header_bytes;

	//size_t per_buf = buffer_size / 3;
	size_t MAX_RW = get_max_io_size();
	size_t per_buf = std::min(buffer_size, MAX_RW);
	size_t buf_elems = std::max((size_t)1, per_buf / sizeof(double));
	size_t buf_bytes = buf_elems * sizeof(double);
	int nthreads = omp_get_max_threads();

	std::cout << "buffer for each matrix (x3): " << buf_bytes << " (" << (buf_bytes/1024/1024) << "MB)" << std::endl;
	std::cout << "using threads: " << nthreads << std::endl;
	//std::cout << "ssize_max: " << SSIZE_MAX << std::endl;

	// if use vectors
	//std::vector<double> buf1(buf_elems), buf2(buf_elems), buf_out(buf_elems);

	double *buf1, *buf2, *buf_out;
	if (posix_memalign((void**)&buf1,    64, buf_bytes) ||
			posix_memalign((void**)&buf2,    64, buf_bytes) ||
			posix_memalign((void**)&buf_out, 64, buf_bytes)) {
		std::cerr << "Error allocating buffers\n";
		// cleanup
		close(fd1); close(fd2); close(fd_out);
		return;
	}

	// stream through data region
	off_t offset = header_bytes;
	while ((size_t)offset < file_size) {
		size_t to_read = std::min(buf_bytes,
				file_size - (size_t)offset);
		size_t nelem   = to_read / sizeof(double);

		// read inputs
		ssize_t r = pread(fd1, buf1, to_read, offset);
		if (r < 0) {
			perror("pread fd1");
			break;
		}
		if ((size_t)r != to_read) {
			std::cerr << "Short read: got " << r 
				<< " of " << to_read << " bytes\n";
			break;
		}
		r = pread(fd2, buf2, to_read, offset);
		if (r < 0) {
			perror("pread fd2");
			break;
		}
		if ((size_t)r != to_read) {
			std::cerr << "Short read: got " << r 
				<< " of " << to_read << " bytes\n";
			break;
		}

		// multiply into buf_out
#pragma omp parallel for simd schedule(static) num_threads(nthreads) aligned(buf1,buf2,buf_out:64)
		for (size_t i = 0; i < nelem; ++i) {
			buf_out[i] = buf1[i] * buf2[i];
		}

		// write result
		if (pwrite(fd_out, buf_out, to_read, offset)
				!= (ssize_t)to_read)
		{
			std::cerr << "Error writing at offset "
				<< offset << "\n";
			break;
		}

		offset += to_read;
	}

	// cleanup
	close(fd1);
	close(fd2);
	close(fd_out);
}

void hadamard_product_mmap(const std::string& file1, const std::string& file2, const std::string& output_file, size_t buffer_size) {

	std::array<uint64_t,2> dims = {0,0};
	size_t  file_size = 0;

	// open input files
	// first file initializes rows/cols/file_size
	int fd1 = check_file(file1, dims, file_size);
	if (fd1 < 0) return;  

	// second file verifies same dims & size
	int fd2 = check_file(file2, dims, file_size);
	if (fd2 < 0) { close(fd1); return; }

	int fd_out = open(output_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
	// ensure output file has correct size (dimensions + matrix data)
	if (ftruncate(fd_out, file_size) == -1) {
		std::cerr << "Error extending output file! " << strerror(errno) << std::endl;
		close(fd1);
		close(fd2);
		close(fd_out);
		return;
	}

	// write dimensions at the start of the output file
	if (pwrite(fd_out, dims.data(), sizeof(dims), 0) != sizeof(dims)) {
		std::cerr << "Error writing dimensions to output file!" << std::endl;
		close(fd1);
		close(fd2);
		close(fd_out);
		return;
	}

	uint64_t num_elements = dims[0] * dims[1];

	// single big mmap instead of chunked loop
	size_t header_bytes = sizeof(dims);
	size_t data_bytes   = file_size - header_bytes;
	int    nthreads     = omp_get_max_threads();
	uint64_t N          = num_elements;  // dims[0]*dims[1]

	std::cout << "using threads: " << nthreads << std::endl;


	// align mmap offset down to a page boundary
	long   page_size    = sysconf(_SC_PAGESIZE);
	off_t  map_off      = header_bytes & ~(page_size - 1);   // → 0
	size_t map_diff     = header_bytes - map_off;           // → 16
	size_t map_len      = data_bytes + map_diff;            // cover header + data

	void* m1 = mmap(nullptr, map_len,    PROT_READ,  MAP_PRIVATE, fd1,    map_off);
	void* m2 = mmap(nullptr, map_len,    PROT_READ,  MAP_PRIVATE, fd2,    map_off);
	void* mo = mmap(nullptr, map_len,    PROT_WRITE, MAP_SHARED,  fd_out, map_off);
	if (m1==MAP_FAILED || m2==MAP_FAILED || mo==MAP_FAILED) {
		perror("mmap");
		close(fd1); close(fd2); close(fd_out);
		return;
	}

	// bump each base pointer past the 16-byte header
	double* A = reinterpret_cast<double*>((char*)m1 + map_diff);
	double* B = reinterpret_cast<double*>((char*)m2 + map_diff);
	double* O = reinterpret_cast<double*>((char*)mo + map_diff);

#pragma omp parallel for simd schedule(static) num_threads(nthreads)
	for (uint64_t i = 0; i < N; ++i) {
		O[i] = A[i] * B[i];
	}

	// cleanup
	munmap(m1, data_bytes);
	munmap(m2, data_bytes);
	munmap(mo, data_bytes);
	close(fd1);
	close(fd2);
	close(fd_out);

}

void hadamard_product_mmap_old(const std::string& file1, const std::string& file2, const std::string& output_file, size_t buffer_size) {
	// Open input files

	std::array<uint64_t,2> dims = {0,0};
	size_t  file_size = 0;

	// first file initializes rows/cols/file_size
	int fd1 = check_file(file1, dims, file_size);
	if (fd1 < 0) return;  

	// second file verifies same dims & size
	int fd2 = check_file(file2, dims, file_size);
	if (fd2 < 0) { close(fd1); return; }

	int fd_out = open(output_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
	// Ensure output file has correct size (dimensions + matrix data)
	if (ftruncate(fd_out, file_size) == -1) {
		std::cerr << "Error extending output file! " << strerror(errno) << std::endl;
		close(fd1);
		close(fd2);
		close(fd_out);
		return;
	}

	// Write dimensions at the start of the output file
	if (pwrite(fd_out, dims.data(), sizeof(dims), 0) != sizeof(dims)) {
		std::cerr << "Error writing dimensions to output file!" << std::endl;
		close(fd1);
		close(fd2);
		close(fd_out);
		return;
	}

	uint64_t num_elements = dims[0] * dims[1];

	// Setup memory alignment
	size_t page_size = sysconf(_SC_PAGE_SIZE);
	size_t data_offset = sizeof(dims);  // Ensure matrix data starts after the 16-byte dimensions

	// Process file in chunks
	// Ensure at least 1MB per chunk and not bigger than file
	size_t chunk_size = std::max(std::min(buffer_size / sizeof(double), static_cast<size_t>(num_elements)), static_cast<size_t>(1024 * 1024));
	int nthreads = omp_get_max_threads();
	size_t count_doubles = buffer_size / sizeof(double);
	size_t per_thread  = count_doubles / nthreads;
	size_t doubles_per_chunk = per_thread * nthreads;
	size_t chunk_bytes = doubles_per_chunk * sizeof(double);

	// Debugging: Print file sizes and alignment details
	std::cerr << "File1: fd=" << fd1 << ", size=" << lseek(fd1, 0, SEEK_END) << " bytes\n";
	std::cerr << "File2: fd=" << fd2 << ", size=" << lseek(fd2, 0, SEEK_END) << " bytes\n";
	std::cerr << "Output File: fd=" << fd_out << ", size=" << lseek(fd_out, 0, SEEK_END) << " bytes\n";
	std::cerr << "Trying to mmap with:\nchunk_bytes=" << chunk_bytes
		<< "\n num_elements=" << num_elements
		<< "\n count_doubles=" << count_doubles
		<< "\n per_thread=" << per_thread
		<< "\n doubles_per_chunk=" << doubles_per_chunk
		<< "\n chunk_bytes=" << chunk_bytes
		<< "\n nthreads=" << nthreads
		<< "\n page_size=" << page_size << std::endl;


	off_t offset = data_offset;  // Start after dimensions


	while (offset < file_size) {
		// figure out how many bytes to map this iteration
		size_t remaining    = file_size - offset;
		//size_t want_bytes   = chunk_size * sizeof(double);
		//size_t map_size     = (remaining < want_bytes ? remaining : want_bytes);
		size_t map_size = (remaining < chunk_bytes ? remaining : chunk_bytes);

		// align the mmap() start to the system page boundary
		off_t map_start     = (offset / page_size) * page_size;
		size_t page_diff    = offset - map_start;     // bytes from map_start to our real data
		size_t map_len      = page_diff + map_size;   // map enough to cover the chunk

		// maps for each file
		void* m1 = mmap(NULL, map_len, PROT_READ,  MAP_PRIVATE, fd1, map_start);
		void* m2 = mmap(NULL, map_len, PROT_READ,  MAP_PRIVATE, fd2, map_start);
		void* mo = mmap(NULL, map_len, PROT_WRITE, MAP_SHARED,  fd_out, map_start);
		if (m1 == MAP_FAILED || m2 == MAP_FAILED || mo == MAP_FAILED) {
			std::cerr << "mmap failed at offset " << offset << "\n";
			// … clean up fds, unmap any successes, then exit …
		}

		//compute pointers exactly at the first double of this chunk
		double* data1_ptr = reinterpret_cast<double*>(
				reinterpret_cast<char*>(m1) + page_diff);
		double* data2_ptr = reinterpret_cast<double*>(
				reinterpret_cast<char*>(m2) + page_diff);
		double* out_ptr   = reinterpret_cast<double*>(
				reinterpret_cast<char*>(mo) + page_diff);

		//run the Hadamard product
		size_t n = map_size / sizeof(double);
		//#pragma omp parallel for simd schedule(guided, chunk_per_thread)
		//#pragma omp parallel for simd schedule(guided)
#pragma omp parallel for simd schedule(static) num_threads(nthreads)
		for (size_t i = 0; i < n; ++i) {
			out_ptr[i] = data1_ptr[i] * data2_ptr[i];
		}

		// 6) unmap exactly what we mapped
		munmap(m1, map_len);
		munmap(m2, map_len);
		munmap(mo, map_len);

		// 7) advance to the next chunk
		offset += map_size;
	}

	// Close files
	close(fd1);
	close(fd2);
	close(fd_out);
}

static void print_usage(const char* prog_name) {
	std::cout
		<< "Usage: " << prog_name << " [options] <base_name>\n"
		<< "\n"
		<< "Options:\n"
		<< "  -h, --help           Show this help message and exit\n"
		<< "  -b, --buffer SIZE    Buffer size (e.g. 8M, 1G). Default: " << DBUFF_
		<< "M\n                       will use 3*buffer memroy for the 2 input and 1 output matrices\n"
		<< "  -t, --threads N      Number of threads (default: system max)\n"
		<< "\n"
		<< "Operation filters (if none given, all are run):\n"
		<< "  --aaa                Compute only “.g.AAA” output\n"
		<< "  --add                Compute only “.g.ADD” output\n"
		<< "  --aad                Compute only “.g.AAD” output\n"
		<< "  --ddd                Compute only “.g.DDD” output\n";
}

int main(int argc, char* argv[]) {
	auto t0 = std::chrono::high_resolution_clock::now();
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <base_name> [buffer_size (MB/GB)] [num_threads]" << std::endl;
		return 1;
	}

	const char* prog = argv[0];

	if (argc < 2) {
		print_usage(prog);
		return 1;
	}

	std::string base_name;
	size_t buffer_size = -1;
	int num_threads   = -1;
	std::vector<int> selected_ops;

	// Map option name → index in operations array
	const std::vector<std::pair<std::string,int>> op_map = {
		{"aaa", 0},
		{"add", 1},
		{"aad", 2},
		{"ddd", 3}
	};

	// Parse arguments
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "-h" || arg == "--help") {
			print_usage(prog);
			return 0;
		}
		else if (arg == "-b" || arg == "--buffer") {
			if (i+1 >= argc) {
				std::cerr << "ERROR: “" << arg << "” requires an argument\n";
				return 1;
			}
			buffer_size = parse_buffer_size(argv[++i]);
		}
		else if (arg == "-t" || arg == "--threads") {
			if (i+1 >= argc) {
				std::cerr << "ERROR: “" << arg << "” requires an argument\n";
				return 1;
			}
			num_threads = std::stoi(argv[++i]);
		}
		else if (arg.rfind("--", 0) == 0) {
			std::string name = arg.substr(2);
			bool found = false;
			for (auto &p : op_map) {
				if (p.first == name) {
					selected_ops.push_back(p.second);
					found = true;
					break;
				}
			}
			if (!found) {
				std::cerr << "ERROR: unknown option “" << arg << "”\n";
				return 1;
			}
		}
		else if (base_name.empty()) {
			base_name = arg;
		}
		else {
			std::cerr << "ERROR: unexpected argument “" << arg << "”\n";
			return 1;
		}
	}

	if (base_name.empty()) {
		std::cerr << "ERROR: missing <base_name>\n";
		print_usage(prog);
		return 1;
	}

	// If no operation flags, run all
	if (selected_ops.empty()) {
		selected_ops = {0, 1, 2, 3};
	}

    long L3 = sysconf(_SC_LEVEL3_CACHE_SIZE);
	std::cout << "Your L3: " << L3 << " (" << (L3/1024/1024) << "MB)" << std::endl;


	configure_threads(num_threads);
    num_threads = omp_get_max_threads();
    // new default: based on L3 cache
    if (buffer_size == -1) {
        buffer_size = (L3 * num_threads) / 3;
    }
	std::cout << "Using buffer size: " << buffer_size << " (" << (buffer_size/1024/1024) << "MB)" << std::endl;


	// Define input and output file names
	std::array<std::tuple<std::string, std::string, std::string>, 4> operations = {{
		{base_name + ".g.A", base_name + ".g.AA", base_name + ".g.AAA"},
			{base_name + ".g.A", base_name + ".g.DD", base_name + ".g.ADD"},
			{base_name + ".g.D", base_name + ".g.AA", base_name + ".g.AAD"},
			{base_name + ".g.D", base_name + ".g.DD", base_name + ".g.DDD"}
	}};


	// Perform Hadamard products with timing
	for (int idx : selected_ops) {
		const auto& op = operations[idx];
		std::string in1 = std::get<0>(op);
		std::string in2 = std::get<1>(op);
		std::string out = std::get<2>(op);

		std::cout << "\nProcessing “" << in1
			<< " * " << in2
			<< " -> " << out << "”\n";

		auto t1 = std::chrono::high_resolution_clock::now();
		hadamard_product(in1, in2, out, buffer_size);
		auto t2 = std::chrono::high_resolution_clock::now();

		std::cout << "done in "
			<< std::chrono::duration<double>(t2-t1).count()
			<< " sec\n";
	}
	auto t3 = std::chrono::high_resolution_clock::now();

	std::cout << "\nTotal elapsed: "
		<< std::chrono::duration<double>(t3-t0).count()
		<< " sec\n";
	return 0;
}
