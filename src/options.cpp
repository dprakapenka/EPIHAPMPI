#include "options.h"
#include "read.h"  // for str_split()

#include <getopt.h>     // long options
#include <iostream>     // For std::cout, std::endl
#include <vector>       // For std::vector
#include <string>       // For std::string
#include <unistd.h>     // For getopt, optarg

#include <iomanip>      // for std::setw, std::left

void parse_options(int argc, char *argv[], int myrank_mpi, Options &options) {

    // get num mkl threads
    options.num_threads = mkl_get_max_threads();

    //options
    enum {
        OPT_VAA,
        OPT_VAD,
        OPT_VDD,
        OPT_VAAA,
        OPT_VAAD,
        OPT_VADD,
        OPT_VDDD,
        OPT_MISSING,
        OPT_TOL,
        OPT_HTOL
    };

    // long options
    static struct option long_opts[] = {
        { "input",     required_argument, nullptr, 'i' },
        { "missing",   required_argument, nullptr, 'm' },
        { "pheno",     required_argument, nullptr, 'p' },
        { "blocks",    required_argument, nullptr, 'n' },
        { "iterations",     required_argument, nullptr, 'I' },
        { "tol",       required_argument, nullptr, 't' },
        { "htol",      required_argument, nullptr, 'b' },
        { "ai_start",     required_argument, nullptr, 's' },
        { "va",     required_argument, nullptr, 'a' },
        { "vd",     required_argument, nullptr, 'd' },
        { "ve",     required_argument, nullptr, 'e' },
        { "vah",     required_argument, nullptr, 'h' },
        { "vaa",     required_argument, nullptr, OPT_VAA },
        { "vad",     required_argument, nullptr, OPT_VAD },
        { "vdd",     required_argument, nullptr, OPT_VDD },
        { "vaaa",     required_argument, nullptr, OPT_VAAA },
        { "vaad",     required_argument, nullptr, OPT_VAAD },
        { "vadd",     required_argument, nullptr, OPT_VADD },
        { "vddd",     required_argument, nullptr, OPT_VDDD },
        { "out",       required_argument, nullptr, 'o' },
        { "flag",      no_argument,       nullptr, 'f' },
        { nullptr,     0,                 nullptr,  0  }
    };

    int opt;
    std::vector<std::string> pheno_name;
    if (myrank_mpi == 0) printf("Usage: ./reml CHANGE\n");

    //while((opt = getopt(argc, argv, ":i:p:r:c:n:l:t:b:s:a:d:o:h:e:f")) != -1) 
    while((opt = getopt_long(argc, argv, ":i:p:r:c:n:m:I:t:b:s:a:d:o:h:e:f", long_opts, nullptr)) != -1) 
    { 
        std::cout << "processing option " << opt << " --- " << optarg << std::endl;
        switch(opt) 
        { 

            case 'n': // block size
                if (myrank_mpi == 0) printf("set block size: %i\n", atoi(optarg)); 
                options.mb = options.nb = atoi(optarg);
                break;
            case 'r': // processor rows
                if (myrank_mpi == 0) printf("override num proc rows: %i\n", atoi(optarg)); 
                options.nprow = atoi(optarg);
                break; 
            case 'c': // processor cols
                if (myrank_mpi == 0) printf("override num proc cols: %i\n", atoi(optarg)); 
                options.npcol = atoi(optarg);
                break; 
            case 'i': //input data (load) project name
                if (myrank_mpi == 0) printf("proj_name: %s\n", optarg); 
                options.load_name = optarg;
                break; 
            case 'm': //missing phenotype number (default: -9999111)
                if (myrank_mpi == 0) printf("missing: %s\n", optarg); 
                options.missing_phenotype = atof(optarg);
                break; 
            case 'o': //output path/project name
                if (myrank_mpi == 0) printf("out_name: %s\n", optarg); 
                options.save_name = optarg;
                break; 
            case 'p': //pheno file and column "2,phenotype.txt"
                if (myrank_mpi == 0) printf("pheno trait_col,file: %s\n", optarg); 
                pheno_name = str_split(",", optarg);
                if (pheno_name.size() == 2) {
                    try {
                        options.trait_col = std::stoi(pheno_name[0]);
                        options.pheno_file = pheno_name[1];
                        if (myrank_mpi == 0) { 
                            std::cout << "phenotype file: " << options.pheno_file << std::endl;
                            std::cout << "trait column: " << options.trait_col << std::endl;
                        }
                    } catch (const std::exception &e) {
                        if (myrank_mpi == 0) {
                            std::cerr << "Error converting trait column to integer: " << e.what() << std::endl;
                        }
                    }
                } else {
                    std::cerr << "Invalid phenotype input format: " << optarg << std::endl;
                }
                break; 
            case 'f': 
                if (myrank_mpi == 0) printf("f flag set: %c\n", opt); 
                break;
            case 'I': //limit of iterations
                options.iter_n = atoi(optarg);
                if (myrank_mpi == 0) std::cout << "set iterations: " << options.iter_n << std::endl;
                break;
            case 'a': //additive variance
                options.var_snp_a = atof(optarg);
                options.variances["a"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set a variance: " << options.var_snp_a << std::endl;
                break;
            case 'd': //dominance variance
                options.var_snp_d = atof(optarg);
                options.variances["d"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set d variance: " << options.var_snp_d << std::endl;
                break;
            case 'h': //haplotype additive variance
                options.var_hap_a = atof(optarg);
                options.variances["h"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set ah variance: " << options.var_hap_a << std::endl;
                break;
            case OPT_VAA: //additive x additive variance
                options.variances["aa"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set aa variance: " << optarg << std::endl;
                break;
            case OPT_VAD: //additive x additive variance
                options.variances["ad"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set ad variance: " << optarg << std::endl;
                break;
            case OPT_VDD: //additive x additive variance
                options.variances["dd"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set dd variance: " << optarg << std::endl;
                break;
            case OPT_VAAA: //additive x additive variance
                options.variances["aaa"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set aaa variance: " << optarg << std::endl;
                break;
            case OPT_VAAD: //additive x additive variance
                options.variances["aad"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set aad variance: " << optarg << std::endl;
                break;
            case OPT_VADD: //additive x additive variance
                options.variances["add"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set add variance: " << optarg << std::endl;
                break;
            case OPT_VDDD: //additive x additive variance
                options.variances["ddd"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set ddd variance: " << optarg << std::endl;
                break;
            case 'e': //residual variance
                options.var_e = atof(optarg);
                options.variances["e"] = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set residual: " << options.var_e << std::endl;
                break;
            case 't': //tolerance
                options.tolerance = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set variance tolerance: " << options.tolerance << std::endl;
                break;
            case 'b': //heritability tolerance (negative to skip)
                options.htolerance = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set heritability tolerance: " << options.htolerance << std::endl;
                break;
            case 's': //ai start at iteration (negative to skip)
                options.ai_start = atof(optarg);
                if (myrank_mpi == 0) std::cout << "set ai_start at iteration: " << options.ai_start << std::endl;
                break;
            case ':':
                if (myrank_mpi == 0)
                    std::cerr << "Missing value for option: " << argv[optind - 1] << std::endl;
                break;
            case '?': 
                if (myrank_mpi == 0) printf("unknown option: %c\n", optopt);
                break; 
        } 
    } 

    // optind is for the extra arguments
    // which are not parsed
    for(; optind < argc; optind++) {     
        if (options.load_name.empty()) {
            options.load_name = argv[optind];
        } else {
            if (myrank_mpi == 0) printf("extra arguments: %s\n", argv[optind]); 
        }
    }
}
