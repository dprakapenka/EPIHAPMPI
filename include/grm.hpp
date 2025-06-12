#include <iomanip>
#include <unordered_map>
#include <array>
#include <vector>

// define hash
namespace std{
    template<>
        struct hash<array<int, 2> >
        {
            typedef array<int, 2> argument_type;
            typedef size_t result_type;

            result_type operator()(const argument_type& a)const
            {
                hash<int> hasher;
                result_type h = 0;
                for (result_type i = 0; i < 2; ++i){
                    h = h * 31 + hasher(a[i]);
                }
                return h;
            }
        };

}

typedef struct _IDXCOUNT{int idx;double count;} idxcount;
typedef std::unordered_map<int,idxcount> allele_list;
typedef std::vector<allele_list> hap_block_allele; 
typedef std::array<int,2> geno_code;
typedef std::unordered_map<geno_code,idxcount> het_geno_list;
typedef std::vector<het_geno_list> hap_block_het_geno;
hap_block_het_geno  hb_het_geno;
hap_block_allele    hb_allele;
std::vector<int>   hap_max_key;

std::vector<std::string>   var_hap_name;
