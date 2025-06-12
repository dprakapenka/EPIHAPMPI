#!/bin/bash
# by Denis Prakapenka
# creates parameter file and sbatch files with some defaults
# and user options
# usage: ./make-params-sbatch.sh -t fatpct -hs snp4 -p ../pheno_propct_fatpct.txt -g ../geno -m ../../map.txt -a 7 -d 7 -ha 7 --save
#
# TODO:
#   -add factors
#   -add covars

# parameter defaults
#bin_file=greml_qm
#bin_file=gvchap
bin_file=/project/dairyxbreed/bin/gvchap
hap_prefix=../hap.geno

var_e=3

iterations=40000

tol="1.0E-08"
tol_h="1.0E-06"

missing_p=-9999
ai_reml=3

# log direcotory
slogs="slogs/"
# sbatch run files
sbatch_folder="sbatch/"
# save directory
save_folder="saved/"
# parameters directory
params_folder="params/"

# sbatch defaults
num_threads=20
# memory in gigabytes
memory=256
# time in hours
runtime=336
email="email@domain.edu"
account="account_name"
partition="partition_name"

positional_args=()

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -t|--trait)
        trait=$2
        shift # past argument
        shift # past value
        ;;
        -p|--pheno_file)
        pheno_file=$2
        shift # past argument
        shift # past value
        ;;
        -g|--geno)
        geno_path=$2
        shift # past argument
        shift # past value
        ;;
        -m|--map)
        map_file=$2
        shift # past argument
        shift # past value
        ;;
        -hs|--hap_size)
        hap_size=$2
        shift # past argument
        shift # past value
        ;;
        -hp|--hap_prefix)
        hap_prefix=$2
        shift # past argument
        shift # past value
        ;;
        --hap)
        hap_path=$2
        shift # past argument
        shift # past value
        ;;
        --save)
        g_save=true
        shift # past argument
        ;;
        --load)
        g_load=true
        shift # past argument
        ;;
        --load-mpi)
        g_load_mpi=true
        shift # past argument
        ;;
        -t|--tolerance)
        tol=$2
        shift # past argument
        shift # past value
        ;;
        -th|--tolerance_her)
        tol_h=$2
        shift # past argument
        shift # past value
        ;;
        --no_skip_rel)
        noskiprel=true
        shift # past argument
        ;;
        -i|--iter)
        iterations=$2
        shift # past argument
        shift # past value
        ;;
        -n|--num_inds)
        num_inds=$2
        shift # past argument
        shift # past value
        ;;
        -a)
        var_a=$2
        shift # past argument
        shift # past value
        ;;
        -d)
        var_d=$2
        shift # past argument
        shift # past value
        ;;
        -e)
        var_e=$2
        shift # past argument
        shift # past value
        ;;
        -ha)
        var_ha=$2
        shift # past argument
        shift # past value
        ;;
        -b|--bin)
        bin_file=$2
        shift # past argument
        shift # past value
        ;;
        --threads)
        num_threads=$2
        shift # past argument
        shift # past value
        ;;
        --mem|--memory)
        memory=$2
        shift # past argument
        shift # past value
        ;;
        --time|--runtime)
        runtime=$2
        shift # past argument
        shift # past value
        ;;
        -c|--trait_col)
        trait_col=$2
        shift # past argument
        shift # past value
        ;;
        --parameter)
        new_para_file=$2
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        positional_args+=($1) # save it in an array for later
        shift # past argument
        ;;
    esac
done

set -- "${positional_args[@]}" # restore positional parameters
if (( ${#positional_args[@]} )); then
    echo
    echo "*** unknown option(s): ${positional_args[*]}"
    echo
fi

identifier=${trait}

if [ -z "$hap_path" ]; then
    if [ -n "$hap_size" ]; then
        hap_path="${hap_prefix}/${trait}/${hap_size}"
    fi
fi

if [ -n "$hap_size" ]; then
    identifier=${identifier}.${hap_size}
fi

save_name=$identifier

if [ -n "$geno_path" ]; then
    if [ -d ${geno_path} ]; then
        geno_path="${geno_path}/${trait}"
        geno_files=( $geno_path/* )
    else
        echo "*** genotype path doesn't exist ${geno_path}"
    fi
fi

if [ -n "${var_a}" ]; then
    identifier=${identifier}.sa
fi
if [ -n "${var_d}" ]; then
    identifier=${identifier}.sd
fi

if [ -n "$hap_path" ]; then
    if [ -d "$hap_path" ]; then
        identifier=${identifier}.ha
        hap_files=( $hap_path/* )
        if [ -z "$var_ha" ]; then
            echo "*** path to haplotype genotypes is set but no starting value (ha) specified"
            echo "*** setting to 3"
            var_ha=3
        fi
    else
        echo "*** haplotype genotype path doesn't exist ${hap_path}"
    fi
fi

#TODO: add flag for trait name and trait is the actual (validation) column name
# check phenotype input file
if [ -n "$pheno_file" ]; then
    if [ -f "$pheno_file" ]; then
        if [ -z "$trait_col" ]; then
            read -r pheno_header<${pheno_file}
            #echo ${pheno_header}
            i=0
            for col in $pheno_header; do 
               ((i++)); [[ "$col" = "$trait" ]] && trait_col=$i; 
            done

            if [ -z "$trait_col" ]; then
                echo "couldn't find trait in phenotype file"
                echo "setting trait column to 1"
                trait_col=1
            fi
        fi
    else
        echo "*** can't find phenotype file ${pheno_file}"
        exit 1
    fi

else
    echo "*** no phenotype file given"
    exit 1
fi

if [ -n "$trait_col" ]; then
    identifier=${identifier}.${trait_col}
fi


echo "IDENTIFIER: ${identifier}"
if [ -z "$new_para_file" ]; then
    new_para_file=params.${identifier}.txt
fi

# put new parameter file into a folder
mkdir -p ${params_folder}
new_para_file=$params_folder$new_para_file

echo "parameter file: $new_para_file"
# start writing paramter file
echo "##--- ${identifier} --- $(date) ---##" > ${new_para_file}

# set phenotype input file and trait column
echo "phenotype ${pheno_file}" >> ${new_para_file}
echo "traits_col ${trait_col}" >> ${new_para_file}

# set save or load
#save_name=$(echo $identifier | cut -d'.' -f2-)
#save_name=$identifier
#save_name=${save_name}.${trait_col}
if [ -n "$g_save" ] && [ -z "$g_load" ]; then
    mkdir -p ${save_folder}
    echo "save_g ${save_folder}${save_name}" >> ${new_para_file}
fi
if [ -n "$g_load" ] && [ -z "$g_save" ]; then
    mkdir -p ${save_folder}
    echo "load_g ${save_folder}${save_name}" >> ${new_para_file}
fi
if [ -n "$g_load_mpi" ] && [ -z "$g_load" ] && [ -z "$g_save" ]; then
    mkdir -p ${save_folder}
    echo "load_g_mpi ${save_folder}${save_name}" >> ${new_para_file}
fi

# set number of iterations
echo "num_iter ${iterations}" >> ${new_para_file}

#echo "$geno_path $hap_path $var_a $var_d $var_ha $var_hd"

# set starting value for snps
if [ -n "$geno_path" ]; then
    if [ -n "$var_a" ] || [ -n "$var_d" ]; then
        if [ -n "$var_a" ]; then
            echo "var_snp_a ${var_a}" >> ${new_para_file}
        fi
        if [ -n "$var_d" ]; then
            echo "var_snp_d ${var_d}" >> ${new_para_file}
        fi
    else
        echo "*** genotype path specified but no starting values"
        exit 1
    fi
fi
echo "var_snp_e ${var_e}" >> ${new_para_file}

# set starting value for haps
if [ -n "${hap_path}" ]; then
    if [ -n "${var_ha}" ]; then
        echo "var_hap_a ${var_ha}" >> ${new_para_file}
    fi
    if [ -n "${var_hd}" ]; then
        echo "var_hap_d ${var_hd}" >> ${new_para_file}
    fi
fi

# set map file
if [ -n "$map_file" ]; then
    echo "map_file ${map_file}" >> ${new_para_file}
fi

# a pain to format, remove
#trait_col=$(( sed -n $'1s/\t/\\\n/gp' ${pheno_file} | grep -nx ${trait} ))

# write number of individuals, get from geno or hap if none specified
if [ -z ${num_inds} ]; then
    if [ -n "${geno_files}" ]; then
        num_inds=`expr $(wc -l ${geno_files} | awk '{print $1}') - 1`
    elif [ -n "${hap_files}" ]; then
        num_inds=`expr $(wc -l ${hap_files} | awk '{print $1}') - 1`
    else
        echo "*** can't get number of individuals, exiting"
        exit 1
    fi
fi
echo "num_ind ${num_inds}" >> ${new_para_file}

# set tolerances
echo "tolerance ${tol}" >> ${new_para_file}
echo "tolerance_her ${tol_h}" >> ${new_para_file}

# set threads
echo "num_threads ${num_threads}" >> ${new_para_file}

# set missing phenotype value
echo "missing_phen_val ${missing_p}" >> ${new_para_file}

# set ai_reml
echo "use_ai_reml ${ai_reml}" >> ${new_para_file}

# add output files
mkdir -p ${trait}
echo "out_greml ./${trait}/greml-${identifier}.txt" >> ${new_para_file}
echo "out_gblup ./${trait}/gblup-${identifier}.txt" >> ${new_para_file}

echo "output_fixed_effect ./${trait}/fixed-effect-${identifier}.txt" >> ${new_para_file}
if [ -n "${var_a}" ] || [ -n "${var_d}" ]; then
    echo "output_mrk_effect ./${trait}/mrk-effect-${identifier}.txt" >> ${new_para_file}
fi
if [ -n "${var_ha}" ] || [ -n "${var_hd}" ]; then
    echo "output_hap_effect ./${trait}/hap-effect-${identifier}.txt" >> ${new_para_file}
fi

# fill in genotypes
if [ -n "$geno_path" ]; then
    count-snps.py -i ${geno_path}/* >> ${new_para_file}
fi

# fill in haplotypes
if [ -n "$hap_path" ]; then
    count-haps.py -i ${hap_path}/* >> ${new_para_file}
fi


# start sbatch file
# rename
#sbatch_file="${new_para_file//params/run}"
#sbatch_file="${sbatch_folder}/${sbatch_file//txt/sbatch}"
#echo "sbatch file: $sbatch_file"

mkdir -p ${sbatch_folder}
sbatch_file="./${sbatch_folder}/run.${identifier}.sbatch"
echo "sbatch file: $sbatch_file"

start_line='#SBATCH '
name_bin_file="$(basename bin_file)"
job_name="${start_line}--job-name ${identifier}-${name_bin_file}"
echo "#!/bin/bash" > ${sbatch_file}
echo ${job_name} >> ${sbatch_file}
echo "${start_line} --account ${account}" >> ${sbatch_file}
echo "${start_line} -p ${partition}" >> ${sbatch_file}
echo "${start_line} -N 1" >> ${sbatch_file}
echo "${start_line} --mem=${memory}G" >> ${sbatch_file}
echo "${start_line} -n ${num_threads}" >> ${sbatch_file}
echo "${start_line} -t ${runtime}":00:00 >> ${sbatch_file}
echo "${start_line} --mail-user=${email}" >> ${sbatch_file}
#echo "${start_line} --mail-type=BEGIN,END,FAIL" >> ${sbatch_file}
echo "${start_line} --mail-type=FAIL" >> ${sbatch_file}
mkdir -p ${slogs}
#echo "${start_line} -o \"${slogs}o.${identifier}.%j.%N.txt\"" >> ${sbatch_file}
#echo "${start_line} -o \"${slogs}e.${identifier}.%j.%N.txt\"" >> ${sbatch_file}
echo "${start_line} -o \"../${slogs}o.${identifier}.%j.txt\"" >> ${sbatch_file}
echo "${start_line} -o \"../${slogs}e.${identifier}.%j.txt\"" >> ${sbatch_file}
echo "date" >> ${sbatch_file}
echo "cd .." >> ${sbatch_file}
echo "${bin_file} ${new_para_file} > ${trait}/${identifier}.log.txt" >> ${sbatch_file}
echo "date" >> ${sbatch_file}

