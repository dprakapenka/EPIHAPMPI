#ifndef GBLUP_H
#define GBLUP_H

#include "options.h"
#include "mpi.h"
#include "mkl_scalapack.h"
#include <vector>
#include <string>
#include <map>
#include <set>

/**  
 * Calculate GBLUP for each variance component (A, D, …).  
 *
 * @param options            Final variances are read from options.variances  
 * @param ictxt              BLACS context  
 * @param mkl_num_ind        Global # individuals (Z.rows())  
 * @param ve_i               # of variance components (excluding residual)  
 * @param Z,descZ            Distributed Z matrix  
 * @param PY,descPY          Distributed P*Y vector (N×1)  
 * @param Glocal             Map[name]→local block of G_name (N×N)  
 * @param descG              Descriptor for each G matrix  
 * @param u_effects_global   OUT: on root, u_effects_global[name][i] = GEBV_i for component ‘name’  
 * @param individual_ids     Global ordering of IDs (size = mkl_num_ind)  
 * @param validation_ids     Set of IDs in the validation set  
 */
void calculate_gblup(
    const Options                  &options,
    MKL_INT                          ictxt,
    MKL_INT                          mkl_num_ind,
    double                          *Z,       const MKL_INT *descZ,
    double                          *PY,      const MKL_INT *descY,
    std::map<std::string,double*>   &Glocal,  const MKL_INT *descG,
    std::map<std::string, std::vector<double>> &u_effects_global
);

#endif // GBLUP_H
