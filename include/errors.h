#ifndef ERRORS_H
#define ERRORS_H

#include <string>
#include "mkl.h"

/// Print info on allocation failure, exit
void allocation_error(const std::string& name, MKL_INT mem_req, MKL_INT ictxt);

#endif // ERRORS_H
