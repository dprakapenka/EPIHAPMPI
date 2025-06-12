#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "mkl.h"

// MPI root process rank
static const int MPI_ROOT_PROC_ = 0;

// Length of BLACS descriptor arrays
static const MKL_INT DESC_LEN_ = 9;

// BLACS matrix operation flags
static       char CHAR_RIGHT_            = 'R';
static       char CHAR_LEFT_             = 'L';
static       char CHAR_LOWER_            = 'L';
static       char CHAR_UPPER_            = 'U';
static       char CHAR_TRANS_            = 'T';
static       char CHAR_NOTRANS_          = 'N';
static       char CHAR_LAYOUT_           = 'R';
static       char CHAR_VCHAR_            = 'V';
static       char CHAR_ACHAR_            = 'A';

// Integer scalars for ScaLAPACK
static            MKL_INT IZERO_         = 0;
static            MKL_INT IONE_          = 1;
static const      MKL_INT INEGONE_       = -1;

// Double scalars for ScaLAPACK
static double DZERO_               = 0.0;
static double       DONE_                = 1.0;
static const double DNEGONE_             = -1.0;

// Miscellaneous constants
static const int    NDIMS_               = 2;
static const int    SKIP_BYTES_          = 2 * sizeof(MKL_INT);
static const int    BUF_SIZE_            = 1024 * 1024;
static const        MKL_INT LWORK_QUERY_ = -1;

#endif // CONSTANTS_H
