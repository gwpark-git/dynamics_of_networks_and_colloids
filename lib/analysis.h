#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <iostream>
#include "matrix.h"
#include "omp.h"

namespace POST_PROCESSING
{
  double
  compute_autocorrelation(MATRIX& time_series_data, const MKL_LONG index,
                          MATRIX& result, const MKL_LONG index_result);

  double
  compute_autocorrelation_OMP(MATRIX& time_series_data, const MKL_LONG index,
                              const MKL_LONG N_THREADS,
                              MATRIX& result, const MKL_LONG index_result);
    
}

#endif // ANALYSIS_H
