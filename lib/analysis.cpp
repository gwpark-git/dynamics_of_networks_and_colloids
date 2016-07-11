#include "analysis.h"

using namespace std;

double
POST_PROCESSING::
compute_autocorrelation(MATRIX& time_series_data, const MKL_LONG index,
                        MATRIX& result, const MKL_LONG index_result)
{
  MKL_LONG Nt_corr = time_series_data.rows/2;
  // MATRIX result(Nt_corr, 1, 0.);
  for(MKL_LONG t=0; t<Nt_corr; t++)
    {
      for(MKL_LONG s=0; s<Nt_corr -1; s++)
        {
          result(t) += time_series_data(s, index)*time_series_data(s + t, index);
        }
      result(t) /= (double)(Nt_corr - 1);
    }
  return 0.;
}

double
POST_PROCESSING::
compute_autocorrelation_OMP(MATRIX& time_series_data, const MKL_LONG index,
                            const MKL_LONG N_THREADS,
                            MATRIX& result, const MKL_LONG index_result)
{
  MKL_LONG Nt_corr = time_series_data.rows/2;
  // MATRIX result(Nt_corr, 1, 0.);
#pragma omp parallel for default(none) if (N_THREADS > 1)    \
  shared(time_series_data, result, index, Nt_corr)          \
  num_threads(N_THREADS)
  for(MKL_LONG t=0; t<Nt_corr; t++)
    {
      for(MKL_LONG s=0; s<Nt_corr -1; s++)
        {
          result(t) += time_series_data(s, index)*time_series_data(s + t, index);
        }
      result(t) /= (double)(Nt_corr - 1);
    }
  return 0.;
}
