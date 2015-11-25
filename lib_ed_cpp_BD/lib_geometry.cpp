
#include "lib_geometry.h"

double GEOMETRY::get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec)
{
  // this comentation is related with the assumes that
  // the given dimensionality is properly working
  // therefore, it is assumed that there is no mistake
  // and no internal check for the dimensionality.
  GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, index_j, given_vec);
  return given_vec.norm();
}

double GEOMETRY::return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j)
{
  MATRIX rel_vec(TRAJ.dimension, 1, 0.);
  GEOMETRY::get_minimum_distance(TRAJ, index_t, index_i, index_j, rel_vec);
  return rel_vec.norm();
}

MKL_LONG GEOMETRY::minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t)
{
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      for (MKL_LONG k=0; k<TRAJ.dimension; k++)
        {
          double diff = TRAJ(target_t, i, k) - 0.5*TRAJ.box_dimension[k];
          double sign = diff/fabs(diff);
          if (fabs(diff) > 0.5*TRAJ.box_dimension[k])
            {
              TRAJ(target_t, i, k) -= sign*TRAJ.box_dimension[k];
            }
        }
    }
  return 0;
}

MKL_LONG GEOMETRY::get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
{
  for(MKL_LONG k=0; k<TRAJ.dimension; k++)
    {
      given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ(index_t, given_index, k), TRAJ(index_t, target_index, k), TRAJ.box_dimension[k]);
    }
  return 0;
}

MKL_LONG GEOMETRY::get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
{
  GEOMETRY::get_minimum_distance_pos_vector(TRAJ, index_t, given_index, target_index, given_vec);
  for(MKL_LONG k=0; k<TRAJ.dimension; k++)
    {
      // direction convention:
      // +: direction to the given bead
      // -: direction to the target bead
      given_vec(k) -= TRAJ(index_t, given_index, k);
    }
  return 0;
}

double UTIL_ARR::get_minimum_image_k_from_x(double x, double k, double dimension)
 {
   double kd[3] = {k-dimension - x, k - x, k+dimension - x};
   double re= kd[get_index_minimum_abs(kd, 3)] + x;
   return re;
 }

MKL_LONG UTIL_ARR::get_index_minimum_abs(double *k, MKL_LONG N)
 {
   MKL_LONG re = 0;
   for(MKL_LONG i=1; i<N; i++)
     {
       if (fabs(k[i]) < fabs(k[re]))
         re = i;
     }
   return re;
 }

double GEOMETRY::get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i)
{
  double distance = 0.;
  for(MKL_LONG k=0; k<TRAJ.dimension; k++)
    {
      distance += pow(TRAJ(index_t, index_i, k) - TRAJ(index_t, index_i, k), 2.0);
    }
  return sqrt(distance);

}
