
#include "lib_geometry.h"

double GEOMETRY::get_minimum_distance(TRAJECTORY& TRAJ, long index_t, long index_i, long index_j, MATRIX& given_vec)
{
  // this comentation is related with the assumes that
  // the given dimensionality is properly working
  // therefore, it is assumed that there is no mistake
  // and no internal check for the dimensionality.
  // if(given_vec.size != TRAJ.dimension)
  //   {
  //     given_vec.initial(TRAJ.dimension, 1, 0.);
  //   }
  // printf("--> FLAG_GET_MINIMUM_DISTANCE_1\n");
  GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, index_j, given_vec);
  // printf("--> FLAG_GET_MINIMUM_DISTANCE_2\n");
  return given_vec.norm();
}

double GEOMETRY::return_minimum_distance(TRAJECTORY& TRAJ, long index_t, long index_i, long index_j)
{
  MATRIX rel_vec(TRAJ.dimension, 1, 0.);
  GEOMETRY::get_minimum_distance(TRAJ, index_t, index_i, index_j, rel_vec);
  return rel_vec.norm();
}



long GEOMETRY::minimum_image_convention(TRAJECTORY& TRAJ, long target_t)
{
  // target_t = target_t%TRAJ.Nt;
  for (long i=0; i<TRAJ.Np; i++)
    {
      for (long k=0; k<TRAJ.dimension; k++)
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

long GEOMETRY::get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, long index_t, long given_index, long target_index, MATRIX& given_vec)
{
  // if(given_vec.size != TRAJ.dimension)
  //   {
  //     given_vec.initial(TRAJ.dimension, 1, 0.);
  //   }
  for(long k=0; k<TRAJ.dimension; k++)
    {
      // printf("    --> FLAG_GET_MINIMUM_DISTANCE_POS_VECTOR_1\t%ld\t%ld\t%ld\t%6.3e\t%6.3e\n", index_t, given_index, target_index, 0.0, 0.0);
      
      given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ(index_t, given_index, k), TRAJ(index_t, target_index, k), TRAJ.box_dimension[k]);
      // printf("    --> FLAG_GET_MINIMUM_DISTANCE_POS_VECTOR_2\n");      
    }
  return 0;
}

long GEOMETRY::get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, long index_t, long given_index, long target_index, MATRIX& given_vec)
{
  // printf("  --> FLAG_GET_MINIMUM_DISTANCE_REL_VECTOR_1\n");
  GEOMETRY::get_minimum_distance_pos_vector(TRAJ, index_t, given_index, target_index, given_vec);
  // printf("  --> FLAG_GET_MINIMUM_DISTANCE_REL_VECTOR_2\n");
  for(long k=0; k<TRAJ.dimension; k++)
    {
      // direction convention:
      // +: direction to the given bead
      // -: direction to the target bead
      given_vec(k) -= TRAJ(index_t, given_index, k);
      // re(k) -= *(TRAJ.R_ref[index_t](given_index, k));
    }
  return 0;
}

double UTIL_ARR::get_minimum_image_k_from_x(double x, double k, double dimension)
 {
   // printf("      --> FLAG_GET_MINIMUM_IMAGE_K_FROM_K_0\n");
   double kd[3] = {k-dimension - x, k - x, k+dimension - x};
   // printf("      --> FLAG_GET_MINIMUM_IMAGE_K_FROM_K_1\n");
   double re= kd[get_index_minimum_abs(kd, 3)] + x;
   // printf("      --> FLAG_GET_MINIMUM_IMAGE_K_FROM_K_2\n");
   
   return re;
 }

long UTIL_ARR::get_index_minimum_abs(double *k, long N)
 {
   long re = 0;
   for(long i=1; i<N; i++)
     {
       if (fabs(k[i]) < fabs(k[re]))
         re = i;
     }
   return re;
 }

double GEOMETRY::get_simple_distance(TRAJECTORY& TRAJ, long index_t, long index_i, long i)
{
  double distance = 0.;
  for(long k=0; k<TRAJ.dimension; k++)
    {
      distance += pow(TRAJ(index_t, index_i, k) - TRAJ(index_t, index_i, k), 2.0);
    }
  return sqrt(distance);

}
