
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include "trajectory.h"
#include "read_file_condition.h"

class RDIST : public CLIST
{
 public:

  /* double** Rmin; */
  MATRIX** Rvec; // Rvec[i][j] will be the relative vector between i- and j-th particles
  MATRIX* Rsca; // Rsca[i](j) will be the relative distance between i- and j-th particles. The form is differ from the original.
  
  RDIST()
    {
      std::cout << "There is no empty constructor for RDIST class\n" << std::endl;
      return -1;
    }
 RDIST(COND& given_condition) : CLIST(given_condition)
    {
      Rvec = (MATRIX**)mkl_malloc(Np*sizeof(MATRIX*), BIT);
      Rsca = (MATRIX*)mkl_malloc(Np*sizeof(MATRIX), BIT);
      for(MKL_LONG i=0; i<Np; i++)
        {
          Rvec[i] = (MATRIX*)mkl_malloc(Np*sizeof(MATRIX), BIT);
          for(MKL_LONG j=0; j<Np; j++)
            {
              Rvec[i][j].initial(N_dimension, 1, 0.);
            }
          Rsca[i].initial(Np, 1, 0.);
        }
    }

};

namespace GEOMETRY
{
  double get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec);
  double get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost);  

  double return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j);
  MKL_LONG minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);

  double get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i);
}

namespace UTIL_ARR
{
  double get_minimum_image_k_from_x(double x, double k, double dimension);
  MKL_LONG get_index_minimum_abs(double *k, MKL_LONG N);
}


#endif
