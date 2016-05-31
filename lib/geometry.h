
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include "trajectory.h"
/* #include "read_file_condition.h" */
#include "file_IO.h"
#include "cell_list.h"
#include <mkl.h>

class RDIST : public CLIST
{
 public:

  /* double** Rmin; */
  MATRIX** Rvec; // Rvec[i][j] will be the relative vector between i- and j-th particles
  MATRIX* Rsca; // Rsca[i](j) will be the relative distance between i- and j-th particles. The form is differ from the original.
  RDIST()
    {
      std::cout << "There is no empty constructor for RDIST class\n" << std::endl;
    }
  RDIST(COND& given_condition);
  virtual ~RDIST()
    {
      if(INITIALIZATION)
        {
          /* mkl_free(Rsca); */
          delete[] Rsca;
          for(MKL_LONG i=0; i<Np; i++)
            {
              /* mkl_free(Rvec[i]); */
              delete[] Rvec[i];
            }
          /* mkl_free(Rvec); */
          delete[] Rvec;
        }
    }

  // member function
  MKL_LONG compute_RDIST_particle(const MKL_LONG index_particle, TRAJECTORY& TRAJ, MKL_LONG index_t);
};

namespace GEOMETRY
{
  // the original set
  double get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec);
  /* double get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check); */
  double get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost); 
  double return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j);
  double minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t);
  /* MKL_LONG minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t, MKL_LONG index_particle); */ // inlined
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  /* MKL_LONG get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check); */
  double get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i);

  // with newly defined trajectory_HDF5
  /* double get_minimum_distance(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec); */
  /* double get_minimum_distance_cell_list(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check); */
  /* double get_minimum_distance_for_particle(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost);  */
  /* double return_minimum_distance(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j); */
  /* MKL_LONG minimum_image_convention(TRAJECTORY_HDF5& TRAJ, MKL_LONG target_t); */
  /* MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec); */
  /* MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec); */
  /* MKL_LONG get_minimum_distance_rel_vector(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec); */
  /* MKL_LONG get_minimum_distance_rel_vector_cell_list(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check); */

  /* double get_simple_distance(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i); */

  

  double get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i);

  // inlined
  /* double get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check); */
  /* MKL_LONG get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check); */

  inline MKL_LONG get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check)
  {
    // for(MKL_LONG k=0; k<TRAJ.dimension; k++)
    //   {
    //     // direction convention:
    //     // +: direction to the given bead
    //     // -: direction to the target bead
    //     given_vec(k) -= TRAJ(index_t, given_index, k);
    //   }
    for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
      {
        double PBC_coord_target = TRAJ(index_t, target_index, k) + (double)beyond_box_check[k]*TRAJ.box_dimension[k];
        given_vec(k) = PBC_coord_target - TRAJ(index_t, given_index, k);
      }
    return 0;
  }
  
  inline double get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check)
  {
    get_minimum_distance_rel_vector_cell_list(TRAJ, index_t, index_i, index_j, given_vec, beyond_box_check);

    return given_vec.norm();
    /* return cblas_dnrm2(given_vec.size, given_vec.data, 1); // 1 denote increment for x */
    /* the cblas_dnrm2 is slower than inlined norm */
  }
  
  inline MKL_LONG minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t, MKL_LONG index_particle)
  {
    for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
      {
        double diff = TRAJ(target_t, index_particle, k) - 0.5*TRAJ.box_dimension[k];
        double sign = diff/fabs(diff);
        if (fabs(diff) > 0.5*TRAJ.box_dimension[k])
          {
            TRAJ(target_t, index_particle, k) -= sign*TRAJ.box_dimension[k];
          }
      }
    return 0;
  }



}

namespace UTIL_ARR
{
  double get_minimum_image_k_from_x(double x, double k, double dimension);
  MKL_LONG get_index_minimum_abs(double *k, MKL_LONG N);
}


#endif
