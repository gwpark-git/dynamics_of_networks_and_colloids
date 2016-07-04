
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
  
  
  /* MKL_LONG (*compute_RDIST_particle)(RDIST&, const MKL_LONG, TRAJECTORY&, MKL_LONG); */
  bool SIMPLE_SHEAR;
  MKL_LONG shear_axis;
  MKL_LONG shear_grad_axis;
  double map_to_central_box_image;
  /* double* shear_variables; // it will store information related with simple shear flow */
  
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

  // member function case dependency
  // function point to measure minimum distance between pair of micelle

  // in order to use pointer-to-function rather than pointer-to-member_function, the measuring_minimum_distance have the type argument with itself: RDIST&. Note that the allocation is related with inline wrapping function defined in GEOMETRY namespace
  double (*measure_minimum_distance)(RDIST&, TRAJECTORY&, const MKL_LONG, const MKL_LONG, const MKL_LONG, MKL_LONG*);


  // the following are related with pointer-to-member_function.
  // in consequence of poingter-to-member_function, the fellowing have low readability which is not satisfactory by myself. Hence, the fellowing methods are deleted and left as recorded history to explain. (in future, these will be deleted completely) 

  // by some reason, the pointer-to-member_function is slightly differ from pointer-to-function.
  // it is not sure what exactly differ each other, but if this allocation step comapred with member function for POTENTIAL_SET, it uses typical allocation for pointer and uses it as name (as reference variable)
  // while the case on here is just works as pointer.
  // In consequence, allocation uses reference type rather than typical pointer type.
  // in future, the explanation will be revisited

  
  /* double (RDIST::*measure_minimum_distance)(TRAJECTORY&, const MKL_LONG, const MKL_LONG, const MKL_LONG, MKL_LONG*); */
  /* // measure_minimum_distance will be allocatged between following member function of RDIST. Hence, the type for function pointer should have RDIST:: notation.   */

  /* // the following are related with specific cases */
  /* // note that the definition inside class declaration will be inlined, which will not occur overhead for mapping functions. */
  /* double measure_minimum_distance_default(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check); */
  /* // { */
  /* //   // note that the beyond_box_check is not necessary for default measuring_minimum_distance function */
  /* //   // however, it is set in order to achive compatibility with cell-list implementation (for function pointer) */
  /* //   return GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, index_target, Rvec[index_particle][index_target]); */
  /* // } */

  /* double measure_minimum_distance_cell_list(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check); */
  /* // { */
  /* //   return GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t, index_particle, index_target, Rvec[index_particle][index_target], beyond_box_check); */
  /* // } */

  /* double measure_minimum_distance_simple_shear(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check); */
  /* // { */
  /* //   // here again, beyond_box_check is not necessary */
  /* //   // just added because of argument design */
  /* //   return GEOMETRY::get_minimum_distance_simple_shear(TRAJ, index_t, index_particle, index_target, Rvec[index_particle][index_target], shear_aixs, shear_grad_axis, map_to_central_box_image); */
  /* // } */
  
};




namespace UTIL_ARR
{
  double get_minimum_image_k_from_x(double x, double k, double dimension);
  MKL_LONG get_index_minimum_abs(double *k, MKL_LONG N);
}


namespace GEOMETRY
{


  
  MKL_LONG compute_RDIST_particle(RDIST& R_boost, const MKL_LONG index_particle, TRAJECTORY& TRAJ, MKL_LONG index_t);
  // the original set
  double get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec);
  /* double get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check); */
  double get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost); 
  double return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j);
  double minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t);
  double minimum_image_convention_loop(TRAJECTORY& TRAJ, MKL_LONG target_t);
  double apply_shear_boundary_condition(TRAJECTORY& TRAJ, MKL_LONG target_t, const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis, const double shift_factor);

  
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  double get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i);

  double get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i);


  // inlined functions
  double measure_minimum_distance_default(RDIST& R_boost, TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check);
  double measure_minimum_distance_cell_list(RDIST& R_boost, TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check);
  double measure_minimum_distance_simple_shear(RDIST& R_boost, TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check);  

  
  MKL_LONG get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check);
  double get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check);
  double get_minimum_distance_pos_vector_simple_shear(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_i, const MKL_LONG index_j, MATRIX& given_vec, const MKL_LONG shear_axis, const double map_to_central_box_image);
  double get_minimum_distance_rel_vector_simple_shear(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_i, const MKL_LONG index_j, MATRIX& given_vec, const MKL_LONG shear_axis, const double map_to_central_box_image);
  double get_minimum_distance_simple_shear(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis, const double map_to_central_box_image);
  MKL_LONG minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t, MKL_LONG index_particle);

} // namespace GEOMETRY



// the following are related with inlined function
// it is extracted from above for readability
inline double UTIL_ARR::get_minimum_image_k_from_x(double x, double k, double dimension)
{
  double kd[3] = {k-dimension - x, k - x, k+dimension - x};
  double re= kd[get_index_minimum_abs(kd, 3)] + x;
  return re;
}

inline MKL_LONG UTIL_ARR::get_index_minimum_abs(double *k, MKL_LONG N)
{
  MKL_LONG re = 0;
  for(MKL_LONG i=1; i<N; i++)
    {
      if (fabs(k[i]) < fabs(k[re]))
        re = i;
    }
  return re;
}

inline MKL_LONG GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check)
{
  // direction convention:
  // +: direction to the given bead
  // -: direction to the target bead

  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      // the scheme is differ from the default minimum distance definition
      // since there are beyond_box_check will identify a pair of particles have different beloning for PBC box
      // the cell list already have information related with a pair of micelles wheather the micelles are belong to different PBC box or not
      // in equilibrium simulation, beyond_box_check is constant, so there are no meaning to minimize distance for given all possible images of the pair of micelles.
      double PBC_coord_target = TRAJ(index_t, target_index, k) + (double)beyond_box_check[k]*TRAJ.box_dimension[k];
      given_vec(k) = PBC_coord_target - TRAJ(index_t, given_index, k);
    }
  return 0;
}
  
inline double GEOMETRY::get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check)
{
  GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJ, index_t, index_i, index_j, given_vec, beyond_box_check);

  return given_vec.norm();
  /* return cblas_dnrm2(given_vec.size, given_vec.data, 1); // 1 denote increment for x */
  /* the cblas_dnrm2 is slower than inlined norm */
}

  
inline double GEOMETRY::get_minimum_distance_pos_vector_simple_shear(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_i, const MKL_LONG index_j, MATRIX& given_vec, const MKL_LONG shear_axis, const double map_to_central_box_image)
{
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      double perturbed_image_coordinate = TRAJ(index_t, index_i, k);
      if(k == shear_axis)
        perturbed_image_coordinate += map_to_central_box_image;
      given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(perturbed_image_coordinate, TRAJ(index_t, index_j, k), TRAJ.box_dimension[k]);
    }
  return 0;
}

inline double GEOMETRY::get_minimum_distance_rel_vector_simple_shear(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_i, const MKL_LONG index_j, MATRIX& given_vec, const MKL_LONG shear_axis, const double map_to_central_box_image)
{
  GEOMETRY::get_minimum_distance_pos_vector_simple_shear(TRAJ, index_t, index_i, index_j, given_vec, shear_axis, map_to_central_box_image);
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      // direction convention:
      // +: direction to the given bead
      // -: direction to the target bead
      given_vec(k) -= TRAJ(index_t, index_i, k);
    }
  return 0;
}
  
inline double GEOMETRY::get_minimum_distance_simple_shear(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis, const double map_to_central_box_image)
{
  GEOMETRY::get_minimum_distance_rel_vector_simple_shear(TRAJ, index_t, index_i, index_j, given_vec, shear_axis, map_to_central_box_image);
  return given_vec.norm();
}

  
inline MKL_LONG GEOMETRY::minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t, MKL_LONG index_particle)
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

inline double GEOMETRY::get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec)
{
  // this comentation is related with the assumes that
  // the given dimensionality is properly working
  // therefore, it is assumed that there is no mistake
  // and no internal check for the dimensionality.
  GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, index_j, given_vec);
  return given_vec.norm();
}

inline double GEOMETRY::return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j)
{
  MATRIX rel_vec(TRAJ.N_dimension, 1, 0.);
  GEOMETRY::get_minimum_distance(TRAJ, index_t, index_i, index_j, rel_vec);
  return rel_vec.norm();
}

inline MKL_LONG GEOMETRY::get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
{
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ(index_t, given_index, k), TRAJ(index_t, target_index, k), TRAJ.box_dimension[k]);
    }
  return 0;
}

inline MKL_LONG GEOMETRY::get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
{
  GEOMETRY::get_minimum_distance_pos_vector(TRAJ, index_t, given_index, target_index, given_vec);
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      // direction convention:
      // +: direction to the given bead
      // -: direction to the target bead
      given_vec(k) -= TRAJ(index_t, given_index, k);
    }
  return 0;
}



inline double GEOMETRY::get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i)
{
  double distance = 0.;
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      distance += pow(TRAJ(index_t, index_i, k) - TRAJ(index_t, index_i, k), 2.0);
    }
  return sqrt(distance);

}


// inline double RDIST::measure_minimum_distance_default(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check)
// {
//   // note that the beyond_box_check is not necessary for default measuring_minimum_distance function
//   // however, it is set in order to achive compatibility with cell-list implementation (for function pointer)
//   return GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, index_target, Rvec[index_particle][index_target]);
// }

// inline double RDIST::measure_minimum_distance_cell_list(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check)
// {
//   return GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t, index_particle, index_target, Rvec[index_particle][index_target], beyond_box_check);
// }

// inline double RDIST::measure_minimum_distance_simple_shear(TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check)
// {
//   // here again, beyond_box_check is not necessary
//   // just added because of argument design
//   return GEOMETRY::get_minimum_distance_simple_shear(TRAJ, index_t, index_particle, index_target, Rvec[index_particle][index_target], shear_axis, shear_grad_axis, map_to_central_box_image);
// }

inline double GEOMETRY::measure_minimum_distance_default(RDIST& R_boost, TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
}

inline double GEOMETRY::measure_minimum_distance_cell_list(RDIST& R_boost, TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target], beyond_box_check);
}

inline double GEOMETRY::measure_minimum_distance_simple_shear(RDIST& R_boost, TRAJECTORY& TRAJ, const MKL_LONG index_t, const MKL_LONG index_particle, const MKL_LONG index_target, MKL_LONG* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance_simple_shear(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.shear_axis, R_boost.shear_grad_axis, R_boost.map_to_central_box_image);
}

#endif
