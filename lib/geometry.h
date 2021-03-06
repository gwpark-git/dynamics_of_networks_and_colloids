
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

  MATRIX** Rvec; // Rvec[i][j] will be the relative vector between i- and j-th particles
  MATRIX* Rsca; // Rsca[i](j) will be the relative distance between i- and j-th particles. The form is differ from the original.
  
  
  RDIST()
  {
    std::cout << "There is no empty constructor for RDIST class\n" << std::endl;
  }
  RDIST(COND& given_condition);
  virtual
    ~RDIST()
  {
    if(INITIALIZATION)
      {
        delete[] Rsca;
        for(MKL_LONG i=0; i<Np; i++)
          {
            delete[] Rvec[i];
          }
        delete[] Rvec;
      }
  }

  // member function case dependency
  // function point to measure minimum distance between pair of micelle

  // in order to use pointer-to-function rather than pointer-to-member_function, the measuring_minimum_distance have the type argument with itself: RDIST&. Note that the allocation is related with inline wrapping function defined in GEOMETRY namespace
  double
  (*measure_minimum_distance)(RDIST&, TRAJECTORY&, const MKL_LONG, const MKL_LONG, const MKL_LONG, MKL_LONG*);


  // the following are related with pointer-to-member_function.
  // in consequence of poingter-to-member_function, the fellowing have low readability which is not satisfactory by myself. Hence, the fellowing methods are deleted and left as recorded history to explain. (in future, these will be deleted completely) 

  // by some reason, the pointer-to-member_function is slightly differ from pointer-to-function.
  // it is not sure what exactly differ each other, but if this allocation step comapred with member function for POTENTIAL_SET, it uses typical allocation for pointer and uses it as name (as reference variable)
  // while the case on here is just works as pointer.
  // In consequence, allocation uses reference type rather than typical pointer type.
  // in future, the explanation will be revisited

  
};

namespace UTIL_ARR
{
  double
    get_minimum_image_k_from_x(double x, double k, double dimension);

  MKL_LONG
    get_index_minimum_abs(double *k, MKL_LONG N);
}


namespace GEOMETRY
{
  
  MKL_LONG
    compute_RDIST_particle(RDIST& R_boost, const MKL_LONG index_particle,
			   TRAJECTORY& TRAJ, MKL_LONG index_t);
  // the original set
  double
    get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t,
			 MKL_LONG index_i, MKL_LONG index_j,
			 MATRIX& given_vec);
  double
    get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t,
				      MKL_LONG index_particle,
				      MATRIX& R_minimum_boost_particle,
				      MATRIX** R_minimum_vec_boost); 
  double
    return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t,
			    MKL_LONG index_i, MKL_LONG index_j);
  double
    minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t);
  double
    minimum_image_convention_loop(TRAJECTORY& TRAJ, MKL_LONG target_t);
  double
  apply_shear_boundary_condition(TRAJECTORY& TRAJ, MKL_LONG target_t,
				   const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis, const double shift_factor);

  double
  set_box_shift_factor(double &shear_PBC_shift,
                       double const maximum_displacement,
                       double const box_dimension_shear_axis);
  
  double
  set_box_shift_factor(double &shear_PBC_shift,
                       RDIST& R_boost,
                       double const maximum_displacement,
                       double const box_dimension_shear_axis);
  
  double
  apply_step_shear(TRAJECTORY& TRAJ, MKL_LONG target_t,
                   const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis,
                   const double gamma_0, const double box_dimension);
  
  MKL_LONG
    get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t,
				    MKL_LONG given_index, MKL_LONG target_index,
				    MATRIX& given_vec);
  MKL_LONG
    get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t,
				    MKL_LONG given_index, MKL_LONG target_index,
				    MATRIX& given_vec);
  MKL_LONG
    get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t,
				    MKL_LONG given_index, MKL_LONG target_index,
				    MATRIX& given_vec);
  double
    get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t,
			MKL_LONG index_i, MKL_LONG i);

  double
    get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t,
			MKL_LONG index_i, MKL_LONG i);


  // inlined functions
  double
    measure_minimum_distance_default(RDIST& R_boost,
				     TRAJECTORY& TRAJ, const MKL_LONG index_t,
				     const MKL_LONG index_particle, const MKL_LONG index_target,
				     MKL_LONG* beyond_box_check);
  double
    measure_minimum_distance_cell_list(RDIST& R_boost,
				       TRAJECTORY& TRAJ, const MKL_LONG index_t,
				       const MKL_LONG index_particle, const MKL_LONG index_target,
				       MKL_LONG* beyond_box_check);
  double
    measure_minimum_distance_simple_shear_fixed_axis(RDIST& R_boost,
						     TRAJECTORY& TRAJ, const MKL_LONG index_t,
						     const MKL_LONG index_particle, const MKL_LONG index_target,
						     MKL_LONG* beyond_box_check);  

  
  MKL_LONG
    get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t,
					      MKL_LONG given_index, MKL_LONG target_index,
					      MATRIX& given_vec,
					      MKL_LONG* beyond_box_check);
  double
    get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t,
				   MKL_LONG index_i, MKL_LONG index_j,
				   MATRIX& given_vec,
				   MKL_LONG* beyond_box_check);
  double
    get_minimum_distance_rel_vector_simple_shear_fixed_axis(TRAJECTORY& TRAJ, const MKL_LONG index_t,
							    const MKL_LONG index_i, const MKL_LONG index_j,
							    MATRIX& given_vec,
							    const double map_to_central_box_image);
  double
    get_minimum_distance_simple_shear_fixed_axis(TRAJECTORY& TRAJ, MKL_LONG index_t,
						 MKL_LONG index_i, MKL_LONG index_j,
						 MATRIX& given_vec,
						 const double map_to_central_box_image);
  MKL_LONG
    minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t,
				      MKL_LONG index_particle);
  double
    get_minimum_distance_rel_vector_equilibrium_test(TRAJECTORY& TRAJ, const MKL_LONG index_t,
						     const MKL_LONG index_i, const MKL_LONG index_j,
						     MATRIX& given_vec);
  double
    get_minimum_distance_equilibrium_test(TRAJECTORY& TRAJ, MKL_LONG index_t,
					  MKL_LONG index_i, MKL_LONG index_j,
					  MATRIX& given_vec);
  MKL_LONG
    minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t,
				      MKL_LONG index_particle);
  
} // namespace GEOMETRY

// the following are related with inlined function
// it is extracted from above for readability
inline double
UTIL_ARR::
get_minimum_image_k_from_x(double x, double k, double dimension)
{
  double kd[3] = {k-dimension - x, k - x, k+dimension - x};
  double re= kd[get_index_minimum_abs(kd, 3)] + x;
  return re;
}

inline MKL_LONG
UTIL_ARR::
get_index_minimum_abs(double *k, MKL_LONG N)
{
  MKL_LONG re = 0;
  for(MKL_LONG i=1; i<N; i++)
    {
      if (fabs(k[i]) < fabs(k[re]))
        re = i;
    }
  return re;
}

inline MKL_LONG
GEOMETRY::
get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t,
					  MKL_LONG given_index, MKL_LONG target_index,
					  MATRIX& given_vec,
					  MKL_LONG* beyond_box_check)
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
  
inline double
GEOMETRY::
get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t,
			       MKL_LONG index_i, MKL_LONG index_j,
			       MATRIX& given_vec,
			       MKL_LONG* beyond_box_check)
{
  GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJ, index_t, index_i, index_j, given_vec, beyond_box_check);

  return given_vec.norm();
  /* return cblas_dnrm2(given_vec.size, given_vec.data, 1); // 1 denote increment for x */
  /* the cblas_dnrm2 is slower than inlined norm */
}


inline double
GEOMETRY::
get_minimum_distance_rel_vector_equilibrium_test(TRAJECTORY& TRAJ, const MKL_LONG index_t,
						 const MKL_LONG index_i, const MKL_LONG index_j,
						 MATRIX& given_vec)
{
  double min_norm = TRAJ.box_dimension[0];
  MKL_LONG min_sf_0 = 0, min_sf_1 = 0, min_sf_2 = 0;
  for(MKL_LONG sf_0=-1; sf_0<=+1; sf_0++)
    {
      for(MKL_LONG sf_1=-1; sf_1<=+1; sf_1++)
        {
          for(MKL_LONG sf_2=-1; sf_2<=+1; sf_2++)
            {
              double rpos_ij_0 = TRAJ(index_t, index_j, 0) + sf_0*TRAJ.box_dimension[0] - TRAJ(index_t, index_i, 0);
              double rpos_ij_1 = TRAJ(index_t, index_j, 1) + sf_1*TRAJ.box_dimension[1] - TRAJ(index_t, index_i, 1);
              double rpos_ij_2 = TRAJ(index_t, index_j, 2) + sf_2*TRAJ.box_dimension[2] - TRAJ(index_t, index_i, 2);
              double norm_pos = sqrt(rpos_ij_0*rpos_ij_0 + rpos_ij_1*rpos_ij_1 + rpos_ij_2*rpos_ij_2);
              if(norm_pos < min_norm)
                {
                  min_sf_0 = sf_0;
                  min_sf_1 = sf_1;
                  min_sf_2 = sf_2;
                  min_norm = norm_pos;
                }
            } // it is already assumed for 3-dimensional case
        }
    }
  given_vec(0) = TRAJ(index_t, index_j, 0) + min_sf_0*TRAJ.box_dimension[0] - TRAJ(index_t, index_i, 0);
  given_vec(1) = TRAJ(index_t, index_j, 1) + min_sf_1*TRAJ.box_dimension[1] - TRAJ(index_t, index_i, 1);
  given_vec(2) = TRAJ(index_t, index_j, 2) + min_sf_2*TRAJ.box_dimension[2] - TRAJ(index_t, index_i, 2);
  return min_norm;
}



inline double
GEOMETRY::
get_minimum_distance_equilibrium_test(TRAJECTORY& TRAJ, const MKL_LONG index_t,
				      const MKL_LONG index_i, const MKL_LONG index_j,
				      MATRIX& given_vec)
{
  return GEOMETRY::get_minimum_distance_rel_vector_equilibrium_test(TRAJ, index_t, index_i, index_j, given_vec);
}

inline double
GEOMETRY::
get_minimum_distance_rel_vector_simple_shear_fixed_axis(TRAJECTORY& TRAJ, const MKL_LONG index_t,
							const MKL_LONG index_i, const MKL_LONG index_j,
							MATRIX& given_vec,
							const double map_to_central_box_image)
{
  /*
    Even if it slows, it have worth to try compare with the existing scheme.
    Note that the equivalent with get_minimum_distance is confirmed.
   */
  double min_norm = TRAJ.box_dimension[0];
  MKL_LONG min_sf_0 = 0, min_sf_1 = 0, min_sf_2 = 0;
  for(MKL_LONG sf_2=-1; sf_2<=+1; sf_2++)
    {
      double rpos_ij_2 = TRAJ(index_t, index_j, 2) + sf_2*TRAJ.box_dimension[2]
        - TRAJ(index_t, index_i, 2);
      for(MKL_LONG sf_1=-1; sf_1<=+1; sf_1++)
        {
          double rpos_ij_1 = TRAJ(index_t, index_j, 1) + sf_1*TRAJ.box_dimension[1]
            - TRAJ(index_t, index_i, 1);
          
          for(MKL_LONG sf_0=-1; sf_0<=+1; sf_0++)
            {
              double rpos_ij_0 = TRAJ(index_t, index_j, 0) + sf_0*TRAJ.box_dimension[0] + sf_1*map_to_central_box_image
                - TRAJ(index_t, index_i, 0);
          
              double norm_pos = sqrt(rpos_ij_0*rpos_ij_0 + rpos_ij_1*rpos_ij_1 + rpos_ij_2*rpos_ij_2);
              if(norm_pos < min_norm)
                {
                  min_sf_0 = sf_0;
                  min_sf_1 = sf_1;
                  min_sf_2 = sf_2;
                  min_norm = norm_pos;
                }
            } // it is already assumed for 3-dimensional case
        }
    }
  given_vec(0) = TRAJ(index_t, index_j, 0) + min_sf_0*TRAJ.box_dimension[0] + min_sf_1*map_to_central_box_image - TRAJ(index_t, index_i, 0);
  given_vec(1) = TRAJ(index_t, index_j, 1) + min_sf_1*TRAJ.box_dimension[1] - TRAJ(index_t, index_i, 1);
  given_vec(2) = TRAJ(index_t, index_j, 2) + min_sf_2*TRAJ.box_dimension[2] - TRAJ(index_t, index_i, 2);
  return min_norm;
}

inline double
GEOMETRY::
get_minimum_distance_simple_shear_fixed_axis(TRAJECTORY& TRAJ, MKL_LONG index_t,
					     MKL_LONG index_i, MKL_LONG index_j,
					     MATRIX& given_vec,
					     const double map_to_central_box_image)
{
  double distance = GEOMETRY::get_minimum_distance_rel_vector_simple_shear_fixed_axis(TRAJ, index_t, index_i, index_j, given_vec, map_to_central_box_image);
  return distance;
}

  
inline MKL_LONG
GEOMETRY::
minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t,
				  MKL_LONG index_particle)
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

inline double
GEOMETRY::
get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t,
		     MKL_LONG index_i, MKL_LONG index_j,
		     MATRIX& given_vec)
{
  // this comentation is related with the assumes that
  // the given dimensionality is properly working
  // therefore, it is assumed that there is no mistake
  // and no internal check for the dimensionality.
  GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, index_j, given_vec);
  return given_vec.norm();
}

inline double
GEOMETRY::
return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t,
			MKL_LONG index_i, MKL_LONG index_j)
{
  MATRIX rel_vec(TRAJ.N_dimension, 1, 0.);
  GEOMETRY::get_minimum_distance(TRAJ, index_t, index_i, index_j, rel_vec);
  return rel_vec.norm();
}

inline MKL_LONG
GEOMETRY::
get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t,
				MKL_LONG given_index, MKL_LONG target_index,
				MATRIX& given_vec)
{
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ(index_t, given_index, k), TRAJ(index_t, target_index, k), TRAJ.box_dimension[k]);
    }
  return 0;
}

inline MKL_LONG
GEOMETRY::
get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t,
				MKL_LONG given_index, MKL_LONG target_index,
				MATRIX& given_vec)
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

inline double
GEOMETRY::
get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t,
		    MKL_LONG index_i, MKL_LONG i)
{
  double distance = 0.;
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      distance += pow(TRAJ(index_t, index_i, k) - TRAJ(index_t, index_i, k), 2.0);
    }
  return sqrt(distance);

}

inline double
GEOMETRY::
measure_minimum_distance_default(RDIST& R_boost,
				 TRAJECTORY& TRAJ, const MKL_LONG index_t,
				 const MKL_LONG index_particle, const MKL_LONG index_target,
				 MKL_LONG* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
  // the following function is reference (and very slow).
  // to be on safe side, it still alive (just commented)
  // return GEOMETRY::get_minimum_distance_equilibrium_test(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
}

inline double
GEOMETRY::
measure_minimum_distance_cell_list(RDIST& R_boost,
				   TRAJECTORY& TRAJ, const MKL_LONG index_t,
				   const MKL_LONG index_particle, const MKL_LONG index_target,
				   MKL_LONG* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target], beyond_box_check);
}

inline double
GEOMETRY::
measure_minimum_distance_simple_shear_fixed_axis(RDIST& R_boost,
						 TRAJECTORY& TRAJ, const MKL_LONG index_t,
						 const MKL_LONG index_particle, const MKL_LONG index_target,
						 MKL_LONG* beyond_box_check)
{
  return GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(TRAJ, index_t, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.map_to_central_box_image);
}

#endif
