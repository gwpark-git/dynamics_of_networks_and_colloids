
#include "geometry.h"

MKL_LONG
GEOMETRY::
compute_RDIST_particle(RDIST& R_boost, const MKL_LONG index_particle,
		       TRAJECTORY& TRAJ, MKL_LONG index_t)
{
  MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
  for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
    {
      MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
      for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
        {
          MKL_LONG index_target = R_boost(cell_index_neighbor, p);
          double distance = R_boost.measure_minimum_distance(R_boost, TRAJ, index_t, index_particle, index_target, R_boost.BEYOND_BOX[cell_index_particle][k]);
          R_boost.Rsca[index_particle](index_target) = distance;
        }
    }
  return 0;
}


RDIST::
RDIST(COND& given_condition)
  : CLIST(given_condition)
{
  Rvec = new MATRIX* [Np];
  Rsca = new MATRIX [Np];
    
  for(MKL_LONG i=0; i<Np; i++)
    {
      // since the space complexity is not the matter for our simulation (at this moment),
      // the Rvec have Np*Np array that is much larger when we used cell-list approaches
      Rvec[i] = new MATRIX [Np];
      for(MKL_LONG j=0; j<Np; j++)
        {
          Rvec[i][j].initial(N_dimension, 1, 0.);
        }
      Rsca[i].initial(Np, 1, 0.);
    }
  // if(SIMPLE_SHEAR) // it already declare in CLIST class
  if(SIMPLE_SHEAR || STEP_SHEAR) 
    {
      // note that both of the cases are applicable based on the following functions
      // hence the name of functions are subject to change from SIMPLE_SHEAR to SHEAR
      if(CELL_LIST_BOOST)
        {
          printf("ERR: CELL LIST with SIMPLE SHEAR FLOW is not implemented. \n");
        }
      if (!(shear_axis == 0 && shear_grad_axis == 1) && !CELL_LIST_BOOST)
        {
          printf("ERR: SIMPLE SHEAR FLOW with different direction without CELL LIST is not implemented\n");
        }
      // note again, the cell list implementation is not work at this moment

      measure_minimum_distance = GEOMETRY::measure_minimum_distance_simple_shear_fixed_axis;
      // check
      printf("DISTANCE is using the function GEOMETRY::measure_minimum_distance_simple_shear_fixed_axis\n");
      
    }
  else
    {
      if(CELL_LIST_BOOST)
        measure_minimum_distance = GEOMETRY::measure_minimum_distance_cell_list;
      else
        measure_minimum_distance = GEOMETRY::measure_minimum_distance_default;
    }
}

double
GEOMETRY::
get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t,
				  MKL_LONG index_particle,
				  MATRIX& R_minimum_boost_particle,
				  MATRIX** R_minimum_vec_boost)
{
  // this function generating minimum relative vector into R_minimum_vec_boost,
  // and distance of it into R_minimum_boost.
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      R_minimum_boost_particle(i) = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, i, R_minimum_vec_boost[index_particle][i]);
    }
  return 0;
}

double
GEOMETRY::
minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t)
{
  double time_st = dsecnd();
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          double coord = TRAJ(target_t, i, k);
          if (coord < 0 || coord >= TRAJ.box_dimension[k]) // check left and right boundary
            {
              double sign = coord/fabs(coord);
              TRAJ(target_t, i, k) -= sign*TRAJ.box_dimension[k];
            }
          
        }
    }
  return dsecnd() - time_st;
}

double
GEOMETRY::
minimum_image_convention_loop(TRAJECTORY& TRAJ, MKL_LONG target_t)
{
  double time_st = dsecnd();
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          double coord = TRAJ(target_t, i, k);
          while(coord < 0 || coord >= TRAJ.box_dimension[k])
            {
              double sign = coord/fabs(coord);
              TRAJ(target_t, i, k) -= sign*TRAJ.box_dimension[k];
              coord = TRAJ(target_t, i, k);
            }
        }
    }
  return dsecnd() - time_st;

}

double
GEOMETRY::
apply_shear_boundary_condition(TRAJECTORY& TRAJ, MKL_LONG target_t,
                               const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis, const double shift_factor)
{
  double time_st = dsecnd();
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      double coord = TRAJ(target_t, i, shear_grad_axis);
      if(coord < 0 || coord >= TRAJ.box_dimension[shear_grad_axis])
        {
          double sign = coord/fabs(coord);
          TRAJ(target_t, i, shear_axis) -= sign*shift_factor;
          
        }
    }
  return dsecnd() - time_st;
}

double
GEOMETRY::
apply_step_shear(TRAJECTORY& TRAJ, MKL_LONG target_t,
                       const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis,
                       const double gamma_0, const double box_dimension)
{
  double time_st = dsecnd();
  printf("APPLY STEP DEFORMATION (STEP_SHEAR is on)\n");
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      // apply initial step deformation
      // note that shift in x-direction is given by: y*gamma_0 <= 1
      TRAJ(target_t, i, shear_axis) += gamma_0*TRAJ(target_t, i, shear_grad_axis);
    }
  // // the following is not necessary at this moment since the initial box deformation only apply for shear direction (shear gradient direction still in minimum image convention)
  // GEOMETRY::
  //   apply_shear_boundary_condition(TRAJ, target_t, shear_axis, shear_grad_axis, gamma_0*box_dimension);

  // following will apply minimum image convention (for shear direction)
  GEOMETRY::
    minimum_image_convention_loop(TRAJ, target_t);
  return dsecnd() - time_st;
}

double
GEOMETRY::
set_box_shift_factor(double &shear_PBC_shift,
                     double const maximum_displacement,
                     double const box_dimension_shear_axis)
{
  double time_st = dsecnd();
  shear_PBC_shift = fmod(maximum_displacement, box_dimension_shear_axis);
  return dsecnd() - time_st;
}

double
GEOMETRY::
set_box_shift_factor(double &shear_PBC_shift,
                     RDIST& R_boost,
                     double const maximum_displacement,
                     double const box_dimension_shear_axis)
{
  /* // Original Form
     VAR.shear_PBC_shift = fmod(VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis], TRAJ.box_dimension[VAR.shear_axis]);
     // VAR.shear_PBC_shift = VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis];
     R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
     MKL_LONG central_standard = (MKL_LONG)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
     R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
  */
  double time_st = dsecnd();
  shear_PBC_shift = fmod(maximum_displacement, box_dimension_shear_axis);
  R_boost.map_to_central_box_image = fmod(shear_PBC_shift, box_dimension_shear_axis);
  MKL_LONG central_standard = (MKL_LONG)(2*R_boost.map_to_central_box_image/box_dimension_shear_axis);
  R_boost.map_to_central_box_image -= box_dimension_shear_axis*(double)central_standard;
  return dsecnd() - time_st;
}
