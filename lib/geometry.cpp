
#include "geometry.h"

MKL_LONG GEOMETRY::compute_RDIST_particle(RDIST& R_boost, const MKL_LONG index_particle, TRAJECTORY& TRAJ, MKL_LONG index_t)
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


RDIST::RDIST(COND& given_condition) : CLIST(given_condition)
{
  Rvec = new MATRIX* [Np];
  Rsca = new MATRIX [Np];
    
  for(MKL_LONG i=0; i<Np; i++)
    {
      // since the space complexity is not the matter for our simulation (at this moment),
      // the Rvec have Np*Np array that is much larger when we used cell-list approaches
      // Rvec[i] = (MATRIX*)mkl_malloc(Np*sizeof(MATRIX), BIT);
      Rvec[i] = new MATRIX [Np];
      for(MKL_LONG j=0; j<Np; j++)
        {
          Rvec[i][j].initial(N_dimension, 1, 0.);
        }
      Rsca[i].initial(Np, 1, 0.);
    }
  if(given_condition("SIMPLE_SHEAR") == "TRUE")
    {
      if(CELL_LIST_BOOST)
        {
          printf("ERR: CELL LIST with SIMPLE SHEAR FLOW is not implemented. \n");
        }
      SIMPLE_SHEAR = TRUE;
      shear_axis = atoi(given_condition("shear_axis").c_str());
      shear_grad_axis = atoi(given_condition("shear_grad_axis").c_str());
      map_to_central_box_image = 0.; // started with zero (equilibrium PBC box)
      // note again, the cell list implementation is not work at this moment

      measure_minimum_distance = GEOMETRY::measure_minimum_distance_simple_shear;
    }
  else
    {
      if(CELL_LIST_BOOST)
        measure_minimum_distance = GEOMETRY::measure_minimum_distance_cell_list;
      else
        measure_minimum_distance = GEOMETRY::measure_minimum_distance_default;
    }
}

double GEOMETRY::get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost)
{
  // this function generating minimum relative vector into R_minimum_vec_boost,
  // and distance of it into R_minimum_boost.
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      R_minimum_boost_particle(i) = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, i, R_minimum_vec_boost[index_particle][i]);
    }
  return 0;
}



double GEOMETRY::minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t)
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

double GEOMETRY::minimum_image_convention_loop(TRAJECTORY& TRAJ, MKL_LONG target_t)
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

double GEOMETRY::apply_shear_boundary_condition(TRAJECTORY& TRAJ, MKL_LONG target_t, const MKL_LONG shear_axis, const MKL_LONG shear_grad_axis, const double shift_factor)
{
  double time_st = dsecnd();
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      double coord = TRAJ(target_t, i, shear_grad_axis);
      if(coord < 0 || coord >= TRAJ.box_dimension[shear_grad_axis])
        {
          double sign = coord/fabs(coord);
          // TRAJ(target_t, i, shear_axis) -= fmod(sign*shift_factor, TRAJ.box_dimension[shear_axis]);
          TRAJ(target_t, i, shear_axis) -= sign*shift_factor;
          
        }
    }
  return dsecnd() - time_st;
}


