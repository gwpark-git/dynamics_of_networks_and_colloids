
#include "time_evolution.h"

// double INTEGRATOR::time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next, const MKL_LONG index_bead, MATRIX** contributions, const MKL_LONG N_forces) 
// {
//   // contributions[0] represent the Euler
//   // contributions[1] represent mechanical perturbation (if exist)
//   // contributions[2] represent repulsive potential (if exist)
//   // contributions[3] represent interaction through association (if exist)
//   for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//     {
//       TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k); // inheritance the current positions
//       for(MKL_LONG ind_f=0; ind_f<N_forces; ind_f++)
//         {
//           // TRAJ(index_t_next, i, k) += TRAJ.dt*force_array[ind_f][i](k);
//         }
//     }
//   return 0;
// }

double INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost)
{
  double time_st = dsecnd();
  given_vec.set_value(0.);
  for (MKL_LONG j=1; j<CONNECT.TOKEN[given_index]; j++)
    {
      MKL_LONG target_index = (MKL_LONG)CONNECT.HASH[given_index](j);
      MATRIX& rel_vector = R_minimum_vec_boost[given_index][target_index];
      double distance = R_minimum_distance_boost[given_index](target_index);
      double force = CONNECT.weight[given_index](j)*POTs.f_connector(distance, POTs.force_variables);
      // a*x + y -> y
      cblas_daxpy(given_vec.size,
                  force/distance, // a
                  rel_vector.data, // x
                  1,
                  given_vec.data, // y
                  1);
      
    }
  return dsecnd() - time_st;
}

double INTEGRATOR::EULER::cal_repulsion_force_R_boost(POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_particle, RDIST& R_boost)
{
  double time_st = dsecnd();
  given_vec.set_value(0.);
  MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
  for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
    {
      MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
      for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
        {
          MKL_LONG index_target = R_boost(cell_index_neighbor, p);
          double distance = R_boost.Rsca[index_particle](index_target);
          if (index_target != index_particle)
            {
              double repulsion = POTs.f_repulsion(distance, POTs.force_variables);
              // for(MKL_LONG d=0; d<R_boost.N_dimension; d++)
              //   {
              //     given_vec(d) += repulsion*R_boost.Rvec[index_particle][index_target](d);
              //   }
	      // the following make overhead.
              cblas_daxpy(given_vec.size,
                          repulsion/distance,
                          R_boost.Rvec[index_particle][index_target].data,
                          1,
                          given_vec.data,
                          1);
            }
        }
    }
  return dsecnd() - time_st;
}


double INTEGRATOR::EULER::cal_random_force_boost(POTENTIAL_SET& POTs, MATRIX& given_vec, gsl_rng* r_boost)
{
  double time_st = dsecnd();
  RANDOM::single_random_vector_generator_variance_boost(given_vec, 1.0, r_boost);
  matrix_mul(given_vec, sqrt(2.0));
  POTs.scale_random(given_vec, POTs.force_variables);
  return dsecnd() - time_st;
}

double INTEGRATOR::EULER::cal_random_force_boost_simplified(POTENTIAL_SET& POTs, MATRIX& given_vec, gsl_rng* r_boost)
{
  // this is duplicate one with INTEGRATOR::EULER::cal_random_force_boost
  // however, it simplified and reduce the overhead
  double time_st = dsecnd();
  double variance = 1.0;
  double pre_factor = sqrt(3.0)*sqrt(2.0);
  for(MKL_LONG k=0; k<given_vec.size; k++)
    {
      given_vec(k) = pre_factor*variance*(2.0*(gsl_rng_uniform(r_boost) - 0.5));
    }
  return dsecnd() - time_st;
}

double ANALYSIS::cal_total_energy_R_boost(POTENTIAL_SET& POTs, RDIST& R_boost)
{
  return cal_potential_energy_R_boost(POTs, R_boost);// + cal_kinetic_energy(TRAJ, POTs, index_t);
}

double ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy_R_boost(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& tmp_vec, RDIST& R_boost)
{
  double energy = 0.;

  for(MKL_LONG i=0; i<CONNECT.Np; i++)
    {
      for(MKL_LONG j=0; j<CONNECT.TOKEN[i]; j++)
        {
          // energy += CONNECT.weight[i](j)*POTs.e_connector(GEOMETRY::get_minimum_distance(TRAJ, index_t, i, (MKL_LONG)CONNECT.HASH[i](j), tmp_vec), POTs.force_variables);
          energy += CONNECT.weight[i](j)*POTs.e_connector(R_boost.Rsca[i]((MKL_LONG)CONNECT.HASH[i](j)), POTs.force_variables);
        }
    }
  return energy + ANALYSIS::cal_potential_energy_R_boost(POTs, R_boost);
}

double ANALYSIS::cal_potential_energy_R_boost(POTENTIAL_SET& POTs, RDIST& R_boost)
{
  double energy = 0.;
  MATRIX tmp_vec(R_boost.N_dimension, 1, 0.);
  for(MKL_LONG index_particle=0; index_particle<R_boost.Np; index_particle++)
    {
      MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
      // printf("cell_index[%d] = %d\n", index_particle, cell_index_particle);
      for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
        {
          MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
	  // printf("NEIGHBOR_CELLS[%d][%d] = %d\n", cell_index_particle, k, cell_index_neighbor);
          for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
            {
              MKL_LONG index_target = R_boost(cell_index_neighbor, p);
              double distance = R_boost.Rsca[index_particle](index_target);
              double U_ij = POTs.e_repulsion(distance, POTs.force_variables);
              energy += U_ij;
            }
        }
    }
  return energy;
}
double ANALYSIS::CAL_ENERGY_BROWNIAN(POTENTIAL_SET& POTs, MATRIX& mat_energy, double time)
{
  double time_st = dsecnd();
  mat_energy(0) = time;
  mat_energy(2) = 0.; // this is temporal setting
  // mat_energy(2) = cal_potential_energy_R_boost(POTs, R_boost);
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return dsecnd() - time_st;
}


double ANALYSIS::CAL_ENERGY_R_boost(POTENTIAL_SET& POTs, MATRIX& mat_energy, double time, RDIST& R_boost)
{
  double time_st = dsecnd();
  mat_energy(0) = time;
  mat_energy(2) = cal_potential_energy_R_boost(POTs, R_boost);
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return dsecnd() - time_st;
}

double ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY_R_boost(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, double time, MATRIX& tmp_vec, RDIST& R_boost)
{
  // mat_energy(0) = (TRAJ.c_t - 1)*TRAJ.dt;
  double time_st = dsecnd();
  mat_energy(0) = time;
  mat_energy(2) = ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy_R_boost(POTs, CONNECT, tmp_vec, R_boost);
  
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return dsecnd() - time_st;
}

