
#include "lib_evolution.h"

long ANALYSIS::GET_dCDF_POTENTIAL_target(TRAJECTORY& TRAJ, long index_t, POTENTIAL_SET& POTs, long& index_particle, long& index_target, double& INDEX_dCDF_U_ij, double& dCDF_U_ij, MATRIX& vec_boost_ordered_pdf_ij)
{
  INDEX_dCDF_U_ij = (double)index_target;
  double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, index_target, vec_boost_ordered_pdf_ij);
  dCDF_U_ij = POTs.PDF_connector(distance, POTs.force_variables);
  return dCDF_U_ij;
}


long ANALYSIS::GET_dCDF_POTENTIAL(TRAJECTORY& TRAJ, long index_t, POTENTIAL_SET& POTs, long index_particle, MATRIX& INDEX_dCDF_U, MATRIX& dCDF_U, MATRIX& vec_boost_ordered_pdf)
{
  // note that the index_itself is not extracted in order to avoid if-phrase
  // even if it is not the case, the given distance will be zero
  // that gave us the potential as zero
  // then, Boltzmann distribution becomes unity, which is the maximum value
  for(long i=0; i<TRAJ.Np; i++)
    {
      double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, i, vec_boost_ordered_pdf);
      INDEX_dCDF_U(i) = i;
      dCDF_U(i) = POTs.PDF_connector(distance, POTs.force_variables);
    }
  return 0;
}

long ANALYSIS::GET_ORDERED_dCDF_POTENTIAL(TRAJECTORY& TRAJ, long index_t, POTENTIAL_SET& POTs, long index_particle, MATRIX& INDEX_dCDF_U, MATRIX& dCDF_U, MATRIX& vec_boost_ordered_pdf, double& time_dCDF, double& time_SORT)
{
  double time_st = dsecnd();
  GET_dCDF_POTENTIAL(TRAJ, index_t, POTs, index_particle, INDEX_dCDF_U, dCDF_U, vec_boost_ordered_pdf);
  double time_end_pdf = dsecnd();
  dCDF_U.sort2(INDEX_dCDF_U);
  double time_end_sort = dsecnd();
  time_dCDF += time_end_pdf - time_st;
  time_SORT += time_end_sort - time_end_pdf;
  return 0;
}

long INTEGRATOR::EULER::cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t, long given_index)
{
  given_vec.set_value(0.);
  return 0;
}

long INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, long index_t, long given_index, MATRIX& vec_boost_Nd)
{
  given_vec.set_value(0.);
  for (long j=1; j<CONNECT.TOKEN[given_index]; j++)
    {
      long target_index = (long)CONNECT.HASH[given_index](j);
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, given_index, target_index, vec_boost_Nd);
      double distance = cblas_dnrm2(vec_boost_Nd.size,
                                    vec_boost_Nd.data, // given vec
                                    1);
      double force = CONNECT.weight[given_index](j)*POTs.f_connector(distance, POTs.force_variables);
      // a*x + y -> y
      cblas_daxpy(given_vec.size,
                  force/distance, // a
                  vec_boost_Nd.data, // x
                  1,
                  given_vec.data, // y
                  1);
      
      // printf("force = %lf\n", force);
      // make_unit_vector(vec_boost_Nd);
      // matrix_mul(vec_boost_Nd, force);
      // cblas_daxpy(given_vec.size, 1.0, vec_boost_Nd.data, 1, given_vec.data, 1);
      // given_vec += vec_boost_Nd;
      
    }
  return 0;
}


long INTEGRATOR::EULER_ASSOCIATION::cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, long index_t, long given_index)
{
  given_vec.set_value(0.);
  MATRIX tmp_vec(TRAJ.Np, 1, 0.);
  for (long j=1; j<CONNECT.TOKEN[given_index]; j++)
    {
      long target_index = (long)CONNECT.HASH[given_index](j);
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, given_index, target_index, tmp_vec);
      double distance = tmp_vec.norm();
      double force = CONNECT.weight[given_index](j)*POTs.f_connector(distance, POTs.force_variables);
      make_unit_vector(tmp_vec);
      matrix_mul(tmp_vec, force);
      given_vec += tmp_vec;
    }
  return 0;
}

long INTEGRATOR::EULER::cal_repulsion_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t, long index_i)
{
  given_vec.set_value(0.);
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);

  for(long i=0; i<TRAJ.Np; i++)
    {
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, i, tmp_vec);
      double distance = tmp_vec.norm();
      if (i != index_i)
      // if ((i != index_i) && (i != GEOMETRY::get_connected_bead(TRAJ, index_i))) // this is temporally commented. need reinforcement
        {
          double repulsion = POTs.f_repulsion(distance, POTs.force_variables);
          make_unit_vector(tmp_vec);
          matrix_mul(tmp_vec, repulsion); 
          given_vec.add(tmp_vec);
        }
    }
  return 0;
}

long INTEGRATOR::EULER::cal_repulsion_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t, long index_i, MATRIX& vec_boost_Nd)
{
  given_vec.set_value(0.);
  // MATRIX tmp_vec(TRAJ.dimension, 1, 0.);

  for(long i=0; i<TRAJ.Np; i++)
    {
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, i, vec_boost_Nd);
      double distance = vec_boost_Nd.norm();
      if (i != index_i)
      // if ((i != index_i) && (i != GEOMETRY::get_connected_bead(TRAJ, index_i))) // this is temporally commented. need reinforcement
        {
          double repulsion = POTs.f_repulsion(distance, POTs.force_variables);
          make_unit_vector(vec_boost_Nd);
          matrix_mul(vec_boost_Nd, repulsion);
          given_vec.add(vec_boost_Nd);
        }
    }
  return 0;
}


long INTEGRATOR::EULER::cal_random_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t)
{
  RANDOM::single_random_vector_generator_variance(given_vec, 1.0);
  matrix_mul(given_vec, sqrt(2.0));
  POTs.scale_random(given_vec, POTs.force_variables);
  return 0;
}

long INTEGRATOR::EULER::simple_Euler(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, long index_t_now)
{
  long index_t_next = (index_t_now + 1)%TRAJ.Nt;
  TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;
  MATRIX force_spring(TRAJ.dimension, 1, 0.);
  MATRIX force_repulsion(TRAJ.dimension, 1, 0.);
  MATRIX force_random(TRAJ.dimension, 1, 0.);


  for (long i=0; i<TRAJ.Np; i++)
    {

      // from its print functionality, MATRIX objects are cleary working properly.
      // The reason is that with parallization, all the matrix object showed different address
      cal_connector_force(TRAJ, POTs, force_spring, index_t_now, i);
      cal_repulsion_force(TRAJ, POTs, force_repulsion, index_t_now, i);
      cal_random_force(TRAJ, POTs, force_random, index_t_now);
      for (long k=0; k<TRAJ.dimension; k++)
        {
          // printf("fc = %6.3e, fr = %6.3e, fR = %6.3e\n", force_spring(k), force_repulsion(k), force_random(k));
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*(force_spring(k) + force_repulsion(k)) + sqrt(TRAJ.dt)*force_random(k);
        }
    }
  return 0;
}



double ANALYSIS::cal_total_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, long index_t)
{
  return cal_potential_energy(TRAJ, POTs, index_t);// + cal_kinetic_energy(TRAJ, POTs, index_t);
}

double ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long index_t)
{
  double energy = 0.;

  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  for(long i=0; i<CONNECT.Np; i++)
    {
      for(long j=0; j<CONNECT.TOKEN[i]; j++)
        {
          energy += CONNECT.weight[i](j)*POTs.e_connector(GEOMETRY::get_minimum_distance(TRAJ, index_t, i, (long)CONNECT.HASH[i](j), tmp_vec), POTs.force_variables);
        }
    }
  return energy + ANALYSIS::cal_potential_energy(TRAJ, POTs, index_t);
}


double ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long index_t, MATRIX& tmp_vec)
{
  double energy = 0.;

  for(long i=0; i<CONNECT.Np; i++)
    {
      for(long j=0; j<CONNECT.TOKEN[i]; j++)
        {
          energy += CONNECT.weight[i](j)*POTs.e_connector(GEOMETRY::get_minimum_distance(TRAJ, index_t, i, (long)CONNECT.HASH[i](j), tmp_vec), POTs.force_variables);
        }
    }
  return energy + ANALYSIS::cal_potential_energy(TRAJ, POTs, index_t);
}


double ANALYSIS::cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, long index_t)
{
  double energy = 0.;

  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  for(long i=0; i<TRAJ.Np - 1; i++)
    {
      for(long j=i+1; j<TRAJ.Np; j++)
        {
          double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, i, j, tmp_vec);
          double U_ij = POTs.e_repulsion(distance, POTs.force_variables);
          energy += U_ij;
        }
    }
  return energy;
}

long ANALYSIS::cal_detail_repulsion(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, const char* fn, long index_t)
{
  ofstream FILE_ENERGY_INFO(fn, std::ios_base::app);
  double effective_range = POTs.force_variables[1];
  double range_limit = effective_range; // the previously setted 3 is temporal used, which is fixed now
  for(long i=0; i<TRAJ.Np -1; i++)
    {
      MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
      
      for(long j=i+1; j<TRAJ.Np; j++)
        {
          double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, i, j, tmp_vec);
          double U_ij = POTs.e_repulsion(distance, POTs.force_variables);
          if (distance <= range_limit)
            {
              FILE_ENERGY_INFO << std::scientific << TRAJ.c_t << '\t';
              FILE_ENERGY_INFO << i << '\t' << j << '\t';

              for(long k=0; k<TRAJ.dimension; k++)
                {
                  FILE_ENERGY_INFO << std::scientific << tmp_vec(k) << '\t';
                }
              FILE_ENERGY_INFO << std::scientific << distance << '\t' << U_ij << '\n';
            }
        }
    }
  FILE_ENERGY_INFO.close();
  return 0;
}

 
long ANALYSIS::CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& mat_energy, long index_t)
{
  mat_energy(0) = (TRAJ.c_t - 1)*TRAJ.dt;
  mat_energy(2) = cal_potential_energy(TRAJ, POTs, index_t);
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  // mat_energy(3) = cal_kinetic_energy(TRAJ, index_t);
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return 0;
}

long ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, long index_t, MATRIX& tmp_vec)
{
  mat_energy(0) = (TRAJ.c_t - 1)*TRAJ.dt;
  mat_energy(2) = ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy_boost(TRAJ, POTs, CONNECT, index_t, tmp_vec);
  
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return 0;
}

