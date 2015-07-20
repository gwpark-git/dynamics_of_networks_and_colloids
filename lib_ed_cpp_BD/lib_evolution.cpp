
#include "lib_evolution.h"

MKL_LONG ANALYSIS::GET_PDF_POTENTIAL_target(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG& index_particle, MKL_LONG& index_target, double& INDEX_PDF_U_ij, double& PDF_U_ij, MATRIX& vec_boost_ordered_pdf_ij)
{
  INDEX_PDF_U_ij = (double)index_target;
  double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, index_target, vec_boost_ordered_pdf_ij);
  PDF_U_ij = exp(-POTs.e_connector(distance, POTs.force_variables));
  return PDF_U_ij;
}


MKL_LONG ANALYSIS::GET_PDF_POTENTIAL(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG index_particle, MATRIX& INDEX_PDF_U, MATRIX& PDF_U, MATRIX *vec_boost_ordered_pdf)
{
// #pragma omp parallel for default(none) shared(TRAJ, INDEX_PDF_U, PDF_U, POTs, index_t, index_particle, vec_boost_ordered_pdf) schedule(static)
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, i, vec_boost_ordered_pdf[i]);
      INDEX_PDF_U(i) = i;
      PDF_U(i) = exp(-POTs.e_connector(distance, POTs.force_variables)); // that gave us probability
    }
  return 0;
}

MKL_LONG ANALYSIS::GET_ORDERED_PDF_POTENTIAL(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG index_particle, MATRIX& INDEX_PDF_U, MATRIX& PDF_U, MATRIX *vec_boost_ordered_pdf, double& time_PDF, double& time_SORT)
{
  double time_st = dsecnd();
  GET_PDF_POTENTIAL(TRAJ, index_t, POTs, index_particle, INDEX_PDF_U, PDF_U, vec_boost_ordered_pdf);
  double time_end_pdf = dsecnd();
  PDF_U.sort2(INDEX_PDF_U);
  double time_end_sort = dsecnd();
  time_PDF += time_end_pdf - time_st;
  time_SORT += time_end_sort - time_end_pdf;
  return 0;
}

MKL_LONG INTEGRATOR::EULER::cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index)
{
  given_vec.set_value(0.);
  return 0;
}

MKL_LONG INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index, MATRIX vec_boost_Nd)
{
  given_vec.set_value(0.);
  // MATRIX tmp_vec(TRAJ.Np, 1, 0.);
// #pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, given_vec, index_t, given_index, vec_boost_Nd) schedule(dynamic)
  for (MKL_LONG j=1; j<CONNECT.TOKEN(given_index); j++)
    {
      MKL_LONG target_index = CONNECT.HASH(given_index, j);
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, given_index, target_index, vec_boost_Nd);
      double distance = vec_boost_Nd.norm();
      double force = (double)CONNECT.weight(given_index,j)*POTs.f_connector(distance, POTs.force_variables);
      // printf("force = %lf\n", force);
      make_unit_vector(vec_boost_Nd);
      matrix_mul(vec_boost_Nd, force);
      given_vec += vec_boost_Nd;
    }
  return 0;
}


MKL_LONG INTEGRATOR::EULER_ASSOCIATION::cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index)
{
  given_vec.set_value(0.);
  MATRIX tmp_vec(TRAJ.Np, 1, 0.);
// #pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, given_vec, index_t, given_index) firstprivate(tmp_vec) 
  for (MKL_LONG j=1; j<CONNECT.TOKEN(given_index); j++)
    {
      MKL_LONG target_index = CONNECT.HASH(given_index, j);
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, given_index, target_index, tmp_vec);
      double distance = tmp_vec.norm();
      double force = (double)CONNECT.weight(given_index,j)*POTs.f_connector(distance, POTs.force_variables);
      // printf("force = %lf\n", force);
      make_unit_vector(tmp_vec);
      matrix_mul(tmp_vec, force);
      given_vec += tmp_vec;
    }
  return 0;
}

MKL_LONG INTEGRATOR::EULER::cal_repulsion_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG index_i)
{
  given_vec.set_value(0.);
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);

// #pragma omp parallel for default(none) shared(TRAJ, POTs, index_t, index_i, given_vec) firstprivate(tmp_vec)
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, i, tmp_vec);
      double distance = tmp_vec.norm();
      if (i != index_i)
      // if ((i != index_i) && (i != GEOMETRY::get_connected_bead(TRAJ, index_i))) // this is temporally commented. need reinforcement
        {
          double repulsion = POTs.f_repulsion(distance, POTs.force_variables);
          // printf("repulsion::repulsion(distance(t, i, j)) = r(d(%6.3e, %ld, %ld)) = r(%6.3e) = %6.3e\n", TRAJ(index_t, 0), index_i, i, distance, repulsion);
          make_unit_vector(tmp_vec);
          matrix_mul(tmp_vec, repulsion);
          given_vec.add(tmp_vec);
        }
    }
  return 0;
}

MKL_LONG INTEGRATOR::EULER::cal_repulsion_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG index_i, MATRIX vec_boost_Nd)
{
  given_vec.set_value(0.);
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);

// #pragma omp parallel for default(none) shared(TRAJ, POTs, index_t, index_i, given_vec, vec_boost_Nd) schedule(static)
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, i, vec_boost_Nd);
      double distance = vec_boost_Nd.norm();
      if (i != index_i)
      // if ((i != index_i) && (i != GEOMETRY::get_connected_bead(TRAJ, index_i))) // this is temporally commented. need reinforcement
        {
          double repulsion = POTs.f_repulsion(distance, POTs.force_variables);
          // printf("repulsion::repulsion(distance(t, i, j)) = r(d(%6.3e, %ld, %ld)) = r(%6.3e) = %6.3e\n", TRAJ(index_t, 0), index_i, i, distance, repulsion);
          make_unit_vector(vec_boost_Nd);
          matrix_mul(vec_boost_Nd, repulsion);
          given_vec.add(vec_boost_Nd);
        }
    }
  return 0;
}


MKL_LONG INTEGRATOR::EULER::cal_random_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t)
{
  RANDOM::single_random_vector_generator_variance(given_vec, 1.0);
  matrix_mul(given_vec, sqrt(2.0));
  POTs.scale_random(given_vec, POTs.force_variables);
  return 0;
}

// MKL_LONG INTEGRATOR::EULER_ASSOCIATION::simple_Euler(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_t_now)
// {
  

// }

MKL_LONG INTEGRATOR::EULER::simple_Euler(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t_now)
{
  MKL_LONG index_t_next = (index_t_now + 1)%TRAJ.Nt;
  TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;
  // std::cout << TRAJ(index_t_next) << std::endl;
  // TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;
  MATRIX force_spring(TRAJ.dimension, 1, 0.);
  MATRIX force_repulsion(TRAJ.dimension, 1, 0.);
  MATRIX force_random(TRAJ.dimension, 1, 0.);


#pragma omp parallel for default(none) shared(TRAJ, POTs, index_t_now, index_t_next) firstprivate(force_spring, force_repulsion, force_random) schedule(static) // firstprivate called copy-constructor while private called default constructor
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {

      // from its print functionality, MATRIX objects are cleary working properly.
      // The reason is that with parallization, all the matrix object showed different address
      cal_connector_force(TRAJ, POTs, force_spring, index_t_now, i);
      cal_repulsion_force(TRAJ, POTs, force_repulsion, index_t_now, i);
      cal_random_force(TRAJ, POTs, force_random, index_t_now);
      for (MKL_LONG k=0; k<TRAJ.dimension; k++)
        {
          // printf("fc = %6.3e, fr = %6.3e, fR = %6.3e\n", force_spring(k), force_repulsion(k), force_random(k));
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*(force_spring(k) + force_repulsion(k)) + sqrt(TRAJ.dt)*force_random(k);
        }
    }
  return 0;
}

// MKL_LONG INTEGRATOR::EULER::simple_Euler_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t_now, MATRIX *vec_boost_Nd_connector, MATRIX *vec_boost_Nd_repulsion, MATRIX *vec_boost_Nd_random)
// {
//   MKL_LONG index_t_next = (index_t_now + 1)%TRAJ.Nt;
//   TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;
//   // std::cout << TRAJ(index_t_next) << std::endl;
//   // TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;
//   MATRIX force_spring(TRAJ.dimension, 1, 0.);
//   MATRIX force_repulsion(TRAJ.dimension, 1, 0.);
//   MATRIX force_random(TRAJ.dimension, 1, 0.);

// #pragma omp parallel for default(none) shared(TRAJ, POTs, index_t_now, index_t_next, vec_boost_Nd_connector, vec_boost_Nd_repulsion, vec_boost_Nd_random)
//   // note that becase of overhead for the firstprivate, everything included on the shared
//   // // firstprivate called copy-constructor while private called default constructor
//   for (MKL_LONG i=0; i<TRAJ.Np; i++)
//     {
//       // from its print functionality, MATRIX objects are cleary working properly.
//       // The reason is that with parallization, all the matrix object showed different address
//       cal_connector_force(TRAJ, POTs, force_spring, index_t_now, i);
//       cal_repulsion_force(TRAJ, POTs, force_repulsion, index_t_now, i);
//       cal_random_force(TRAJ, POTs, force_random, index_t_now);
//       for (MKL_LONG k=0; k<TRAJ.dimension; k++)
//         {
//           // printf("fc = %6.3e, fr = %6.3e, fR = %6.3e\n", force_spring(k), force_repulsion(k), force_random(k));
//           TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*(force_spring(k) + force_repulsion(k)) + sqrt(TRAJ.dt)*force_random(k);
//         }
//     }
//   return 0;
// }




double ANALYSIS::cal_total_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t)
{
  return cal_potential_energy(TRAJ, POTs, index_t);// + cal_kinetic_energy(TRAJ, POTs, index_t);
}

double ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_t)
{
  double energy = 0.;

  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  for(MKL_LONG i=0; i<CONNECT.Np; i++)
    {
      for(MKL_LONG j=0; j<CONNECT.TOKEN(i); j++)
        {
          energy += (double)CONNECT.weight(i,j)*POTs.e_connector(GEOMETRY::get_minimum_distance(TRAJ, index_t, i, CONNECT.HASH(i,j), tmp_vec), POTs.force_variables);
        }
    }
  return energy + ANALYSIS::cal_potential_energy(TRAJ, POTs, index_t);
}

double ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_t, MATRIX *vec_boost_Nd)
{
  double energy = 0.;

  // MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
#pragma omp parallel for default(none) shared(energy, TRAJ, POTs, CONNECT, index_t, vec_boost_Nd) schedule(static)
  for(MKL_LONG i=0; i<CONNECT.Np; i++)
    {
      for(MKL_LONG j=0; j<CONNECT.TOKEN(i); j++)
        {
          energy += (double)CONNECT.weight(i,j)*POTs.e_connector(GEOMETRY::get_minimum_distance(TRAJ, index_t, i, CONNECT.HASH(i,j), vec_boost_Nd[i]), POTs.force_variables);
        }
    }
  return energy + ANALYSIS::cal_repulsion_energy_boost(TRAJ, POTs, index_t, vec_boost_Nd);
}

double ANALYSIS::cal_repulsion_energy_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t, MATRIX *vec_boost_Nd)
{
  double energy = 0.;

  // MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
// #pragma omp parallel for default(none) shared(energy, TRAJ, POTs, index_t) firstprivate(tmp_vec)
#pragma omp parallel for default(none) shared(energy, TRAJ, POTs, index_t, vec_boost_Nd) schedule(static)
  for(MKL_LONG i=0; i<TRAJ.Np - 1; i++)
    {
      for(MKL_LONG j=i+1; j<TRAJ.Np; j++)
        {
          double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, i, j, vec_boost_Nd[i]);
          double U_ij = POTs.e_repulsion(distance, POTs.force_variables);
          energy += U_ij;
        }
    }
  return energy;
}


double ANALYSIS::cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t)
{
  double energy = 0.;

  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
// #pragma omp parallel for default(none) shared(energy, TRAJ, POTs, index_t) firstprivate(tmp_vec)
  for(MKL_LONG i=0; i<TRAJ.Np - 1; i++)
    {
      for(MKL_LONG j=i+1; j<TRAJ.Np; j++)
        {
          double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, i, j, tmp_vec);
          double U_ij = POTs.e_repulsion(distance, POTs.force_variables);
          energy += U_ij;
        }
    }
  return energy;
}

MKL_LONG ANALYSIS::cal_detail_repulsion(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, const char* fn, MKL_LONG index_t)
{
  ofstream FILE_ENERGY_INFO(fn, std::ios_base::app);
  // double effective_range = TRAJ.force_variables[1];
  double range_limit = 3.;
  for(MKL_LONG i=0; i<TRAJ.Np -1; i++)
    {
      MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
      
      for(MKL_LONG j=i+1; j<TRAJ.Np; j++)
        {
          double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t, i, j, tmp_vec);
          double U_ij = POTs.e_repulsion(distance, POTs.force_variables);
          if (distance <= range_limit)
          // if (U_ij != 0.)
            {
              FILE_ENERGY_INFO << std::scientific << TRAJ.c_t << '\t';
              FILE_ENERGY_INFO << i << '\t' << j << '\t';
              // std::cout << index_t << '\t' << i << '\t' << j << '\t';

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

 
MKL_LONG ANALYSIS::CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& mat_energy, MKL_LONG index_t)
{
  // index_t = index_t%TRAJ.Nt; // no more needed
  mat_energy(0) = (TRAJ.c_t - 1)*TRAJ.dt;
  mat_energy(2) = cal_potential_energy(TRAJ, POTs, index_t);
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  // mat_energy(3) = cal_kinetic_energy(TRAJ, index_t);
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return 0;
}

MKL_LONG ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, MKL_LONG index_t)
{
  // index_t = index_t%TRAJ.Nt; // no more needed
  mat_energy(0) = (TRAJ.c_t - 1)*TRAJ.dt;
  mat_energy(2) = ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy(TRAJ, POTs, CONNECT, index_t);
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  // mat_energy(3) = cal_kinetic_energy(TRAJ, index_t);
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return 0;
}

