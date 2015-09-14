#include <iostream>
#include "lib_ed_cpp_BD/matrix_ed.h"
#include "lib_ed_cpp_BD/lib_traj.h"
#include "lib_ed_cpp_BD/lib_evolution.h"
#include "lib_ed_cpp_BD/lib_association.h"
#include "lib_ed_cpp_BD/lib_potential.h"
// #include "lib_ed_cpp_BD/potential_definition.h"
#include <string>
// #include <istream>
using namespace std;

int help()
{
  cout << "USAGE of Brownian Dynamics with Quadratic Repulsive Potential\n";
  cout << "argv[1] == condition files including variables:\n";
  cout << "\tLONG INTEGER format: N_dimension, Nt, Np\n";
  cout << "\tDOUBLE FLOAT format: dt, repulsion_coefficient, effective_distance\n";
  cout << "\tSTRING format: filename_trajectory, filename_energy, filename_energy_detail\n";
  return 0;
}

MKL_LONG main_PURE_BROWNIAN(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, COND& given_condition);
MKL_LONG main_NAPLE_REPULSION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, COND& given_condition);
MKL_LONG main_NAPLE_ASSOCIATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, COND& given_condition);
// MKL_LONG main_NAPLE_ASSOCIATION_TEST(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, COND& given_condition);

int main(int argc, char* argv[])
{
  if(argc==1)
    {
      help();
      return 0;
    }
  else
    {
      COND given_condition(argv[1]);
      given_condition.cond_print();

      MKL_LONG N_basic = 2;
      if (given_condition("Integrator") == "Euler") 
        {
          // at the moment, the only method available is simple Euler integrator
          // for this reason, this condition will not be varied during condition test
          // however, the explicit description should be needed for condition data
          // that is because to prevent potential errors when we implement higher order methods.
          N_basic = 2;
        }

      TRAJECTORY TRAJ(given_condition, N_basic);
      GENERATOR::random_position_generator(TRAJ);

      // // MATRIX energy_info(Np*((MKL_LONG)Np/2 + 1), 2*TRAJ.dimension + 2, 0.);
      // if (given_condition("INITIAL")=="LOAD")
      //   {
      //     // should load trajectory file
      //   }
      // else
      //   {
      //   }
      

      POTENTIAL_SET POTs;
      if(given_condition("Method") == "NAPLE_REPULSION")
        {
          FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_set(POTs, given_condition);
          main_NAPLE_REPULSION(TRAJ, POTs, given_condition);
        }
      // else if(given_condition("Method") == "NAPLE_ASSOCIATION_TEST")
      //   {
      //     MKL_LONG Np = atol(given_condition("Np").c_str());
      //     MKL_LONG N_chains_per_particle = atol(given_condition("N_chains_per_particle").c_str());
      //     MKL_LONG TOL_connection = atol(given_condition("tolerance_allowing_connections").c_str());
      //     bool allowing_multiple_connections = FALSE;
      //     if (given_condition("allowing_multiple_connections") == "TRUE")
      //       allowing_multiple_connections = TRUE;
      //     ASSOCIATION CONNECT(Np, N_chains_per_particle, TOL_connection, allowing_multiple_connections);
      //     FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition);
      //     main_NAPLE_ASSOCIATION_TEST(TRAJ, POTs, CONNECT, given_condition);
      //   }
      else if(given_condition("Method") == "NAPLE_ASSOCIATION")
        {
          // MKL_LONG Np = atol(given_condition("Np").c_str());
          // MKL_LONG N_chains_per_particle = atol(given_condition("N_chains_per_particle").c_str());
          // MKL_LONG TOL_connection = atol(given_condition("tolerance_allowing_connections").c_str());
          // bool allowing_multiple_connections = FALSE;
          // if (given_condition("allowing_multiple_connections") == "TRUE")
          //   allowing_multiple_connections = TRUE;
          // ASSOCIATION CONNECT(Np, N_chains_per_particle, TOL_connection, allowing_multiple_connections);
          ASSOCIATION CONNECT(TRAJ, given_condition);
          FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition);
          main_NAPLE_ASSOCIATION(TRAJ, POTs, CONNECT, given_condition);
        }
      else if(given_condition("Method") == "PURE_BROWNIAN")
        {
          FORCE::DEFAULT::EMPTY_force_set(POTs, given_condition);
          main_PURE_BROWNIAN(TRAJ, POTs, given_condition);
        }
    }
  return 0;
}

MKL_LONG main_PURE_BROWNIAN(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, COND& given_condition)
{
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_trajectory")).c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_energy")).c_str();
  string filename_energy_info  = (given_condition("output_path") + '/' + given_condition("filename_energy_info")).c_str();

  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 4, 0.);
  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  TRAJ.fprint_row(filename_trajectory.c_str(), 0);
  MKL_LONG N_basic = TRAJ.rows;
  for(MKL_LONG t=0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      INTEGRATOR::EULER::simple_Euler(TRAJ, POTs, index_t_now); 
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next);
      ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, index_t_next);
      if(t%N_skip==0)
        {
          printf("STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_next), energy(1));
          TRAJ.fprint_row(filename_trajectory.c_str(), index_t_next);
          energy.fprint(filename_energy.c_str());
        }
    }
  return 0;
}

MKL_LONG main_NAPLE_REPULSION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, COND& given_condition)
{
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_trajectory")).c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_energy")).c_str();
  string filename_energy_info  = (given_condition("output_path") + '/' + given_condition("filename_energy_info")).c_str();

  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 4, 0.);
  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  TRAJ.fprint_row(filename_trajectory.c_str(), 0);

  MKL_LONG N_basic = TRAJ.rows;
  for(MKL_LONG t=0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      INTEGRATOR::EULER::simple_Euler(TRAJ, POTs, index_t_now); 
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next);
      ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, index_t_next);

      if(t%N_skip==0)
        {
          printf("STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_next), energy(1));
          TRAJ.fprint_row(filename_trajectory.c_str(), index_t_next);
          energy.fprint(filename_energy.c_str());
        }
      if(t%N_energy_frequency==0)
        {
          ANALYSIS::cal_detail_repulsion(TRAJ, POTs, filename_energy_info.c_str(), index_t_next);
        }
    }
  return 0;
}


MKL_LONG main_NAPLE_ASSOCIATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, COND& given_condition)
{
  double time_st_simulation = dsecnd();
  MKL_LONG N_steps_block = atol(given_condition("N_steps_block").c_str());
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_trajectory")).c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_energy")).c_str();
  string filename_energy_info  = (given_condition("output_path") + '/' + given_condition("filename_energy_info")).c_str();
  string filename_HASH = (given_condition("output_path") + '/' + given_condition("filename_HASH")).c_str();
  string filename_CASE = (given_condition("output_path") + '/' + given_condition("filename_CASE")).c_str();
  string filename_weight = (given_condition("output_path") + '/' + given_condition("filename_weight")).c_str();
  string filename_MC_LOG = (given_condition("output_path") + '/' + given_condition("filename_MC_LOG")).c_str();
  ofstream FILE_LOG;

  MKL_LONG N_THREADS = atol(given_condition("N_THREADS").c_str());
  omp_set_num_threads(N_THREADS);

  double time_MC = 0.;
  double time_LV = 0.;
  double time_AN = 0.;
  double time_file = 0.;

  // the following are the boosting the allocation and dislocation of the MATRIX classes
  MKL_LONG cond_longer = TRAJ.Np > CONNECT.N_max ? TRAJ.Np : CONNECT.N_max;
  MATRIX *vec_boost_Nd_parallel = new MATRIX [cond_longer];
  MATRIX *vec_boost_Np_parallel = new MATRIX [cond_longer];
  MATRIX **vec_boost_Nd_Np_parallel_connector = new MATRIX* [cond_longer];
  MATRIX **vec_boost_Nd_Np_parallel_repulsion = new MATRIX* [cond_longer];
  
  for(MKL_LONG i=0; i<cond_longer; i++)
    {
      vec_boost_Nd_parallel[i].initial(TRAJ.dimension, 1, 0.);
      vec_boost_Np_parallel[i].initial(TRAJ.Np, 1, 0.);
      vec_boost_Nd_Np_parallel_connector[i] = new MATRIX [TRAJ.Np];
      vec_boost_Nd_Np_parallel_repulsion[i] = new MATRIX [TRAJ.Np];
      for(MKL_LONG j=0; j<TRAJ.Np; j++)
        {
          vec_boost_Nd_Np_parallel_connector[i][j].initial(TRAJ.dimension, 1, 0.);
          vec_boost_Nd_Np_parallel_repulsion[i][j].initial(TRAJ.dimension, 1, 0.);
        }
    }
  MATRIX *force_spring = new MATRIX [TRAJ.Np];
  MATRIX *force_repulsion = new MATRIX [TRAJ.Np];
  MATRIX *force_random = new MATRIX [TRAJ.Np];
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      force_spring[i].initial(TRAJ.dimension, 1, 0.);
      force_repulsion[i].initial(TRAJ.dimension, 1, 0.);
      force_random[i].initial(TRAJ.dimension, 1, 0.);
    }
  // MATRIX *vec_boost_Np_parallel_Np = new MATRIX [TRAJ.Np];
  // MATRIX *vec_boost_Np_parallel_Token = new MATRIX [CONNECT.N_max];
  // if (given_condition("MC_LOG") == "TRUE")
  //   {
  //     MKL_LONG cnt_log = 0;
  //     FILE_LOG.open(filename_MC_LOG.c_str(), std::ios_base::out);
  //     fprintf(FILE_LOG, "MC_LOG FILE, column information:\n");
  //     fprintf(FILE_LOG, "%02d, steps: time\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, steps: MC\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, index: itself\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, index: target\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, index: new target\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, hash: itself\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, hash: target\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, hash: new target (k-value)\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, roll: dCDF\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, roll: dCDF_U", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, stat: N_associations\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, stat: functionality (N_association/N_beads)\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, stat: TOKEN (N_diff_connections)\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, stat: N_chain_ends\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, count: add\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, count: del\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, count: move\n", cnt_log++);
  //     fprintf(FILE_LOG, "%02d, count: cancel\n", cnt_log++);
  //   }
  if(given_condition("MC_LOG") == "TRUE")
    {
      FILE_LOG.open(filename_MC_LOG.c_str(), std::ios_base::out);
      FILE_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dCDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << endl;
    }
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 4, 0.);
  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  // TRAJ.fprint_row(filename_trajectory.c_str(), 0);

  MKL_LONG N_basic = TRAJ.rows;

  // double bMSE_tolerance = atof(given_condition("tolerance_MSE_association").c_str());
  double tolerance_association = atof(given_condition("tolerance_association").c_str());
  MKL_LONG N_max_steps = atol(given_condition("N_max_steps").c_str());

  MATRIX *dCDF_U = new MATRIX [TRAJ.Np];
  MATRIX *INDEX_dCDF_U = new MATRIX [TRAJ.Np];
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      dCDF_U[i].initial(TRAJ.Np, 1, 0.);
      INDEX_dCDF_U[i].initial(TRAJ.Np, 1, 0);
    }
  // MATRIX dCDF_U(TRAJ.Np, 1, 0.);
  // MATRIX INDEX_dCDF_U(TRAJ.Np, 1, 0);

  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  CONNECT.initial();
  for(MKL_LONG i=0; i<CONNECT.Np; i++)
    CONNECT.TOKEN(i) = 1;

  double dt_1 = 0., dt_2 = 0., dt_3 = 0., dt_4 = 0., dt_5 = 0., dt_6 = 0., dt_7 = 0.;
  double dt_det_pdf = 0.;
  double dt_pdf = 0., dt_sort = 0.;
  
  double time_MC_1 = 0., time_MC_2 = 0., time_MC_3 = 0., time_MC_4 = 0., time_MC_5 = 0., time_MC_6 = 0., time_MC_7 = 0., time_MC_8 = 0.;
  const gsl_rng_type *T_boost;
  gsl_rng *r_boost;
  gsl_rng_env_setup();
  T_boost = gsl_rng_default;
  r_boost = gsl_rng_alloc(T_boost);
  

  for(MKL_LONG t=0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      TRAJ(index_t_next) = (++TRAJ.c_t) * TRAJ.dt;

      
      tmp_vec.set_value(0.);

      // TRAJ(index_t_next)++;
      
      // while(fabs(block_MSE_next - block_MSE_now) > bMSE_tolerance || ++cnt < N_max_steps)
      MKL_LONG IDENTIFIER_ASSOC = TRUE;
      double max_try_ASSOC = tolerance_association;
      double N_diff = 0.;
      double max_N_diff = 0.;
      MKL_LONG sum_over_MC_steps = 0;
      MKL_LONG count_M = 0;
      MKL_LONG pre_N_associations = 0;
      MKL_LONG N_associations = 0;
      // double max_block_MSE = bMSE_tolerance;
      
      // double block_MSE_now = 0., block_MSE_next = 10000.;
      // double block_MSE_square_mean = 0., block_MSE_mean_square = 0.;
      MKL_LONG cnt = 1;
      if(given_condition("MC_renewal")=="TRUE")
        {
          CONNECT.initial();
          for(MKL_LONG i=0; i<CONNECT.Np; i++)
            CONNECT.TOKEN(i) = 1;
        }

      MKL_LONG cnt_del = 0;
      MKL_LONG cnt_add = 0;
      MKL_LONG cnt_mov = 0;
      MKL_LONG cnt_cancel = 0;
      double time_st_MC = dsecnd();
      gsl_rng_set(r_boost, random());

      // computing the distance map between beads for reducing overhead
      // note that the computing distance map during association spend 80% of computing time
      // That involve computing map from 1 to 3 beads (with number of Monte-Carlo steps)
      // In this case, however, we compute all the distance map, once, then uses this information throughout the Monte-Carlo steps since the configurational information does not change during Monte-Carlo stpes.
      // // Note that NO upper triangle form is generated. All the i-j distance is computed twice. 
      // // Therefore, further reduce is possible to modify this computing.
      // // On here, the feature is ignored for convenience.
      double time_st_det_pdf = dsecnd();
#pragma omp parallel for shared(TRAJ, index_t_now, vec_boost_Nd_parallel, INDEX_dCDF_U, dCDF_U, dt_pdf, dt_sort) num_threads(N_THREADS)
      for(MKL_LONG i=0; i<TRAJ.Np; i++)
        {
          double time_st_pdf = dsecnd();
          ANALYSIS::GET_dCDF_POTENTIAL(TRAJ, index_t_now, POTs, i, INDEX_dCDF_U[i], dCDF_U[i], vec_boost_Nd_parallel[i]);
          double time_end_pdf = dsecnd();
          dCDF_U[i].sort2(INDEX_dCDF_U[i]);
          double time_end_sort = dsecnd();
          dt_pdf += time_end_pdf - time_st_pdf;
          dt_sort += time_end_sort - time_end_pdf;
        }
      double time_end_det_pdf = dsecnd();
      dt_det_pdf += time_end_det_pdf - time_st_det_pdf;

      // if (given_condition("STATUS")=="FULL")
      //   {
      while(IDENTIFIER_ASSOC && ++cnt < N_max_steps)
        {
          // choice for visiting bead   
          MKL_LONG total_chain_ends = CONNECT.N_TOTAL_CONNECTED_ENDS();
          // printf("chain ends attached to beads, %ld, per beads, %ld, (add, mov, del) = (%ld, %ld, %ld)\n", total_chain_ends, total_chain_ends/TRAJ.Np, cnt_add, cnt_mov, cnt_del);
          time_MC_1 = dsecnd();
          // MKL_LONG index_itself = RANDOM::return_LONG_INT_rand(TRAJ.Np);
          MKL_LONG index_itself = RANDOM::return_LONG_INT_rand_boost(r_boost, TRAJ.Np);
          // choice for selected chain end
          double rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(r_boost);
          time_MC_2 = dsecnd();
          // MKL_LONG HASH_INDEX = CONNECT.GET_HASH_FROM_ROLL(index_itself, rolling_P);
          MKL_LONG index_hash_selected_chain = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF);
          MKL_LONG index_other_end_of_selected_chain = CONNECT.HASH(index_itself, index_hash_selected_chain);
          time_MC_3 = dsecnd();
          // choice for behaviour of selected chain end
          double rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(r_boost);
          // the PDF is already computed in the previous map
          // ANALYSIS::GET_ORDERED_PDF_POTENTIAL(TRAJ, index_t_now, POTs, index_other_end_of_selected_chain, INDEX_dCDF_U, dCDF_U, vec_boost_Nd_parallel, dt_pdf, dt_sort);
          time_MC_4 = dsecnd();
          MKL_LONG k=CONNECT.Np-1;
          for(k=CONNECT.Np-1; k>=0; k--)
            {
              // printf("(i, k, P) = (%ld, %ld, %lf)\n", index_other_end_of_selected_chain, k, dCDF_U[index_other_end_of_selected_chain](k));
              if(dCDF_U[index_other_end_of_selected_chain](k) < rolling_dCDF_U)
                {
                  k++;
                  break;
                }
            }
          MKL_LONG index_new_end_of_selected_chain = INDEX_dCDF_U[index_other_end_of_selected_chain](k);
          time_MC_5 = dsecnd();
          // for(MKL_LONG i=CONNECT.Np-10; i<CONNECT.Np; i++)
          //   {
          //     printf("(i, I(i), dPDF(i)) = (%ld, %ld, %6.3f), (k, roll) = (%ld, %6.3f)\n", i, (MKL_LONG) INDEX_dCDF_U(i), dCDF_U(i), k, rolling_dCDF_U);
          //   }

          // printf("(cnt, i_V, i_H, i_T) = (%ld, %ld, %ld, %ld),\t (k, i_k, N_c, TOKEN, CASE) = (%ld, %ld, %ld, %ld, %6.3e)\t(TOL=%lf)\n", cnt, index_itself, index_hash_selected_chain, index_other_end_of_selected_chain, k, index_new_end_of_selected_chain, (MKL_LONG)CONNECT.CASE(index_itself, 0), CONNECT.TOKEN(index_itself), CONNECT.CASE(index_itself, CONNECT.TOKEN(index_itself)-1), fabs(block_MSE_next-block_MSE_now));
          // The following codes are of importance to modification
          // the subject and object for the action is replaced
          // if(CONNECT.DEL_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
          //   {
          //     CONNECT.del_association_hash(index_itself, index_hash_selected_chain);
          //     cnt_del ++;
          //   }
          // else if (CONNECT.NEW_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
          //   {
          //     CONNECT.add_association_INFO(POTs, index_itself, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_itself, index_new_end_of_selected_chain, tmp_vec));
          //     cnt_add ++;
          //   }
          // else if (CONNECT.MOV_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
          //   {
          //     CONNECT.del_association_hash(index_itself, index_hash_selected_chain);
          //     CONNECT.add_association_INFO(POTs, index_itself, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_itself, index_new_end_of_selected_chain, tmp_vec));
          //     cnt_mov ++;
          //   }
          // else
          //   {
          //     cnt_cancel ++;
          //   }
          MKL_LONG index_hash_itself_from_other_end = CONNECT.FIND_HASH_INDEX(index_other_end_of_selected_chain, index_itself);
          if(CONNECT.DEL_ASSOCIATION(CONNECT, index_other_end_of_selected_chain, index_itself, index_new_end_of_selected_chain))
            {
              // CONNECT.del_association_hash(index_itself, index_hash_selected_chain);
              CONNECT.del_association_hash(index_other_end_of_selected_chain, index_hash_itself_from_other_end);
              cnt_del ++;
            }
          else if (CONNECT.NEW_ASSOCIATION(CONNECT, index_other_end_of_selected_chain, index_itself, index_new_end_of_selected_chain))
            {
              // CONNECT.add_association_INFO(POTs, index_itself, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_itself, index_new_end_of_selected_chain, tmp_vec));
              CONNECT.add_association_INFO(POTs, index_other_end_of_selected_chain, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_other_end_of_selected_chain, index_new_end_of_selected_chain, tmp_vec));
              
              cnt_add ++;
            }
          else if (CONNECT.MOV_ASSOCIATION(CONNECT, index_other_end_of_selected_chain, index_itself, index_new_end_of_selected_chain))
            {
              CONNECT.del_association_hash(index_other_end_of_selected_chain, index_hash_itself_from_other_end);
              CONNECT.add_association_INFO(POTs, index_other_end_of_selected_chain, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_other_end_of_selected_chain, index_new_end_of_selected_chain, tmp_vec));
              cnt_mov ++;
            }
          else
            {
              cnt_cancel ++;
            }
          
          time_MC_6 = dsecnd();
          // MKL_LONG k=0;
          MKL_LONG index_set[3] = {index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain};
          // for(MKL_LONG k=0; k<CONNECT.Np; k++)
          for(MKL_LONG i=0; i<3; i++)
            {
              CONNECTIVITY_update_Z_particle(CONNECT, index_set[i]);
              CONNECTIVITY_update_dPDF_particle(CONNECT, index_set[i]);
              CONNECTIVITY_update_dCDF_particle(CONNECT, index_set[i]);
            }
          
          time_MC_7 = dsecnd();

          pre_N_associations = N_associations;
          N_associations = CONNECT.N_TOTAL_ASSOCIATION()/2.;

          if(N_associations != pre_N_associations)
            {
              count_M ++;
              if (count_M%N_steps_block == 0)
                // if (count_M > 1)
                {
                  {
                    N_diff = fabs((double)(1./(count_M*(count_M-1)))*((count_M-1)*N_associations - sum_over_MC_steps));
                    max_N_diff = max_N_diff > N_diff ? max_N_diff : N_diff;
                    if(N_diff/max_N_diff < tolerance_association)
                      {
                        IDENTIFIER_ASSOC = FALSE;
                      }
                    // printf("(count, N_ass, sum_pre, N_diff, max_N_diff, ratio) = (%6ld, %6ld, %6ld, %6.4f, %6.4f, %6.4f)\n", count_M, N_associations, sum_over_MC_steps, N_diff, max_N_diff, N_diff/max_N_diff);
                  }
                }
              sum_over_MC_steps += N_associations;
            }

          
          // block_MSE_square_mean += pow(N_associations,2.);
          // block_MSE_mean_square += N_associations;
          // if(cnt%N_steps_block == 0 && cnt != 0)
          //   {
          //     block_MSE_now = block_MSE_next;
          //     block_MSE_next = fabs(block_MSE_square_mean - pow(block_MSE_mean_square,2.))/(float)N_steps_block;
          //     max_block_MSE = max_block_MSE > block_MSE_next ? max_block_MSE : block_MSE_next;
          //     // printf("bMSE = %6.3e\t bMSE/bMSE_max = %6.3e, maximum = %6.3e\n", block_MSE_next, block_MSE_next/max_block_MSE, max_block_MSE);
          //     // printf("bMSE_now = %6.3e\tbMSE_next = %6.3e\t diff_bMSE = %6.3e\n", fabs(block_MSE_now), fabs(block_MSE_next), fabs(block_MSE_next - block_MSE_now)/float(N_steps_block));
          //     if (block_MSE_next/max_block_MSE < bMSE_tolerance)
          //       {
          //         IDENTIFIER_ASSOC = FALSE;
          //       }
          //     block_MSE_square_mean = 0.;
          //     block_MSE_mean_square = 0.;
          //   }
          time_MC_8 = dsecnd();
          dt_1 += time_MC_2 - time_MC_1;
          dt_2 += time_MC_3 - time_MC_2;
          dt_3 += time_MC_4 - time_MC_3;
          dt_4 += time_MC_5 - time_MC_4;
          dt_5 += time_MC_6 - time_MC_5;
          dt_6 += time_MC_7 - time_MC_6;
          dt_7 += time_MC_8 - time_MC_7;
          if (given_condition("MC_LOG") == "TRUE")
            {
              MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
              FILE_LOG << cnt << '\t' << index_itself << '\t' << rolling_dCDF<< '\t'  << index_hash_selected_chain<< '\t'  << index_other_end_of_selected_chain<< '\t'  << rolling_dCDF_U<< '\t'  << k<< '\t'  << index_new_end_of_selected_chain<< '\t'  << CONNECT.TOKEN(index_itself)<< '\t'<<CONNECT.N_CONNECTED_ENDS(index_itself) << '\t' << CONNECT.weight(index_itself, 0) <<'\t' <<  total_bonds << '\t'  << cnt_add<< '\t'  << cnt_mov<< '\t'  << cnt_del<< '\t'  << cnt_cancel << endl;
            }

          //   }
        } // while
      double time_end_MC = dsecnd();

      // MATRIX force_spring(TRAJ.dimension, 1, 0.);
      // MATRIX force_repulsion(TRAJ.dimension, 1, 0.);
      // MATRIX force_random(TRAJ.dimension, 1, 0.);

// #pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next) firstprivate(force_spring, force_repulsion, force_random) // firstprivate called copy-constructor while private called default constructor
#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next, vec_boost_Nd_parallel, force_spring, force_repulsion, force_random) schedule(static) num_threads(N_THREADS)
      for (MKL_LONG i=0; i<TRAJ.Np; i++)
        {
          force_spring[i].set_value(0);
          force_repulsion[i].set_value(0);
          force_random[i].set_value(0);
          // from its print functionality, MATRIX objects are cleary working properly.
          // The reason is that with parallization, all the matrix object showed different address
          INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(TRAJ, POTs, CONNECT, force_spring[i], index_t_now, i, vec_boost_Nd_parallel[i]);
          INTEGRATOR::EULER::cal_repulsion_force_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, vec_boost_Nd_parallel[i]);
          INTEGRATOR::EULER::cal_random_force(TRAJ, POTs, force_random[i], index_t_now);
          for (MKL_LONG k=0; k<TRAJ.dimension; k++)
            {
              // printf("fc = %6.3e, fr = %6.3e, fR = %6.3e\n", force_spring(k), force_repulsion(k), force_random(k));
              TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k)) + sqrt(TRAJ.dt)*force_random[i](k);
            }
        }
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC
      double time_end_LV = dsecnd();

      // ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, index_t_next);
      ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_now);
      double time_end_AN = dsecnd();
      if(t%N_skip==0)
        {
          printf("STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_now), energy(1));
          printf("time consuming: MC, LV, AN, FILE = %8.6e, %8.6e, %8.6e, %8.6e\n", time_MC, time_LV, time_AN, time_file);
          double total_time = time_MC + time_LV + time_AN + time_file;
          printf("time fraction:  MC, LV, AN, FILE = %6.1f, %6.1f, %6.1f, %6.1f\n", time_MC*100/total_time, time_LV*100/total_time, time_AN*100/total_time, time_file*100/total_time);
          printf("MC step analysis: all pdf = %6.3e, basic_random = %6.3e, getting_hash = %6.3e, det_jump = %6.3e, new_end = %6.3e, action = %6.3e, update = %6.3e, MC_EQ_ident = %6.3e\n", dt_det_pdf, dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7);
          double total_dt = dt_1 + dt_2 + dt_3 + dt_4 + dt_5 + dt_6 + dt_7;
          printf("frac MC step analysis: all pdf = %6.1f, basic_random = %6.1f, getting_hash = %6.1f, det_jump = %6.1f, new_end = %6.1f, action = %6.1f, update = %6.1f, MC_EQ_ident = %6.1f\n", dt_det_pdf*100./total_dt, dt_1*100./total_dt, dt_2*100./total_dt, dt_3*100./total_dt, dt_4*100./total_dt, dt_5*100./total_dt, dt_6*100./total_dt, dt_7*100./total_dt);
          double total_dt_pdf = dt_pdf + dt_sort;
          printf("computing pdf: %6.3e (%3.1f), sorting pdf: %6.3e (%3.1f)\n", dt_pdf, 100.*dt_pdf/total_dt_pdf, dt_sort, dt_sort*100./total_dt_pdf);
          TRAJ.fprint_row(filename_trajectory.c_str(), index_t_now);
          energy.fprint(filename_energy.c_str());
          CONNECT.HASH.fprint(filename_HASH.c_str());
          CONNECT.CASE.fprint(filename_CASE.c_str());
          CONNECT.weight.fprint(filename_weight.c_str());
        }
      
      if(t%N_energy_frequency==0)
        {
          ANALYSIS::cal_detail_repulsion(TRAJ, POTs, filename_energy_info.c_str(), index_t_now);
        }

      double time_end_save = dsecnd();
      time_MC += time_end_MC - time_st_MC;
      time_LV += time_end_LV - time_end_MC;
      time_AN += time_end_AN - time_end_LV;
      time_file += time_end_save - time_end_AN;
    }

  double time_simulation = dsecnd() - time_st_simulation;
  printf("Total simulation time = %6.3e\n", time_simulation);
  if (given_condition("MC_LOG") == "TRUE")
    FILE_LOG.close();
  delete[] vec_boost_Nd_parallel;
  delete[] vec_boost_Np_parallel;
  for(MKL_LONG i=0; i<cond_longer; i++)
    {
      // delete[] vec_boost_Nd_Np_parallel[i];
      delete[] vec_boost_Nd_Np_parallel_connector[i];
      delete[] vec_boost_Nd_Np_parallel_repulsion[i];
    }
  // delete[] vec_boost_Nd_Np_parallel;
  delete[] vec_boost_Nd_Np_parallel_connector;
  delete[] vec_boost_Nd_Np_parallel_repulsion;
  delete[] force_spring;
  delete[] force_repulsion;
  delete[] force_random;
  delete[] dCDF_U;
  delete[] INDEX_dCDF_U;
  gsl_rng_free(r_boost); // for boosting
  return 0;
}



// -ipo option for interprocedure optimization for all files

/*
 * Local variables:
 * compile-command: "icpc -O2 -openmp -Wall -mkl -o Brownian_simulation lib_ed_cpp_BD/lib_traj.cpp lib_ed_cpp_BD/read_file_condition.cpp lib_ed_cpp_BD/lib_evolution.cpp lib_ed_cpp_BD/matrix_ed.cpp lib_ed_cpp_BD/matrix_long_ed.cpp lib_ed_cpp_BD/lib_potential.cpp lib_ed_cpp_BD/lib_connectivity.cpp lib_ed_cpp_BD/lib_association.cpp lib_ed_cpp_BD/lib_geometry.cpp lib_ed_cpp_BD/lib_random.cpp Brownian_simulation.cpp -lgsl -lm"
 * End:
 */



