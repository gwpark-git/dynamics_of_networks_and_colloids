
#include <iostream>
#include "lib_ed_cpp_BD/matrix_ed.h"
#include "lib_ed_cpp_BD/lib_traj.h"
#include "lib_ed_cpp_BD/lib_evolution.h"
#include "lib_ed_cpp_BD/lib_association.h"
#include "lib_ed_cpp_BD/lib_handle_association.h"
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

      POTENTIAL_SET POTs;
      if(given_condition("Method") == "NAPLE_REPULSION")
        {
          FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_set(POTs, given_condition);
          main_NAPLE_REPULSION(TRAJ, POTs, given_condition);
        }
      else if(given_condition("Method") == "NAPLE_ASSOCIATION")
        {
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
  printf("STARTING_MAIN_ROOT\n"); // ERR_TEST
  double time_st_simulation = dsecnd(); // ERR_TEST
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
  printf("THREAD_SETTING: %ld ... ", N_THREADS); //ERR_TEST
  omp_set_num_threads(N_THREADS);
  mkl_set_num_threads(N_THREADS);
  double time_MC = 0.;
  double time_LV = 0.;
  double time_AN = 0.;
  double time_file = 0.;
  printf("DONE\n"); // ERR_TEST
  // the following are the boosting the allocation and dislocation of the MATRIX classes
  printf("GENERATING BOOSTING VECTORS\n");
  MATRIX *vec_boost_Nd_parallel = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(TRAJ.dimension, 1, 0.);
    }
  MKL_LONG index_set[4] = {0};
  MKL_LONG &index_itself = index_set[2], &index_other_end_of_selected_chain = index_set[0], &index_new_end_of_selected_chain = index_set[1]; // this make identify the system
  MKL_LONG &index_hash_selected_chain = index_set[3];
  MKL_LONG (*ACTION_ARR[4])(TRAJECTORY&, MKL_LONG, POTENTIAL_SET&, ASSOCIATION&, MKL_LONG[], MATRIX&);
  ACTION_ARR[0] = ACTION::CANCEL; ACTION_ARR[1] = ACTION::ADD; ACTION_ARR[2] = ACTION::DEL; ACTION_ARR[3] = ACTION::MOV;
  printf("DONE\n"); // ERR_TEST
  printf("FORCE VECTOR GENERATING ... "); // ERR_TEST
  MATRIX *force_spring = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_repulsion = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_random = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      force_spring[i].initial(TRAJ.dimension, 1, 0.);
      force_repulsion[i].initial(TRAJ.dimension, 1, 0.);
      force_random[i].initial(TRAJ.dimension, 1, 0.);
    }
  printf("DONE\n"); // ERR_TEST
  if(given_condition("MC_LOG") == "TRUE")
    {
      FILE_LOG.open(filename_MC_LOG.c_str(), std::ios_base::out);
      FILE_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dCDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << endl;
    }
  printf("SET SIMULATION PARAMETERS ...");
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 4, 0.);
  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  // TRAJ.fprint_row(filename_trajectory.c_str(), 0);

  MKL_LONG N_basic = TRAJ.rows;

  double tolerance_association = atof(given_condition("tolerance_association").c_str());
  MKL_LONG N_max_steps = atol(given_condition("N_max_steps").c_str());
  printf("DONE\n");
  printf("GENERATING CDF and INDEX_CDF VECTORS ...");
  MATRIX *dCDF_U = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *INDEX_dCDF_U = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      dCDF_U[i].initial(TRAJ.Np, 1, 0.);
      INDEX_dCDF_U[i].initial(TRAJ.Np, 1, 0);
    }
  printf("DONE\n");
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  CONNECT.initial();
  for(MKL_LONG i=0; i<CONNECT.Np; i++)
    CONNECT.TOKEN[i] = 1;

  double dt_1 = 0., dt_2 = 0., dt_3 = 0., dt_4 = 0., dt_5 = 0., dt_6 = 0., dt_7 = 0.;
  double dt_det_pdf = 0.;
  double dt_pdf = 0., dt_sort = 0.;
  
  double time_MC_1 = 0., time_MC_2 = 0., time_MC_3 = 0., time_MC_4 = 0., time_MC_5 = 0., time_MC_6 = 0., time_MC_7 = 0., time_MC_8 = 0.;

  printf("GENERATING RANDOM VECTOR BOOST ... ");
  const gsl_rng_type *T_boost;
  gsl_rng *r_boost;
  gsl_rng_env_setup();
  T_boost = gsl_rng_default;
  r_boost = gsl_rng_alloc(T_boost);
  printf("DONE\n");
  printf("START SIMULATION\n");
  for(MKL_LONG t = 0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      TRAJ(index_t_next) = (++TRAJ.c_t) * TRAJ.dt;

      
      tmp_vec.set_value(0.);

      MKL_LONG IDENTIFIER_ASSOC = TRUE;
      double max_try_ASSOC = tolerance_association;
      double N_diff = 0.;
      double max_N_diff = 0.;
      MKL_LONG sum_over_MC_steps = 0;
      MKL_LONG count_M = 0;
      MKL_LONG pre_N_associations = 0;
      MKL_LONG N_associations = 0;

      MKL_LONG cnt = 1;
      if(given_condition("MC_renewal")=="TRUE")
        {
          CONNECT.initial();
          for(MKL_LONG i=0; i<CONNECT.Np; i++)
            CONNECT.TOKEN[i] = 1;
        }
      // MKL_LONG IDENT_CANCEL = 0, IDENT_ADD = 1, IDENT_DEL = 2, IDENT_MOV = 3;
      MKL_LONG cnt_arr[4] = {0};
      MKL_LONG &cnt_cancel = cnt_arr[ACTION::IDX_CANCEL], &cnt_add = cnt_arr[ACTION::IDX_ADD], &cnt_del = cnt_arr[ACTION::IDX_DEL], &cnt_mov = cnt_arr[ACTION::IDX_MOV];

      MKL_LONG N_index_boost_arr[4];
      N_index_boost_arr[ACTION::IDX_CANCEL] = 0; N_index_boost_arr[ACTION::IDX_ADD] = 2; N_index_boost_arr[ACTION::IDX_DEL] = 2; N_index_boost_arr[ACTION::IDX_MOV] = 3;
        
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

      if(given_condition("Step")!="EQUILIBRATION")
        {
          while(IDENTIFIER_ASSOC && ++cnt < N_max_steps)
            {
              // choice for visiting bead   
              // MKL_LONG total_chain_ends = CONNECT.N_TOTAL_CONNECTED_ENDS(); // this is no more importance
              // printf("chain ends attached to beads, %ld, per beads, %ld, (add, mov, del) = (%ld, %ld, %ld)\n", total_chain_ends, total_chain_ends/TRAJ.Np, cnt_add, cnt_mov, cnt_del);
              time_MC_1 = dsecnd();
              index_itself = RANDOM::return_LONG_INT_rand_boost(r_boost, TRAJ.Np);
              // choice for selected chain end
              double rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(r_boost);
              time_MC_2 = dsecnd();
              index_hash_selected_chain = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF); 
              index_other_end_of_selected_chain = CONNECT.HASH[index_itself](index_hash_selected_chain); 
              time_MC_3 = dsecnd();
              // choice for behaviour of selected chain end
              double rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(r_boost);
              // the PDF is already computed in the previous map
              time_MC_4 = dsecnd();
              
              MKL_LONG k = backsearch(dCDF_U[index_itself], rolling_dCDF_U);
              index_new_end_of_selected_chain = INDEX_dCDF_U[index_itself](k);
              time_MC_5 = dsecnd();

              MKL_LONG IDENTIFIER_ACTION = TRUTH_MAP::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, index_set);
              ACTION_ARR[IDENTIFIER_ACTION](TRAJ, index_t_now, POTs, CONNECT, index_set, tmp_vec);
              cnt_arr[IDENTIFIER_ACTION] ++;
              MKL_LONG N_index_boost = N_index_boost_arr[IDENTIFIER_ACTION];
              
              time_MC_6 = dsecnd();
              
              for(MKL_LONG i=0; i<N_index_boost; i++)
                {
                  CONNECTIVITY_update_Z_particle(CONNECT, index_set[i]);
                  CONNECTIVITY_update_dPDF_particle(CONNECT, index_set[i]);
                  CONNECTIVITY_update_dCDF_particle(CONNECT, index_set[i]);
                }
          
              time_MC_7 = dsecnd();

              pre_N_associations = N_associations;
              N_associations = cnt_add - cnt_del;

              if(N_associations != pre_N_associations)
                {
                  count_M ++;
                  if (count_M%N_steps_block == 0)
                    {
                      {
                        N_diff = fabs((double)(1./(count_M*(count_M-1)))*((count_M-1)*N_associations - sum_over_MC_steps));
                        max_N_diff = max_N_diff > N_diff ? max_N_diff : N_diff;
                        if(N_diff/max_N_diff < tolerance_association)
                          {
                            IDENTIFIER_ASSOC = FALSE;
                          }
                      }
                    }
                  sum_over_MC_steps += N_associations;
                }

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
                  FILE_LOG << cnt << '\t' << index_itself << '\t' << rolling_dCDF<< '\t'  << index_hash_selected_chain<< '\t'  << index_other_end_of_selected_chain<< '\t'  << rolling_dCDF_U<< '\t'  << k<< '\t'  << index_new_end_of_selected_chain<< '\t'  << CONNECT.TOKEN[index_itself]<< '\t'<<CONNECT.N_CONNECTED_ENDS(index_itself) << '\t' << CONNECT.weight[index_itself](0) <<'\t' <<  total_bonds << '\t'  << cnt_add<< '\t'  << cnt_mov<< '\t'  << cnt_del<< '\t'  << cnt_cancel << endl;
                }

              //   }
            } // while
        } // equilibration condition check
      double time_end_MC = dsecnd();

      // printf("STARTING MD\n");
// #pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next) firstprivate(force_spring, force_repulsion, force_random) // firstprivate called copy-constructor while private called default constructor
// #pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next, vec_boost_Nd_parallel, force_spring, force_repulsion, force_random) schedule(dynamic) num_threads(N_THREADS)
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
              TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k)) + sqrt(TRAJ.dt)*force_random[i](k);
            }
        }
      // printf("STARTING ANALYSIS\n");
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC
      double time_end_LV = dsecnd();

      ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_now, vec_boost_Nd_parallel[0]);
      double time_end_AN = dsecnd();
      // printf("WRITING DATA\n");
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
          for(MKL_LONG ip=0; ip<TRAJ.Np; ip++)
            {
              CONNECT.HASH[ip].fprint_LONG_transpose(filename_HASH.c_str());
              CONNECT.CASE[ip].fprint_transpose(filename_CASE.c_str());
              CONNECT.weight[ip].fprint_LONG_transpose(filename_weight.c_str());
            }
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
  mkl_free(vec_boost_Nd_parallel);
  mkl_free(force_spring);
  mkl_free(force_repulsion);
  mkl_free(force_random);
  mkl_free(dCDF_U);
  mkl_free(INDEX_dCDF_U);
  gsl_rng_free(r_boost); // for boosting
  return 0;
}



/*
 * Local variables:
 * compile-command: "icpc -ipo -O2 -Wall -mkl -o Brownian_simulation lib_ed_cpp_BD/lib_traj.cpp lib_ed_cpp_BD/read_file_condition.cpp lib_ed_cpp_BD/lib_evolution.cpp lib_ed_cpp_BD/matrix_ed.cpp lib_ed_cpp_BD/lib_potential.cpp lib_ed_cpp_BD/lib_connectivity.cpp lib_ed_cpp_BD/lib_association.cpp lib_ed_cpp_BD/lib_handle_association.cpp lib_ed_cpp_BD/lib_geometry.cpp lib_ed_cpp_BD/lib_random.cpp Brownian_simulation.cpp -lgsl -lm"
 * End:
 */

/*
  -ipo: interprocedure optimization
  -O2: optimization level 2
  -Wall: rigirous compile option (no warning point is allowed)
  -mkl: Math Kernel Librar (Intel)
  -lgsl: GSL library
  -lm: Math library
 */

