
#include <iostream>
#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/association.h"
#include "../lib/handle_association.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"
#include <string>
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
/*

  Think: should I include internal time measurment after update of code?
  The boosting mechanism is of importance for further implementation.
  However, the measurment itself is somehow time consumming.
  Do I bear it until end of my simulation? Think about it.

  struct TIME
  {
  double ST_simulation;
  double ST_MC;
  double ST_PDF, END_PDF;
  double PDF_ST_RDIST, PDF_ST_PDF, PDF_ST_SORT;

  double 
  
  
  }
*/

MKL_LONG main_NAPLE_ASSOCIATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, COND& given_condition);
MKL_LONG main_NAPLE_ASSOCIATION_TRACKING_CHAINS(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, COND& given_condition);

MKL_LONG main_EQUILIBRATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, COND& given_condition);

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
          // Note that the only needed information for Euler is the previous time step. Therefore, we need previous and now which is in 2.
          N_basic = 2;

        }
      TRAJECTORY TRAJ(given_condition, N_basic);
      POTENTIAL_SET POTs;
      if(given_condition("Method") == "NAPLE_ASSOCIATION")
        {
          if(given_condition("Step") == "EQUILIBRATION")
            {
              FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition);
              main_EQUILIBRATION(TRAJ, POTs, given_condition);
            }
          else
            {
              ASSOCIATION CONNECT(TRAJ, given_condition);
              FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition);
              if (given_condition("tracking_individual_chain") == "TRUE")
                {
		  printf("PRE_INIT_CHAIN\n");
                  CHAIN_HANDLE CHAIN(given_condition);
		  printf("PRE_ALLOC_CHAIN\n");
                  CHAIN.allocate_existing_bridges(CONNECT);
		  printf("AFT_CHAIN\n");
                  main_NAPLE_ASSOCIATION_TRACKING_CHAINS(TRAJ, POTs, CONNECT, CHAIN, given_condition);
                }
              else
                {
                  main_NAPLE_ASSOCIATION(TRAJ, POTs, CONNECT, given_condition);
                }
            }
        }
      else
        {
          printf("Method except NAPLE_ASSOCIATION is not clearly defined. The given Method input is %s", given_condition("Method").c_str());
        }
    }
  return 0;
}



MKL_LONG main_NAPLE_ASSOCIATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, COND& given_condition)
{
  printf("STARTING_MAIN_ROOT\n"); // ERR_TEST
  double time_st_simulation = dsecnd(); // ERR_TEST
  // MKL_LONG N_max_steps = atol(given_condition("N_max_steps").c_str());
  MKL_LONG N_steps_block = atol(given_condition("N_steps_block").c_str());
  // MKL_LONG N_max_blocks = N_max_steps/N_steps_block;
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_base") + ".traj").c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_base") + ".ener").c_str();
  string filename_HASH = (given_condition("output_path") + '/' + given_condition("filename_base") + ".hash").c_str();
  string filename_weight = (given_condition("output_path") + '/' + given_condition("filename_base") + ".weight").c_str();
  string filename_MC_LOG = (given_condition("output_path") + '/' + given_condition("filename_base") + ".MC_LOG").c_str();
  // string filename_info = (given_condition("output_path") + '/' + given_condition("file_base") + ".info").c_str();
  
  ofstream FILE_LOG;
  
  MKL_LONG N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
  MKL_LONG N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());
  printf("THREAD_SETTING: %ld ... ", N_THREADS_BD); //ERR_TEST
  mkl_set_num_threads(N_THREADS_BD);
  double time_MC = 0.;
  double time_LV = 0.;
  double time_AN = 0.;
  double time_file = 0.;
  printf("DONE\n");
  printf("GENERATING BOOSTING VECTORS\n");
  LOCK LOCKER(TRAJ.Np);
  MKL_LONG *tmp_index_vec = (MKL_LONG*) mkl_malloc(TRAJ.dimension*sizeof(MKL_LONG), BIT);
  MATRIX *vec_boost_Nd_parallel = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); 
  // MATRIX **R_minimum_vec_boost = (MATRIX**) mkl_malloc(TRAJ.Np*sizeof(MATRIX*), BIT); //RDIST
  // MATRIX *R_minimum_distance_boost = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); // RDIST
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(TRAJ.dimension, 1, 0.);
      // R_minimum_vec_boost[i] = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); // RDIST
      // for(MKL_LONG j=0; j<TRAJ.Np; j++)
      //   R_minimum_vec_boost[i][j].initial(3, 1, 0.); // RDIST
      // R_minimum_distance_boost[i].initial(TRAJ.Np, 1, 0.); // RDIST
    }
  RDIST R_boost(given_condition);
  
  MKL_LONG cnt_arr[5] = {0};
  MKL_LONG &cnt_cancel = cnt_arr[INDEX_MC::CANCEL], &cnt_add = cnt_arr[INDEX_MC::ADD], &cnt_del = cnt_arr[INDEX_MC::OPP_DEL], &cnt_mov = cnt_arr[INDEX_MC::MOV], &cnt_lock = cnt_arr[INDEX_MC::LOCK];
  cnt_add = CONNECT.N_ASSOCIATION;
  // printf("N_ASSOCIATION = %d\n\n", cnt_add);
  // cnt_add = 100;
  
  printf("DONE\n"); 
  printf("FORCE VECTOR GENERATING ... "); 
  MATRIX *force_spring = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_repulsion = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_random = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      force_spring[i].initial(TRAJ.dimension, 1, 0.);
      force_repulsion[i].initial(TRAJ.dimension, 1, 0.);
      force_random[i].initial(TRAJ.dimension, 1, 0.);
    }
  printf("DONE\n"); 
  if(given_condition("MC_LOG") == "TRUE")
    {
      FILE_LOG.open(filename_MC_LOG.c_str(), std::ios_base::out);
      FILE_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dCDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << '\t' << "16_cnt_lock" << endl;
    }
  printf("SET SIMULATION PARAMETERS ...");
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 6, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  // MATRIX energy(1, 12, 0.);

  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  MKL_LONG N_basic = TRAJ.rows;

  double tolerance_association = atof(given_condition("tolerance_association").c_str());
  printf("DONE\n");
  printf("GENERATING CDF and INDEX_CDF VECTORS ...");
  MATRIX *dCDF_U = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MKL_LONG *dCDF_TOKEN = (MKL_LONG*) mkl_malloc(TRAJ.Np*sizeof(MKL_LONG), BIT);
  MATRIX *INDEX_dCDF_U = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      dCDF_U[i].initial(TRAJ.Np, 1, 0.);
      INDEX_dCDF_U[i].initial(TRAJ.Np, 1, 0);
    }
  printf("DONE\n");
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  // The following condition duplicate and may violate the inheritance scheme from the previous association
  // CONNECT.initial();
  // for(MKL_LONG i=0; i<CONNECT.Np; i++)
  //   CONNECT.TOKEN[i] = 1;

  double dt_1 = 0., dt_2 = 0., dt_3 = 0., dt_4 = 0., dt_5 = 0., dt_6 = 0., dt_7 = 0.;
  double dt_det_pdf = 0.;
  double dt_rdist = 0., dt_pdf = 0., dt_sort = 0.;
  double time_MC_1 = 0., time_MC_2 = 0., time_MC_3 = 0., time_MC_4 = 0., time_MC_5 = 0., time_MC_6 = 0., time_MC_7 = 0., time_MC_8 = 0.;

  printf("GENERATING RANDOM VECTOR BOOST ... ");
  /*
    "gsl_rng *r_boost;" is the basic notation. i.e., it is already given by pointer.
    Hence, it is of importance that the dynamic allocation set with its pointer type.
  */
  const gsl_rng_type *T_boost;

  gsl_rng **r_boost_arr = (gsl_rng**)mkl_malloc(N_THREADS_BD*sizeof(gsl_rng*), BIT);
  gsl_rng_env_setup();
  T_boost = gsl_rng_default;
  for(MKL_LONG i=0; i<N_THREADS_BD; i++)
    {
      r_boost_arr[i] = gsl_rng_alloc(T_boost);
      // gsl_rng_set(r_boost_arr[i], random());
      gsl_rng_set(r_boost_arr[i], i);
    }

  gsl_rng **r_boost_arr_SS = (gsl_rng**)mkl_malloc(N_THREADS_SS*sizeof(gsl_rng*), BIT);  
  for(MKL_LONG i=0; i<N_THREADS_SS; i++)
    {
      r_boost_arr_SS[i] = gsl_rng_alloc(T_boost);
      gsl_rng_set(r_boost_arr_SS[i], i+N_THREADS_BD);
    }
  
  INDEX_MC *IDX_ARR = (INDEX_MC*) mkl_malloc(N_THREADS_BD*sizeof(INDEX_MC), BIT);
  for(MKL_LONG i=0; i<N_THREADS_SS; i++)
    {
      IDX_ARR[i].initial();
      IDX_ARR[i].set_initial_variables();
    }
  
  printf("DONE\n");
  printf("START SIMULATION\n");
  MKL_LONG N_associations = 0;

  // MKL_LONG IDENTIFIER_ASSOC = TRUE; // the identification is disabled
  double max_try_ASSOC = tolerance_association;
  double N_diff = 0.;
  MKL_LONG N_tot_associable_chain = TRAJ.Np*atoi(given_condition("N_chains_per_particle").c_str());
  // MKL_LONG count_M = cnt_add, pre_count_M = cnt_add;
  // MKL_LONG count_M = cnt_add;
  
  for(MKL_LONG t = 0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
      
      tmp_vec.set_value(0.);

      MKL_LONG cnt = 1;

      double time_st_rdist = dsecnd();
      R_boost.allocate_cells_from_positions(TRAJ, index_t_now, tmp_index_vec);
      
      // MKL_LONG tmp_index = 397;
      // MKL_LONG cell_tmp_index = R_boost.cell_index[tmp_index];
      // for(MKL_LONG j=0; j<R_boost.TOKEN[tmp_index]; j++)
      // 	{
	  
      // }
      // for(MKL_LONG i=0; i<R_boost.N_cells; i++)
      // 	{
      // 	  for(MKL_LONG j=0; j<R_boost.TOKEN[i]; j++)
      // 	    {
      // 	      printf("C[%ld, %ld] = %ld\n", i, j, R_boost.CELL[i][j]);
      // 	    }
      // 	}
      // #pragma omp parallel default(none) shared(given_condition, TRAJ, index_t_now, t, R_boost, POTs, CONNECT, dCDF_U, INDEX_dCDF_U, vec_boost_Nd_parallel, N_steps_block, cnt_arr, count_M, N_THREADS_BD, dt_rdist, dt_pdf, dt_sort, time_st_rdist) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      //       {
      // #pragma omp for
#pragma omp parallel for default(none) shared(TRAJ, index_t_now, R_boost, N_THREADS_BD) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      /*
        Originally, this parallel regime is designed to use N_cells and TOKEN.
        Because the chunk size is not easily specified, it is used to parallel with index of particles instead of cell-based.
        On this regards, the cell_index array is used which will return the cell index for the subjected particle.
      */
      for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
        {
          MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
          for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
            {
              MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
              for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
                {
                  MKL_LONG index_target = R_boost(cell_index_neighbor, p);
                  // double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
                  // printf("(%4.1e, %4.1e, %4.1e), ", R_boost.Rvec[index_particle][index_target](0), R_boost.Rvec[index_particle][index_target](1), R_boost.Rvec[index_particle][index_target](2));
                  double distance = GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.BEYOND_BOX[cell_index_particle][k]);
                  // printf("(%4.1e, %4.1e, %4.1e)\n ", R_boost.Rvec[index_particle][index_target](0), R_boost.Rvec[index_particle][index_target](1), R_boost.Rvec[index_particle][index_target](2));
                  // printf("(%4.1e, %4.1e, %4.1e, d2 = %4.1e, %4.1e\n", GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target]), GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.BEYOND_BOX[cell_index_particle][k]));
                  R_boost.Rsca[index_particle](index_target) = distance;
                } // p
            } // k
        } // index_particle
      dt_rdist += dsecnd() - time_st_rdist;
      double time_st_MC = dsecnd();
      // if(given_condition("Step")!="EQUILIBRATION" && t%N_steps_block == 0) // including initial time t=0
      if(t%N_steps_block == 0) // the equilibration functionality is disabled at this moment. (it will be seperated for future works)
        {
          // this is rearranged in order to use the previously updated information when MC_renewal is not turned on.
          if(given_condition("MC_renewal")=="TRUE")
            {
              CONNECT.set_initial_condition();
              // reset counting array
              for(MKL_LONG i=0; i<4; i++)
                {
                  cnt_arr[i] = 0;
                }
              // count_M = 0;
            }
          else // MC_renewal check
            {
              
              /*
                The following will check only for the existing brdige.
                Since the bridges will only occurred withing neighbouring cell-list, this bridge chain will not violate the condition for cell-list.
              */	      
#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, R_boost, vec_boost_Nd_parallel) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
              for(MKL_LONG i=0; i<TRAJ.Np; i++)
                {
                  for(MKL_LONG j=0; j<CONNECT.TOKEN[i]; j++)
                    {
                      // CONNECT.HASH[i](j) gave us the index for target
                      // which means we have to compute distance between i and k where k is given by CONNECT.HASH[i](j).
                      // CONNECT.update_CASE_particle_hash_target(POTs, i, j, R_minimum_distance_boost[i](CONNECT.HASH[i](j))); // RDIST
                      // if(R_boost.Rsca[i](CONNECT.HASH[i](j)) > 2.0)
                      // 	{
                      // 	  printf("d=%4.1f\n", R_boost.Rsca[i](CONNECT.HASH[i](j)));
                      // 	}
                      CONNECT.update_CASE_particle_hash_target(POTs, i, j, R_boost.Rsca[i](CONNECT.HASH[i](j)));
                    }
                  CONNECT.update_Z_particle(i);
                  CONNECT.update_dPDF_particle(i);
                  CONNECT.update_dCDF_particle(i);
                  // printf("Z[%ld] = %3.2lf, TOKEN[%ld] = %ld\n", i, CONNECT.Z[i], i, (MKL_LONG)CONNECT.TOKEN[i]);
                } 
            }  // else for MC_renewal check
          double time_st_pdf = dsecnd();
          // #pragma omp for
#pragma omp parallel for default(none) shared(TRAJ, POTs, dCDF_U, INDEX_dCDF_U, dCDF_TOKEN, R_boost, N_THREADS_BD, dt_pdf, dt_sort) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
          for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
            {
              MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
              MKL_LONG count_CDF_TOKEN = 0;
              dCDF_TOKEN[index_particle] = 0;
              INDEX_dCDF_U[index_particle].set_value(-1);
              dCDF_U[index_particle].set_value(0);
              for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
                {
                  MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
                  // if(index_particle==1)
                  //   {
                  //     printf("CELL_INFO: %ld, %ld\n", k, cell_index_neighbor);
                  //   }
                  for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
                    {
                      MKL_LONG index_target = R_boost(cell_index_neighbor, p);
                      double distance = R_boost.Rsca[index_particle](index_target);
                      INDEX_dCDF_U[index_particle](count_CDF_TOKEN) = index_target;
                      dCDF_U[index_particle](count_CDF_TOKEN) = POTs.PDF_connector(distance, POTs.force_variables);
                      if(dCDF_U[index_particle](count_CDF_TOKEN) > 0.0)
                        {
                          dCDF_TOKEN[index_particle] ++;
                        }
                      // if(index_particle==1)
                      // 	{
                      // 	  printf("ts=%ld: Rsca[%ld](%ld) = %4.1e (dCDF=%4.1e), k=%ld, CI=%ld, p=%ld, NCI=%ld, TOKEN_NCI=%ld\n", TRAJ.c_t, index_particle, index_target, R_boost.Rsca[index_particle](index_target), dCDF_U[index_particle](count_CDF_TOKEN), k, cell_index_particle, p, cell_index_neighbor, R_boost.TOKEN[cell_index_neighbor]);
                      // 	}
		      
                      count_CDF_TOKEN ++;
		      
                    } // p
                } // k
              // dCDF_TOKEN[index_particle] = count_CDF_TOKEN;
              dCDF_U[index_particle].sort2(INDEX_dCDF_U[index_particle]);
              for(MKL_LONG k=TRAJ.Np-dCDF_TOKEN[index_particle] + 1; k<TRAJ.Np; k++)
                {
                  dCDF_U[index_particle](k) += dCDF_U[index_particle](k-1);
                }
              for(MKL_LONG k=TRAJ.Np-dCDF_TOKEN[index_particle]; k<TRAJ.Np; k++)
                {
                  dCDF_U[index_particle](k) /= dCDF_U[index_particle](TRAJ.Np - 1);
                }
            } // index_particle, local parallel
          // printf("tmp_check\n");
          // for(MKL_LONG i=0; i<TRAJ.Np; i++)
          //   {
          //     printf("dCDF_U[%ld](399) = %4.1e\n", i, dCDF_U[i](399));
          //   }
#pragma omp parallel for default(none) shared(given_condition, FILE_LOG, TRAJ, POTs, CONNECT, LOCKER, IDX_ARR, index_t_now, vec_boost_Nd_parallel, INDEX_dCDF_U, dCDF_U, dCDF_TOKEN, R_boost, dt_rdist, dt_pdf, dt_sort, cnt_arr, cnt_add, cnt_del, cnt_mov, cnt_cancel, cnt_lock, N_steps_block, r_boost_arr_SS, cnt, N_THREADS_SS, N_associations, N_tot_associable_chain) private(time_MC_1, time_MC_2, time_MC_3, time_MC_4, time_MC_5, time_MC_6, time_MC_7, time_MC_8) num_threads(N_THREADS_SS) if(N_THREADS_SS > 1) reduction(+:dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7)
          // #pragma omp critical(TOPOLOGY_UPDATE)
          // {
          // for(MKL_LONG tp = 0; tp<N_steps_block; tp++)
          for(MKL_LONG tp=0; tp<N_tot_associable_chain; tp++)
            {
              /*
                'it' have the identity number for current thread. Then, the reference variables IDX and r_boost will be used in order to usability and readability. In this case, IDX, r_boost is just reference of existing one, but IDX and r_boost itself is local reference variables which will varied thread to thread
              */
              MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
              MKL_LONG &index_itself = IDX_ARR[it].beads[CONNECT.flag_itself];
              MKL_LONG &index_attached_bead = IDX_ARR[it].beads[CONNECT.flag_other];
              MKL_LONG &index_new_attached_bead = IDX_ARR[it].beads[CONNECT.flag_new];
              MKL_LONG &index_hash_attached_bead = IDX_ARR[it].beads[CONNECT.flag_hash_other];
                  
              time_MC_1 = dsecnd();
              index_itself = RANDOM::return_LONG_INT_rand_boost(r_boost_arr_SS[it], TRAJ.Np);
              // choice for selected chain end
              double rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(r_boost_arr_SS[it]);
              time_MC_2 = dsecnd();
              index_hash_attached_bead = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF); 
              index_attached_bead = CONNECT.HASH[index_itself](index_hash_attached_bead); 
              time_MC_3 = dsecnd();
              // choice for behaviour of selected chain end
              double rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(r_boost_arr_SS[it]);
              // the PDF is already computed in the previous map
              time_MC_4 = dsecnd();
              
              // MKL_LONG k = SEARCHING::backtrace(dCDF_U[index_itself], rolling_dCDF_U);
              MKL_LONG k = SEARCHING::backtrace_cell_list(dCDF_U[index_itself], dCDF_TOKEN[index_itself], rolling_dCDF_U, index_itself, R_boost);
              index_new_attached_bead = INDEX_dCDF_U[index_itself](k);
              // if(k!=TRAJ.Np - 1)
              // 	{
              // 	  printf("k, index=%ld, %ld\n", k, index_new_attached_bead);
              // 	}
              time_MC_5 = dsecnd();
              MKL_LONG IDENTIFIER_ACTION = TRUE; // it can be 1 (IDX.ADD) but just true value
              MKL_LONG IDENTIFIER_LOCKING = FALSE;
              if(N_THREADS_SS > 1)
                {
#pragma omp critical(LOCKING)  // LOCKING is the name for this critical blocks
                  {
                  /*
                    On the omp critical region, the block will work only one thread.
                    If the other thread reaching this reason while there is one thread already working on this block, then the reached thread will wait until finishing the job of the other thread.
                    This benefits to identify the working beads index on this case, since the 
                  */
                  // CHECKING
                  for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++)
                    {
                      if(LOCKER(IDX_ARR[it].beads[I_BEADS]))
                        {
                          IDENTIFIER_ACTION = IDX_ARR[it].CANCEL;
                          IDENTIFIER_LOCKING = TRUE;
                          break;
                        }
                    }
                  // this is LOCKING procedure
                  if(!IDENTIFIER_LOCKING)
                    {
                      cnt++;  // preventing LOCKING affect to the IDENTIFICATION of stochastic balance
                      for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++)
                        {
                          LOCKER(IDX_ARR[it].beads[I_BEADS]) = TRUE;
                        }
                    }
                  else
                    {
                      cnt_lock ++;
                    }
                } // critical(LOCKING)
            } // check the parallel
              time_MC_6 = dsecnd();
              // Note that the critical region only applicable with single thread while the others will be used in parallel regime.
              // In addition, the gap for passing the critical region will tune further gaps, then the computation speed for passing critical region will not be real critical issue.
              double time_MC_pre_ACTION = 0., time_MC_end_ACTION = 0., time_MC_end_UPDATE=0.;
              if(!IDENTIFIER_LOCKING) 
                {
                  // This block only compute when the thread is NOT LOCKED

                  time_MC_pre_ACTION = dsecnd();

                  // double distance_exist_bridge = R_minimum_distance_boost[index_itself](index_attached_bead); // RDIST
                  double distance_exist_bridge = R_boost.Rsca[index_itself](index_attached_bead);
                  double tpa = POTs.transition(distance_exist_bridge, POTs.f_connector(distance_exist_bridge, POTs.force_variables), POTs.force_variables);
                  if (tpa == 1.0)
                    {
                      IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX_ARR[it]);
                    }
                  else
                    {
                      double rolling_transition = RANDOM::return_double_rand_SUP1_boost(r_boost_arr_SS[it]);
                      if (rolling_transition < tpa)
                        {
                          IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX_ARR[it]);
                          // if(index_itself != index_new_attached_bead)
                          //   printf("IDA=%ld, rolled=%4.1e, tpa=%4.1e, (%ld,%ld)\n", IDENTIFIER_ACTION, rolling_transition, tpa, index_itself, index_new_attached_bead);
                        }
                      else
                        IDENTIFIER_ACTION = IDX_ARR[it].CANCEL;
                    }

                  // ACTION::ACT(TRAJ, index_t_now, POTs, CONNECT, IDX_ARR[it], R_minimum_distance_boost, IDENTIFIER_ACTION); // RDIST
                  ACTION::ACT(TRAJ, index_t_now, POTs, CONNECT, IDX_ARR[it], R_boost.Rsca, IDENTIFIER_ACTION);
                  time_MC_end_ACTION = dsecnd();
                  ACTION::UPDATE_INFORMATION(CONNECT, IDX_ARR[it], cnt_arr, IDENTIFIER_ACTION);
                  time_MC_end_UPDATE = dsecnd();

                  // UNLOCKING
                  // The critical directive is no more necessarly since only one thread visited each beads
                  for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++)
                    {
                      LOCKER(IDX_ARR[it].beads[I_BEADS]) = FALSE;
                    }

                  {
                    // update for reduction part
                    // note that the reduction is individually working
                    // when the job in parallel region is finished, it will be summed over all the different threads
                    dt_1 += time_MC_2 - time_MC_1; // basic_random
                    dt_2 += time_MC_3 - time_MC_2; // getting_hash
                    dt_3 += time_MC_4 - time_MC_3; // det_jump
                    dt_4 += time_MC_5 - time_MC_4; // new_end
                    dt_5 += time_MC_6 - time_MC_5; // LOCKING
                    dt_6 += time_MC_end_ACTION - time_MC_pre_ACTION; // ACTION
                    dt_7 += time_MC_end_UPDATE - time_MC_end_ACTION; // UPDATE
                  }

                  /*
                    critical(COUNTING) blocks:
                    This is counting the action information that will be used for the future.
                    Notice that the writing MC_LOG file is inside of this COUNTING critical directive, since all the information should be the same for writing (temporal)
                  */
#pragma omp critical(COUNTING) 
                  {			
			
                    cnt_arr[IDENTIFIER_ACTION]++;
                    N_associations = cnt_add - cnt_del;
		    
                    // count_M += N_associations;
                      
                    if (given_condition("MC_LOG") == "TRUE")
                      {
                        MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
			  
                        // MKL_LONG count_N_associagtions = cnt_add - cnt_del;
                        {
                          FILE_LOG << cnt << '\t' << index_itself << '\t' << setprecision(7) << rolling_dCDF<< '\t'  << index_attached_bead << '\t'  << index_new_attached_bead<< '\t'  << setprecision(7) << rolling_dCDF_U<< '\t'  << k<< '\t'  << index_new_attached_bead << '\t'  << CONNECT.TOKEN[index_itself]<< '\t'<< CONNECT.N_CONNECTED_ENDS(index_itself) << '\t' << CONNECT.weight[index_itself](0) <<'\t' <<  total_bonds << '\t'  << cnt_add<< '\t'  << cnt_mov<< '\t'  << cnt_del<< '\t'  << cnt_cancel << '\t' << cnt_lock << endl;
                        }
                        // FILE_LOG << boost::format("%10d\t%4d\t")
                      } // MC_LOG
			
                  } // critical(COUNTING)
                } // if(!IDENTIFIER_LOCKING)
        } // for loop (ASSOCIATION)
    } // region for topological update
      double time_end_MC = dsecnd();

#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next, R_boost, vec_boost_Nd_parallel, force_spring, force_repulsion, force_random, r_boost_arr, N_THREADS_BD, given_condition) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      for (MKL_LONG i=0; i<TRAJ.Np; i++)
        {
          MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
          
          force_spring[i].set_value(0);
          force_repulsion[i].set_value(0);
          force_random[i].set_value(0);

          // if(given_condition("Step")!="EQUILIBRATION") // EQUILIBRIATION is disabled
          //   {
          // INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(TRAJ, POTs, CONNECT, force_spring[i], index_t_now, i, R_minimum_vec_boost, R_minimum_distance_boost); // RDIST
          INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(TRAJ, POTs, CONNECT, force_spring[i], index_t_now, i, R_boost.Rvec, R_boost.Rsca);
          // }
          // INTEGRATOR::EULER::cal_repulsion_force_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, R_minimum_vec_boost, R_minimum_distance_boost); // RDIST
          // INTEGRATOR::EULER::cal_repulsion_force_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, R_boost.Rvec, R_boost.Rsca);
          INTEGRATOR::EULER::cal_repulsion_force_R_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, R_boost);
          INTEGRATOR::EULER::cal_random_force_boost(TRAJ, POTs, force_random[i], index_t_now, r_boost_arr[it]); 
          for (MKL_LONG k=0; k<TRAJ.dimension; k++)
            {
              TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k)) + sqrt(TRAJ.dt)*force_random[i](k);
            }
        }
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC
      double time_end_LV = dsecnd();
      double time_end_AN = time_end_LV;
      if(t%N_skip==0)
        {
          time_end_LV = dsecnd();
          energy(0) = TRAJ(index_t_now);
          energy(4) = (double)N_associations;
          energy(5) = dsecnd() - time_st_simulation;
          ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_now, vec_boost_Nd_parallel[0]);
          time_end_AN = dsecnd();
          double total_dt = dt_1 + dt_2 + dt_3 + dt_4 + dt_5 + dt_6 + dt_7;
          double total_dt_pdf = dt_rdist + dt_pdf + dt_sort;
          double total_time = time_MC + time_LV + time_AN + time_file + total_dt_pdf;
          double dt_pdf_all = dt_pdf + dt_sort;
          
          printf("##### STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_now), energy(1));
          printf("time consuming: MC, LV, AN, FILE, DIST = %8.6e, %8.6e, %8.6e, %8.6e, %8.6e\n", time_MC, time_LV, time_AN, time_file, total_dt_pdf);
          printf("time fraction:  MC, LV, AN, FILE, DIST = %6.1f, %6.1f, %6.1f, %6.1f, %6.1f\n", time_MC*100/total_time, time_LV*100/total_time, time_AN*100/total_time, time_file*100/total_time, total_dt_pdf*100/total_time);
          printf("MC step analysis: all pdf = %6.3e, basic_random = %6.3e, getting_hash = %6.3e, det_jump = %6.3e, new_end = %6.3e, LOCKING = %6.3e, action = %6.3e, update = %6.3e\n", dt_pdf_all, dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7);
          printf("frac MC step analysis: all pdf = %6.1f, basic_random = %6.1f, getting_hash = %6.1f, det_jump = %6.1f, new_end = %6.1f, LOCKING = %6.3f, action = %6.1f, update = %6.1f\n", dt_pdf_all*100./total_dt, dt_1*100./total_dt, dt_2*100./total_dt, dt_3*100./total_dt, dt_4*100./total_dt, dt_5*100./total_dt, dt_6*100./total_dt, dt_7*100./total_dt);
          printf("computing rdist: %6.3e (%3.1f), computing pdf: %6.3e (%3.1f), sorting pdf: %6.3e (%3.1f)\n", dt_rdist, 100.*dt_rdist/total_dt_pdf, dt_pdf, 100.*dt_pdf/total_dt_pdf, dt_sort, dt_sort*100./total_dt_pdf);
          printf("LAST IDENTIFIER: cnt = %ld, N_diff = %6.3e, N_tot_asso = %ld, ratio = %6.3e, NAS = %ld, fraction=%4.3f, total time=%4.3e ####\n\n", cnt, N_diff, N_tot_associable_chain, N_diff/N_tot_associable_chain, N_associations, N_associations/(double)N_tot_associable_chain, energy(5));
          TRAJ.fprint_row(filename_trajectory.c_str(), index_t_now);
          energy.fprint_row(filename_energy.c_str(), 0);
          for(MKL_LONG ip=0; ip<TRAJ.Np; ip++)
            {
              CONNECT.HASH[ip].fprint_LONG_skip_transpose_LIMROWS(filename_HASH.c_str(), 1, CONNECT.TOKEN[ip]);
              CONNECT.weight[ip].fprint_LONG_skip_transpose_LIMROWS(filename_weight.c_str(), 1, CONNECT.TOKEN[ip]);
            }
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
//   for(MKL_LONG i=0; i<TRAJ.Np; i++)
//     mkl_free(R_minimum_vec_boost[i]); // RDIST
// mkl_free(R_minimum_vec_boost); // RDIST
// mkl_free(R_minimum_distance_boost); // RDIST
mkl_free(tmp_index_vec);
mkl_free(dCDF_U);
mkl_free(INDEX_dCDF_U);
mkl_free(dCDF_TOKEN);  
for(MKL_LONG i=0; i<N_THREADS_BD; i++)
  gsl_rng_free(r_boost_arr[i]); // for boosting

mkl_free(r_boost_arr);
mkl_free(IDX_ARR);
return 0;
}


MKL_LONG main_EQUILIBRATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, COND& given_condition)
{
  printf("STARTING_MAIN_ROOT\n"); // ERR_TEST
  double time_st_simulation = dsecnd(); // ERR_TEST
  // MKL_LONG N_max_steps = atol(given_condition("N_max_steps").c_str());
  MKL_LONG N_steps_block = atol(given_condition("N_steps_block").c_str());
  // MKL_LONG N_max_blocks = N_max_steps/N_steps_block;
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_base") + ".traj").c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_base") + ".ener").c_str();

  ofstream FILE_LOG;
  
  MKL_LONG N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
  printf("THREAD_SETTING: %ld ... ", N_THREADS_BD); //ERR_TEST
  mkl_set_num_threads(N_THREADS_BD);
  double time_MC = 0.;
  double time_LV = 0.;
  double time_AN = 0.;
  double time_file = 0.;
  printf("DONE\n");
  printf("GENERATING BOOSTING VECTORS\n");

  MKL_LONG *tmp_index_vec = (MKL_LONG*) mkl_malloc(TRAJ.dimension*sizeof(MKL_LONG), BIT);
  MATRIX *vec_boost_Nd_parallel = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); 

  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(TRAJ.dimension, 1, 0.);
    }
  RDIST R_boost(given_condition);
  
  printf("DONE\n"); 
  printf("FORCE VECTOR GENERATING ... "); 


  MATRIX *force_repulsion = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_random = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      force_repulsion[i].initial(TRAJ.dimension, 1, 0.);
      force_random[i].initial(TRAJ.dimension, 1, 0.);
    }
  printf("DONE\n"); 
  printf("SET SIMULATION PARAMETERS ...");
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 6, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  // MATRIX energy(1, 12, 0.);

  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  MKL_LONG N_basic = TRAJ.rows;

  double tolerance_association = atof(given_condition("tolerance_association").c_str());
  printf("DONE\n");
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  // The following condition duplicate and may violate the inheritance scheme from the previous association
  // CONNECT.initial();
  // for(MKL_LONG i=0; i<CONNECT.Np; i++)
  //   CONNECT.TOKEN[i] = 1;

  double dt_1 = 0., dt_2 = 0., dt_3 = 0., dt_4 = 0., dt_5 = 0., dt_6 = 0., dt_7 = 0.;
  double dt_det_pdf = 0.;
  double dt_rdist = 0., dt_pdf = 0., dt_sort = 0.;
  double time_MC_1 = 0., time_MC_2 = 0., time_MC_3 = 0., time_MC_4 = 0., time_MC_5 = 0., time_MC_6 = 0., time_MC_7 = 0., time_MC_8 = 0.;

  printf("GENERATING RANDOM VECTOR BOOST ... ");
  /*
    "gsl_rng *r_boost;" is the basic notation. i.e., it is already given by pointer.
    Hence, it is of importance that the dynamic allocation set with its pointer type.
  */
  const gsl_rng_type *T_boost;

  gsl_rng **r_boost_arr = (gsl_rng**)mkl_malloc(N_THREADS_BD*sizeof(gsl_rng*), BIT);
  gsl_rng_env_setup();
  T_boost = gsl_rng_default;
  for(MKL_LONG i=0; i<N_THREADS_BD; i++)
    {
      r_boost_arr[i] = gsl_rng_alloc(T_boost);
      // gsl_rng_set(r_boost_arr[i], random());
      gsl_rng_set(r_boost_arr[i], i);
    }

  
  printf("DONE\n");
  printf("START SIMULATION\n");
  MKL_LONG N_associations = 0;

  // MKL_LONG IDENTIFIER_ASSOC = TRUE; // the identification is disabled
  double max_try_ASSOC = tolerance_association;
  double N_diff = 0.;
  MKL_LONG N_tot_associable_chain = TRAJ.Np*atoi(given_condition("N_chains_per_particle").c_str());
  // MKL_LONG count_M = cnt_add, pre_count_M = cnt_add;
  // MKL_LONG count_M = cnt_add;
  
  for(MKL_LONG t = 0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
      
      tmp_vec.set_value(0.);

      MKL_LONG cnt = 1;

      double time_st_rdist = dsecnd();
      R_boost.allocate_cells_from_positions(TRAJ, index_t_now, tmp_index_vec);
      
#pragma omp parallel for default(none) shared(TRAJ, index_t_now, R_boost, N_THREADS_BD) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      /*
        Originally, this parallel regime is designed to use N_cells and TOKEN.
        Because the chunk size is not easily specified, it is used to parallel with index of particles instead of cell-based.
        On this regards, the cell_index array is used which will return the cell index for the subjected particle.
      */
      for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
        {
          MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
          for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
            {
              MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
              for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
                {
                  MKL_LONG index_target = R_boost(cell_index_neighbor, p);
                  // double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
                  // printf("(%4.1e, %4.1e, %4.1e), ", R_boost.Rvec[index_particle][index_target](0), R_boost.Rvec[index_particle][index_target](1), R_boost.Rvec[index_particle][index_target](2));
                  double distance = GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.BEYOND_BOX[cell_index_particle][k]);
                  // printf("(%4.1e, %4.1e, %4.1e)\n ", R_boost.Rvec[index_particle][index_target](0), R_boost.Rvec[index_particle][index_target](1), R_boost.Rvec[index_particle][index_target](2));
                  // printf("(%4.1e, %4.1e, %4.1e, d2 = %4.1e, %4.1e\n", GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target]), GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.BEYOND_BOX[cell_index_particle][k]));
                  R_boost.Rsca[index_particle](index_target) = distance;
                } // p
            } // k
        } // index_particle
      dt_rdist += dsecnd() - time_st_rdist;
      double time_st_MC = dsecnd();
      // if(given_condition("Step")!="EQUILIBRATION" && t%N_steps_block == 0) // including initial time t=0
      double time_end_MC = dsecnd();

#pragma omp parallel for default(none) shared(TRAJ, POTs, index_t_now, index_t_next, R_boost, vec_boost_Nd_parallel, force_repulsion, force_random, r_boost_arr, N_THREADS_BD, given_condition) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      for (MKL_LONG i=0; i<TRAJ.Np; i++)
        {
          MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
          
          force_repulsion[i].set_value(0);
          force_random[i].set_value(0);

          INTEGRATOR::EULER::cal_repulsion_force_R_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, R_boost);
          INTEGRATOR::EULER::cal_random_force_boost(TRAJ, POTs, force_random[i], index_t_now, r_boost_arr[it]); 
          for (MKL_LONG k=0; k<TRAJ.dimension; k++)
            {
              TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*(force_repulsion[i](k)) + sqrt(TRAJ.dt)*force_random[i](k);
            }
        }
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC
      double time_end_LV = dsecnd();
      double time_end_AN = time_end_LV;
      if(t%N_skip==0)
        {
          time_end_LV = dsecnd();
          energy(0) = TRAJ(index_t_now);
          energy(4) = (double)N_associations;
          energy(5) = dsecnd() - time_st_simulation;
          // ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_now, vec_boost_Nd_parallel[0]);
          time_end_AN = dsecnd();
          double total_dt = dt_1 + dt_2 + dt_3 + dt_4 + dt_5 + dt_6 + dt_7;
          double total_dt_pdf = dt_rdist + dt_pdf + dt_sort;
          double total_time = time_MC + time_LV + time_AN + time_file + total_dt_pdf;
          double dt_pdf_all = dt_pdf + dt_sort;
          
          printf("##### STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_now), energy(1));
          printf("time consuming: MC, LV, AN, FILE, DIST = %8.6e, %8.6e, %8.6e, %8.6e, %8.6e\n", time_MC, time_LV, time_AN, time_file, total_dt_pdf);
          printf("time fraction:  MC, LV, AN, FILE, DIST = %6.1f, %6.1f, %6.1f, %6.1f, %6.1f\n", time_MC*100/total_time, time_LV*100/total_time, time_AN*100/total_time, time_file*100/total_time, total_dt_pdf*100/total_time);
          printf("MC step analysis: all pdf = %6.3e, basic_random = %6.3e, getting_hash = %6.3e, det_jump = %6.3e, new_end = %6.3e, LOCKING = %6.3e, action = %6.3e, update = %6.3e\n", dt_pdf_all, dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7);
          printf("frac MC step analysis: all pdf = %6.1f, basic_random = %6.1f, getting_hash = %6.1f, det_jump = %6.1f, new_end = %6.1f, LOCKING = %6.3f, action = %6.1f, update = %6.1f\n", dt_pdf_all*100./total_dt, dt_1*100./total_dt, dt_2*100./total_dt, dt_3*100./total_dt, dt_4*100./total_dt, dt_5*100./total_dt, dt_6*100./total_dt, dt_7*100./total_dt);
          printf("computing rdist: %6.3e (%3.1f), computing pdf: %6.3e (%3.1f), sorting pdf: %6.3e (%3.1f)\n", dt_rdist, 100.*dt_rdist/total_dt_pdf, dt_pdf, 100.*dt_pdf/total_dt_pdf, dt_sort, dt_sort*100./total_dt_pdf);
          printf("LAST IDENTIFIER: cnt = %ld, N_diff = %6.3e, N_tot_asso = %ld, ratio = %6.3e, NAS = %ld, fraction=%4.3f, total time=%4.3e ####\n\n", cnt, N_diff, N_tot_associable_chain, N_diff/N_tot_associable_chain, N_associations, N_associations/(double)N_tot_associable_chain, energy(5));
          TRAJ.fprint_row(filename_trajectory.c_str(), index_t_now);
          energy.fprint_row(filename_energy.c_str(), 0);
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
  mkl_free(force_repulsion);
  mkl_free(force_random);
  //   for(MKL_LONG i=0; i<TRAJ.Np; i++)
  //     mkl_free(R_minimum_vec_boost[i]); // RDIST
  // mkl_free(R_minimum_vec_boost); // RDIST
  // mkl_free(R_minimum_distance_boost); // RDIST
  mkl_free(tmp_index_vec);
  for(MKL_LONG i=0; i<N_THREADS_BD; i++)
    gsl_rng_free(r_boost_arr[i]); // for boosting

  mkl_free(r_boost_arr);
  return 0;
}

MKL_LONG main_NAPLE_ASSOCIATION_TRACKING_CHAINS(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, COND& given_condition)
{
  printf("STARTING_MAIN_ROOT\n"); // ERR_TEST
  double time_st_simulation = dsecnd(); // ERR_TEST
  // MKL_LONG N_max_steps = atol(given_condition("N_max_steps").c_str());
  MKL_LONG N_steps_block = atol(given_condition("N_steps_block").c_str());
  // MKL_LONG N_max_blocks = N_max_steps/N_steps_block;
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_base") + ".traj").c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_base") + ".ener").c_str();
  string filename_HASH = (given_condition("output_path") + '/' + given_condition("filename_base") + ".hash").c_str();
  string filename_weight = (given_condition("output_path") + '/' + given_condition("filename_base") + ".weight").c_str();
  string filename_MC_LOG = (given_condition("output_path") + '/' + given_condition("filename_base") + ".MC_LOG").c_str();
  // string filename_info = (given_condition("output_path") + '/' + given_condition("file_base") + ".info").c_str();
  
  ofstream FILE_LOG;
  
  MKL_LONG N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
  MKL_LONG N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());
  printf("THREAD_SETTING: %ld ... ", N_THREADS_BD); //ERR_TEST
  mkl_set_num_threads(N_THREADS_BD);
  double time_MC = 0.;
  double time_LV = 0.;
  double time_AN = 0.;
  double time_file = 0.;
  printf("DONE\n");
  printf("GENERATING BOOSTING VECTORS\n");
  LOCK LOCKER(TRAJ.Np);
  MKL_LONG *tmp_index_vec = (MKL_LONG*) mkl_malloc(TRAJ.dimension*sizeof(MKL_LONG), BIT);
  MATRIX *vec_boost_Nd_parallel = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); 
  // MATRIX **R_minimum_vec_boost = (MATRIX**) mkl_malloc(TRAJ.Np*sizeof(MATRIX*), BIT); //RDIST
  // MATRIX *R_minimum_distance_boost = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); // RDIST
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(TRAJ.dimension, 1, 0.);
      // R_minimum_vec_boost[i] = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); // RDIST
      // for(MKL_LONG j=0; j<TRAJ.Np; j++)
      //   R_minimum_vec_boost[i][j].initial(3, 1, 0.); // RDIST
      // R_minimum_distance_boost[i].initial(TRAJ.Np, 1, 0.); // RDIST
    }
  RDIST R_boost(given_condition);
  
  MKL_LONG cnt_arr[5] = {0};
  MKL_LONG &cnt_cancel = cnt_arr[INDEX_MC::CANCEL], &cnt_add = cnt_arr[INDEX_MC::ADD], &cnt_del = cnt_arr[INDEX_MC::OPP_DEL], &cnt_mov = cnt_arr[INDEX_MC::MOV], &cnt_lock = cnt_arr[INDEX_MC::LOCK];
  cnt_add = CONNECT.N_ASSOCIATION;
  // printf("N_ASSOCIATION = %d\n\n", cnt_add);
  // cnt_add = 100;
  
  printf("DONE\n"); 
  printf("FORCE VECTOR GENERATING ... "); 
  MATRIX *force_spring = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_repulsion = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_random = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      force_spring[i].initial(TRAJ.dimension, 1, 0.);
      force_repulsion[i].initial(TRAJ.dimension, 1, 0.);
      force_random[i].initial(TRAJ.dimension, 1, 0.);
    }
  printf("DONE\n"); 
  if(given_condition("MC_LOG") == "TRUE")
    {
      FILE_LOG.open(filename_MC_LOG.c_str(), std::ios_base::out);
      FILE_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dCDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << '\t' << "16_cnt_lock" << endl;
    }
  printf("SET SIMULATION PARAMETERS ...");
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 6, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  // MATRIX energy(1, 12, 0.);

  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  MKL_LONG N_basic = TRAJ.rows;

  double tolerance_association = atof(given_condition("tolerance_association").c_str());
  printf("DONE\n");
  printf("GENERATING CDF and INDEX_CDF VECTORS ...");
  MATRIX *dCDF_U = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MKL_LONG *dCDF_TOKEN = (MKL_LONG*) mkl_malloc(TRAJ.Np*sizeof(MKL_LONG), BIT);
  MATRIX *INDEX_dCDF_U = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      dCDF_U[i].initial(TRAJ.Np, 1, 0.);
      INDEX_dCDF_U[i].initial(TRAJ.Np, 1, 0);
    }
  printf("DONE\n");
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
  // The following condition duplicate and may violate the inheritance scheme from the previous association
  // CONNECT.initial();
  // for(MKL_LONG i=0; i<CONNECT.Np; i++)
  //   CONNECT.TOKEN[i] = 1;

  double dt_1 = 0., dt_2 = 0., dt_3 = 0., dt_4 = 0., dt_5 = 0., dt_6 = 0., dt_7 = 0.;
  double dt_det_pdf = 0.;
  double dt_rdist = 0., dt_pdf = 0., dt_sort = 0.;
  double time_MC_1 = 0., time_MC_2 = 0., time_MC_3 = 0., time_MC_4 = 0., time_MC_5 = 0., time_MC_6 = 0., time_MC_7 = 0., time_MC_8 = 0.;

  printf("GENERATING RANDOM VECTOR BOOST ... ");
  /*
    "gsl_rng *r_boost;" is the basic notation. i.e., it is already given by pointer.
    Hence, it is of importance that the dynamic allocation set with its pointer type.
  */
  const gsl_rng_type *T_boost;

  gsl_rng **r_boost_arr = (gsl_rng**)mkl_malloc(N_THREADS_BD*sizeof(gsl_rng*), BIT);
  gsl_rng_env_setup();
  T_boost = gsl_rng_default;
  for(MKL_LONG i=0; i<N_THREADS_BD; i++)
    {
      r_boost_arr[i] = gsl_rng_alloc(T_boost);
      // gsl_rng_set(r_boost_arr[i], random());
      gsl_rng_set(r_boost_arr[i], i);
    }

  gsl_rng **r_boost_arr_SS = (gsl_rng**)mkl_malloc(N_THREADS_SS*sizeof(gsl_rng*), BIT);  
  for(MKL_LONG i=0; i<N_THREADS_SS; i++)
    {
      r_boost_arr_SS[i] = gsl_rng_alloc(T_boost);
      gsl_rng_set(r_boost_arr_SS[i], i+N_THREADS_BD);
    }
  
  INDEX_MC *IDX_ARR = (INDEX_MC*) mkl_malloc(N_THREADS_BD*sizeof(INDEX_MC), BIT);
  for(MKL_LONG i=0; i<N_THREADS_SS; i++)
    {
      IDX_ARR[i].initial();
      IDX_ARR[i].set_initial_variables();
    }
  
  printf("DONE\n");
  printf("START SIMULATION\n");
  MKL_LONG N_associations = 0;

  // MKL_LONG IDENTIFIER_ASSOC = TRUE; // the identification is disabled
  double max_try_ASSOC = tolerance_association;
  double N_diff = 0.;
  MKL_LONG N_tot_associable_chain = TRAJ.Np*atoi(given_condition("N_chains_per_particle").c_str());
  // MKL_LONG count_M = cnt_add, pre_count_M = cnt_add;
  // MKL_LONG count_M = cnt_add;
  
  for(MKL_LONG t = 0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
      
      tmp_vec.set_value(0.);

      MKL_LONG cnt = 1;

      double time_st_rdist = dsecnd();
      R_boost.allocate_cells_from_positions(TRAJ, index_t_now, tmp_index_vec);
      
      // MKL_LONG tmp_index = 397;
      // MKL_LONG cell_tmp_index = R_boost.cell_index[tmp_index];
      // for(MKL_LONG j=0; j<R_boost.TOKEN[tmp_index]; j++)
      // 	{
	  
      // }
      // for(MKL_LONG i=0; i<R_boost.N_cells; i++)
      // 	{
      // 	  for(MKL_LONG j=0; j<R_boost.TOKEN[i]; j++)
      // 	    {
      // 	      printf("C[%ld, %ld] = %ld\n", i, j, R_boost.CELL[i][j]);
      // 	    }
      // 	}
      // #pragma omp parallel default(none) shared(given_condition, TRAJ, index_t_now, t, R_boost, POTs, CONNECT, dCDF_U, INDEX_dCDF_U, vec_boost_Nd_parallel, N_steps_block, cnt_arr, count_M, N_THREADS_BD, dt_rdist, dt_pdf, dt_sort, time_st_rdist) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      //       {
      // #pragma omp for
#pragma omp parallel for default(none) shared(TRAJ, index_t_now, R_boost, N_THREADS_BD) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      /*
        Originally, this parallel regime is designed to use N_cells and TOKEN.
        Because the chunk size is not easily specified, it is used to parallel with index of particles instead of cell-based.
        On this regards, the cell_index array is used which will return the cell index for the subjected particle.
      */
      for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
        {
          MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
          for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
            {
              MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
              for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
                {
                  MKL_LONG index_target = R_boost(cell_index_neighbor, p);
                  // double distance = GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target]);
                  // printf("(%4.1e, %4.1e, %4.1e), ", R_boost.Rvec[index_particle][index_target](0), R_boost.Rvec[index_particle][index_target](1), R_boost.Rvec[index_particle][index_target](2));
                  double distance = GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.BEYOND_BOX[cell_index_particle][k]);
                  // printf("(%4.1e, %4.1e, %4.1e)\n ", R_boost.Rvec[index_particle][index_target](0), R_boost.Rvec[index_particle][index_target](1), R_boost.Rvec[index_particle][index_target](2));
                  // printf("(%4.1e, %4.1e, %4.1e, d2 = %4.1e, %4.1e\n", GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target]), GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t_now, index_particle, index_target, R_boost.Rvec[index_particle][index_target], R_boost.BEYOND_BOX[cell_index_particle][k]));
                  R_boost.Rsca[index_particle](index_target) = distance;
                } // p
            } // k
        } // index_particle
      dt_rdist += dsecnd() - time_st_rdist;
      double time_st_MC = dsecnd();
      // if(given_condition("Step")!="EQUILIBRATION" && t%N_steps_block == 0) // including initial time t=0
      if(t%N_steps_block == 0) // the equilibration functionality is disabled at this moment. (it will be seperated for future works)
        {
          // this is rearranged in order to use the previously updated information when MC_renewal is not turned on.
          if(given_condition("MC_renewal")=="TRUE")
            {
              CONNECT.set_initial_condition();
              // reset counting array
              for(MKL_LONG i=0; i<4; i++)
                {
                  cnt_arr[i] = 0;
                }
              // count_M = 0;
            }
          else // MC_renewal check
            {
              
              /*
                The following will check only for the existing brdige.
                Since the bridges will only occurred withing neighbouring cell-list, this bridge chain will not violate the condition for cell-list.
              */	      
#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, R_boost, vec_boost_Nd_parallel) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
              for(MKL_LONG i=0; i<TRAJ.Np; i++)
                {
                  for(MKL_LONG j=0; j<CONNECT.TOKEN[i]; j++)
                    {
                      // CONNECT.HASH[i](j) gave us the index for target
                      // which means we have to compute distance between i and k where k is given by CONNECT.HASH[i](j).
                      // CONNECT.update_CASE_particle_hash_target(POTs, i, j, R_minimum_distance_boost[i](CONNECT.HASH[i](j))); // RDIST
                      // if(R_boost.Rsca[i](CONNECT.HASH[i](j)) > 2.0)
                      // 	{
                      // 	  printf("d=%4.1f\n", R_boost.Rsca[i](CONNECT.HASH[i](j)));
                      // 	}
                      CONNECT.update_CASE_particle_hash_target(POTs, i, j, R_boost.Rsca[i](CONNECT.HASH[i](j)));
                    }
                  CONNECT.update_Z_particle(i);
                  CONNECT.update_dPDF_particle(i);
                  CONNECT.update_dCDF_particle(i);
                  // printf("Z[%ld] = %3.2lf, TOKEN[%ld] = %ld\n", i, CONNECT.Z[i], i, (MKL_LONG)CONNECT.TOKEN[i]);
                } 
            }  // else for MC_renewal check
          double time_st_pdf = dsecnd();
          // #pragma omp for
#pragma omp parallel for default(none) shared(TRAJ, POTs, dCDF_U, INDEX_dCDF_U, dCDF_TOKEN, R_boost, N_THREADS_BD, dt_pdf, dt_sort) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
          for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
            {
              MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
              MKL_LONG count_CDF_TOKEN = 0;
              dCDF_TOKEN[index_particle] = 0;
              INDEX_dCDF_U[index_particle].set_value(-1);
              dCDF_U[index_particle].set_value(0);
              for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
                {
                  MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
                  // if(index_particle==1)
                  //   {
                  //     printf("CELL_INFO: %ld, %ld\n", k, cell_index_neighbor);
                  //   }
                  for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
                    {
                      MKL_LONG index_target = R_boost(cell_index_neighbor, p);
                      double distance = R_boost.Rsca[index_particle](index_target);
                      INDEX_dCDF_U[index_particle](count_CDF_TOKEN) = index_target;
                      dCDF_U[index_particle](count_CDF_TOKEN) = POTs.PDF_connector(distance, POTs.force_variables);
                      if(dCDF_U[index_particle](count_CDF_TOKEN) > 0.0)
                        {
                          dCDF_TOKEN[index_particle] ++;
                        }
                      // if(index_particle==1)
                      // 	{
                      // 	  printf("ts=%ld: Rsca[%ld](%ld) = %4.1e (dCDF=%4.1e), k=%ld, CI=%ld, p=%ld, NCI=%ld, TOKEN_NCI=%ld\n", TRAJ.c_t, index_particle, index_target, R_boost.Rsca[index_particle](index_target), dCDF_U[index_particle](count_CDF_TOKEN), k, cell_index_particle, p, cell_index_neighbor, R_boost.TOKEN[cell_index_neighbor]);
                      // 	}
		      
                      count_CDF_TOKEN ++;
		      
                    } // p
                } // k
              // dCDF_TOKEN[index_particle] = count_CDF_TOKEN;
              dCDF_U[index_particle].sort2(INDEX_dCDF_U[index_particle]);
              for(MKL_LONG k=TRAJ.Np-dCDF_TOKEN[index_particle] + 1; k<TRAJ.Np; k++)
                {
                  dCDF_U[index_particle](k) += dCDF_U[index_particle](k-1);
                }
              for(MKL_LONG k=TRAJ.Np-dCDF_TOKEN[index_particle]; k<TRAJ.Np; k++)
                {
                  dCDF_U[index_particle](k) /= dCDF_U[index_particle](TRAJ.Np - 1);
                }
            } // index_particle, local parallel
          // printf("tmp_check\n");
          // for(MKL_LONG i=0; i<TRAJ.Np; i++)
          //   {
          //     printf("dCDF_U[%ld](399) = %4.1e\n", i, dCDF_U[i](399));
          //   }
#pragma omp parallel for default(none) shared(given_condition, FILE_LOG, TRAJ, POTs, CONNECT, CHAIN, LOCKER, IDX_ARR, index_t_now, vec_boost_Nd_parallel, INDEX_dCDF_U, dCDF_U, dCDF_TOKEN, R_boost, dt_rdist, dt_pdf, dt_sort, cnt_arr, cnt_add, cnt_del, cnt_mov, cnt_cancel, cnt_lock, N_steps_block, r_boost_arr_SS, cnt, N_THREADS_SS, N_associations, N_tot_associable_chain) private(time_MC_1, time_MC_2, time_MC_3, time_MC_4, time_MC_5, time_MC_6, time_MC_7, time_MC_8) num_threads(N_THREADS_SS) if(N_THREADS_SS > 1) reduction(+:dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7)
          // #pragma omp critical(TOPOLOGY_UPDATE)
          // {
          // for(MKL_LONG tp = 0; tp<N_steps_block; tp++)
          for(MKL_LONG tp=0; tp<N_tot_associable_chain; tp++)
            {
              /*
                'it' have the identity number for current thread. Then, the reference variables IDX and r_boost will be used in order to usability and readability. In this case, IDX, r_boost is just reference of existing one, but IDX and r_boost itself is local reference variables which will varied thread to thread
              */
              MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
              MKL_LONG &index_itself = IDX_ARR[it].beads[CONNECT.flag_itself];
              MKL_LONG &index_attached_bead = IDX_ARR[it].beads[CONNECT.flag_other];
              MKL_LONG &index_new_attached_bead = IDX_ARR[it].beads[CONNECT.flag_new];
              MKL_LONG &index_hash_attached_bead = IDX_ARR[it].beads[CONNECT.flag_hash_other];
                  
              time_MC_1 = dsecnd();
              index_itself = RANDOM::return_LONG_INT_rand_boost(r_boost_arr_SS[it], TRAJ.Np);
              // choice for selected chain end
              double rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(r_boost_arr_SS[it]);
              time_MC_2 = dsecnd();
              index_hash_attached_bead = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF); 
              index_attached_bead = CONNECT.HASH[index_itself](index_hash_attached_bead); 
              time_MC_3 = dsecnd();
              // choice for behaviour of selected chain end
              double rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(r_boost_arr_SS[it]);
              // the PDF is already computed in the previous map
              time_MC_4 = dsecnd();
              
              // MKL_LONG k = SEARCHING::backtrace(dCDF_U[index_itself], rolling_dCDF_U);
              MKL_LONG k = SEARCHING::backtrace_cell_list(dCDF_U[index_itself], dCDF_TOKEN[index_itself], rolling_dCDF_U, index_itself, R_boost);
              index_new_attached_bead = INDEX_dCDF_U[index_itself](k);
              // if(k!=TRAJ.Np - 1)
              // 	{
              // 	  printf("k, index=%ld, %ld\n", k, index_new_attached_bead);
              // 	}
              time_MC_5 = dsecnd();
              MKL_LONG IDENTIFIER_ACTION = TRUE; // it can be 1 (IDX.ADD) but just true value
              MKL_LONG IDENTIFIER_LOCKING = FALSE;
              if(N_THREADS_SS > 1)
                {
#pragma omp critical(LOCKING)  // LOCKING is the name for this critical blocks
                  {
                    /*
                      On the omp critical region, the block will work only one thread.
                      If the other thread reaching this reason while there is one thread already working on this block, then the reached thread will wait until finishing the job of the other thread.
                      This benefits to identify the working beads index on this case, since the 
                    */
                    // CHECKING
                    for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++)
                      {
                        if(LOCKER(IDX_ARR[it].beads[I_BEADS]))
                          {
                            IDENTIFIER_ACTION = IDX_ARR[it].CANCEL;
                            IDENTIFIER_LOCKING = TRUE;
                            break;
                          }
                      }
                    // this is LOCKING procedure
                    if(!IDENTIFIER_LOCKING)
                      {
                        cnt++;  // preventing LOCKING affect to the IDENTIFICATION of stochastic balance
                        for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++)
                          {
                            LOCKER(IDX_ARR[it].beads[I_BEADS]) = TRUE;
                          }
                      }
                    else
                      {
                        cnt_lock ++;
                      }
                  } // critical(LOCKING)
                } // check the parallel
              time_MC_6 = dsecnd();
              // Note that the critical region only applicable with single thread while the others will be used in parallel regime.
              // In addition, the gap for passing the critical region will tune further gaps, then the computation speed for passing critical region will not be real critical issue.
              double time_MC_pre_ACTION = 0., time_MC_end_ACTION = 0., time_MC_end_UPDATE=0.;
              if(!IDENTIFIER_LOCKING) 
                {
                  // This block only compute when the thread is NOT LOCKED

                  time_MC_pre_ACTION = dsecnd();

                  // double distance_exist_bridge = R_minimum_distance_boost[index_itself](index_attached_bead); // RDIST
                  double distance_exist_bridge = R_boost.Rsca[index_itself](index_attached_bead);
                  double tpa = POTs.transition(distance_exist_bridge, POTs.f_connector(distance_exist_bridge, POTs.force_variables), POTs.force_variables);
                  if (tpa == 1.0)
                    {
                      IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX_ARR[it]);
                    }
                  else
                    {
                      double rolling_transition = RANDOM::return_double_rand_SUP1_boost(r_boost_arr_SS[it]);
                      if (rolling_transition < tpa)
                        {
                          IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX_ARR[it]);
                          // if(index_itself != index_new_attached_bead)
                          //   printf("IDA=%ld, rolled=%4.1e, tpa=%4.1e, (%ld,%ld)\n", IDENTIFIER_ACTION, rolling_transition, tpa, index_itself, index_new_attached_bead);
                        }
                      else
                        IDENTIFIER_ACTION = IDX_ARR[it].CANCEL;
                    }

                  // ACTION::ACT(TRAJ, index_t_now, POTs, CONNECT, IDX_ARR[it], R_minimum_distance_boost, IDENTIFIER_ACTION); // RDIST
                  ACTION::ACT(TRAJ, index_t_now, POTs, CONNECT, IDX_ARR[it], R_boost.Rsca, IDENTIFIER_ACTION);
                  // CHAIN.TRACKING_ACTION(CONNECT, IDENTIFIER_ACTION, IDX_ARR[it]); // it will track individual chain information
                  
                  time_MC_end_ACTION = dsecnd();
                  ACTION::UPDATE_INFORMATION(CONNECT, IDX_ARR[it], cnt_arr, IDENTIFIER_ACTION);
                  time_MC_end_UPDATE = dsecnd();

                  // UNLOCKING
                  // The critical directive is no more necessarly since only one thread visited each beads
                  for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++)
                    {
                      LOCKER(IDX_ARR[it].beads[I_BEADS]) = FALSE;
                    }

                  {
                    // update for reduction part
                    // note that the reduction is individually working
                    // when the job in parallel region is finished, it will be summed over all the different threads
                    dt_1 += time_MC_2 - time_MC_1; // basic_random
                    dt_2 += time_MC_3 - time_MC_2; // getting_hash
                    dt_3 += time_MC_4 - time_MC_3; // det_jump
                    dt_4 += time_MC_5 - time_MC_4; // new_end
                    dt_5 += time_MC_6 - time_MC_5; // LOCKING
                    dt_6 += time_MC_end_ACTION - time_MC_pre_ACTION; // ACTION
                    dt_7 += time_MC_end_UPDATE - time_MC_end_ACTION; // UPDATE
                  }

                  /*
                    critical(COUNTING) blocks:
                    This is counting the action information that will be used for the future.
                    Notice that the writing MC_LOG file is inside of this COUNTING critical directive, since all the information should be the same for writing (temporal)
                  */
#pragma omp critical(COUNTING) 
                  {			
			
                    cnt_arr[IDENTIFIER_ACTION]++;
                    N_associations = cnt_add - cnt_del;
		    
                    // count_M += N_associations;
                      
                    if (given_condition("MC_LOG") == "TRUE")
                      {
                        MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
			  
                        // MKL_LONG count_N_associagtions = cnt_add - cnt_del;
                        {
                          FILE_LOG << cnt << '\t' << index_itself << '\t' << setprecision(7) << rolling_dCDF<< '\t'  << index_attached_bead << '\t'  << index_new_attached_bead<< '\t'  << setprecision(7) << rolling_dCDF_U<< '\t'  << k<< '\t'  << index_new_attached_bead << '\t'  << CONNECT.TOKEN[index_itself]<< '\t'<< CONNECT.N_CONNECTED_ENDS(index_itself) << '\t' << CONNECT.weight[index_itself](0) <<'\t' <<  total_bonds << '\t'  << cnt_add<< '\t'  << cnt_mov<< '\t'  << cnt_del<< '\t'  << cnt_cancel << '\t' << cnt_lock << endl;
                        }
                        // FILE_LOG << boost::format("%10d\t%4d\t")
                      } // MC_LOG
			
                  } // critical(COUNTING)
                } // if(!IDENTIFIER_LOCKING)
            } // for loop (ASSOCIATION)
        } // region for topological update
      double time_end_MC = dsecnd();

#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next, R_boost, vec_boost_Nd_parallel, force_spring, force_repulsion, force_random, r_boost_arr, N_THREADS_BD, given_condition) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      for (MKL_LONG i=0; i<TRAJ.Np; i++)
        {
          MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
          
          force_spring[i].set_value(0);
          force_repulsion[i].set_value(0);
          force_random[i].set_value(0);

          // if(given_condition("Step")!="EQUILIBRATION") // EQUILIBRIATION is disabled
          //   {
          // INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(TRAJ, POTs, CONNECT, force_spring[i], index_t_now, i, R_minimum_vec_boost, R_minimum_distance_boost); // RDIST
          INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(TRAJ, POTs, CONNECT, force_spring[i], index_t_now, i, R_boost.Rvec, R_boost.Rsca);
          // }
          // INTEGRATOR::EULER::cal_repulsion_force_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, R_minimum_vec_boost, R_minimum_distance_boost); // RDIST
          // INTEGRATOR::EULER::cal_repulsion_force_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, R_boost.Rvec, R_boost.Rsca);
          INTEGRATOR::EULER::cal_repulsion_force_R_boost(TRAJ, POTs, force_repulsion[i], index_t_now, i, R_boost);
          INTEGRATOR::EULER::cal_random_force_boost(TRAJ, POTs, force_random[i], index_t_now, r_boost_arr[it]); 
          for (MKL_LONG k=0; k<TRAJ.dimension; k++)
            {
              TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k)) + sqrt(TRAJ.dt)*force_random[i](k);
            }
        }
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC
      double time_end_LV = dsecnd();
      double time_end_AN = time_end_LV;
      if(t%N_skip==0)
        {
          time_end_LV = dsecnd();
          energy(0) = TRAJ(index_t_now);
          energy(4) = (double)N_associations;
          energy(5) = dsecnd() - time_st_simulation;
          ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_now, vec_boost_Nd_parallel[0]);
          time_end_AN = dsecnd();
          double total_dt = dt_1 + dt_2 + dt_3 + dt_4 + dt_5 + dt_6 + dt_7;
          double total_dt_pdf = dt_rdist + dt_pdf + dt_sort;
          double total_time = time_MC + time_LV + time_AN + time_file + total_dt_pdf;
          double dt_pdf_all = dt_pdf + dt_sort;
          
          printf("##### STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_now), energy(1));
          printf("time consuming: MC, LV, AN, FILE, DIST = %8.6e, %8.6e, %8.6e, %8.6e, %8.6e\n", time_MC, time_LV, time_AN, time_file, total_dt_pdf);
          printf("time fraction:  MC, LV, AN, FILE, DIST = %6.1f, %6.1f, %6.1f, %6.1f, %6.1f\n", time_MC*100/total_time, time_LV*100/total_time, time_AN*100/total_time, time_file*100/total_time, total_dt_pdf*100/total_time);
          printf("MC step analysis: all pdf = %6.3e, basic_random = %6.3e, getting_hash = %6.3e, det_jump = %6.3e, new_end = %6.3e, LOCKING = %6.3e, action = %6.3e, update = %6.3e\n", dt_pdf_all, dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7);
          printf("frac MC step analysis: all pdf = %6.1f, basic_random = %6.1f, getting_hash = %6.1f, det_jump = %6.1f, new_end = %6.1f, LOCKING = %6.3f, action = %6.1f, update = %6.1f\n", dt_pdf_all*100./total_dt, dt_1*100./total_dt, dt_2*100./total_dt, dt_3*100./total_dt, dt_4*100./total_dt, dt_5*100./total_dt, dt_6*100./total_dt, dt_7*100./total_dt);
          printf("computing rdist: %6.3e (%3.1f), computing pdf: %6.3e (%3.1f), sorting pdf: %6.3e (%3.1f)\n", dt_rdist, 100.*dt_rdist/total_dt_pdf, dt_pdf, 100.*dt_pdf/total_dt_pdf, dt_sort, dt_sort*100./total_dt_pdf);
          printf("LAST IDENTIFIER: cnt = %ld, N_diff = %6.3e, N_tot_asso = %ld, ratio = %6.3e, NAS = %ld, fraction=%4.3f, total time=%4.3e ####\n\n", cnt, N_diff, N_tot_associable_chain, N_diff/N_tot_associable_chain, N_associations, N_associations/(double)N_tot_associable_chain, energy(5));
          TRAJ.fprint_row(filename_trajectory.c_str(), index_t_now);
          energy.fprint_row(filename_energy.c_str(), 0);
          for(MKL_LONG ip=0; ip<TRAJ.Np; ip++)
            {
              CONNECT.HASH[ip].fprint_LONG_skip_transpose_LIMROWS(filename_HASH.c_str(), 1, CONNECT.TOKEN[ip]);
              CONNECT.weight[ip].fprint_LONG_skip_transpose_LIMROWS(filename_weight.c_str(), 1, CONNECT.TOKEN[ip]);
            }
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
  //   for(MKL_LONG i=0; i<TRAJ.Np; i++)
  //     mkl_free(R_minimum_vec_boost[i]); // RDIST
  // mkl_free(R_minimum_vec_boost); // RDIST
  // mkl_free(R_minimum_distance_boost); // RDIST
  mkl_free(tmp_index_vec);
  mkl_free(dCDF_U);
  mkl_free(INDEX_dCDF_U);
  mkl_free(dCDF_TOKEN);  
  for(MKL_LONG i=0; i<N_THREADS_BD; i++)
    gsl_rng_free(r_boost_arr[i]); // for boosting

  mkl_free(r_boost_arr);
  mkl_free(IDX_ARR);
  return 0;
}
/*
 * Local variables:
 * compile-command: "icpc -openmp -O2 -Wall -mkl -o Brownian_simulation lib/trajectory.cpp lib/read_file_condition.cpp lib/time_evolution.cpp lib/matrix.cpp lib/potential.cpp lib/connectivity.cpp lib/association.cpp lib/handle_association.cpp lib/geometry.cpp lib/random.cpp lib/parallel.cpp src/Brownian_simulation.cpp -L/usr/local/include/ -L/usr/local/lib/ -lgsl -lm"
 * End:
 */

/*
  -openmp: OpenMP option
  -ipo: interprocedure optimization
  -O2: optimization level 2
  -Wall: rigirous compile option (no warning point is allowed)
  -mkl: Math Kernel Librar (Intel)
  -lgsl: GSL library
  -lm: Math library
*/

