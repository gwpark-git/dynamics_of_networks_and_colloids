#include "stochastic_HEUR_flowers.h"


MKL_LONG stochastic_simulation_HEUR_flowers(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, RECORD_DATA& DATA, COND& given_condition)
{

  using namespace std;
  printf("STARTING_MAIN_ROOT\n"); 
  double time_st_simulation = dsecnd(); 
  MKL_LONG N_steps_block = atol(given_condition("N_steps_block").c_str());
  
  MKL_LONG N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
  MKL_LONG N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());
  printf("THREAD_SETTING: %ld ... ", N_THREADS_BD);
  mkl_set_num_threads(N_THREADS_BD);
  double time_MC = 0.;
  double time_LV = 0.;
  double time_AN = 0.;
  double time_file = 0.;
  printf("DONE\n");
  printf("GENERATING BOOSTING VECTORS\n");
  LOCK LOCKER(TRAJ.Np);
  MKL_LONG *tmp_index_vec = (MKL_LONG*) mkl_malloc(TRAJ.N_dimension*sizeof(MKL_LONG), BIT);
  MATRIX *vec_boost_Nd_parallel = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); 
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(TRAJ.N_dimension, 1, 0.);
    }
  RDIST R_boost(given_condition);
  
  MKL_LONG cnt_arr[6] = {0};
  MKL_LONG &cnt_cancel = cnt_arr[INDEX_MC::CANCEL], &cnt_add = cnt_arr[INDEX_MC::ADD], &cnt_del = cnt_arr[INDEX_MC::OPP_DEL], &cnt_mov = cnt_arr[INDEX_MC::MOV], &cnt_lock = cnt_arr[INDEX_MC::LOCK], &cnt_SS = cnt_arr[5];
  cnt_add = CONNECT.N_ASSOCIATION;
  
  printf("DONE\n"); 
  printf("FORCE VECTOR GENERATING ... "); 
  MATRIX *force_spring = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_repulsion = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_random = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      force_spring[i].initial(TRAJ.N_dimension, 1, 0.);
      force_repulsion[i].initial(TRAJ.N_dimension, 1, 0.);
      force_random[i].initial(TRAJ.N_dimension, 1, 0.);
    }
  printf("DONE\n"); 
  if(given_condition("MC_LOG") == "TRUE")
    {
      DATA.MC_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dCDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << '\t' << "16_cnt_lock" << endl;
    }
  printf("SET SIMULATION PARAMETERS ...");
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());

  MATRIX energy(1, 6, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  // MATRIX energy(1, 12, 0.);

  ANALYSIS::CAL_ENERGY_R_boost(POTs, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);
                               
  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  MKL_LONG N_basic = TRAJ.rows;

  double tolerance_association = atof(given_condition("tolerance_association").c_str());
  printf("DONE\n");

  double dt_1 = 0., dt_2 = 0., dt_3 = 0., dt_4 = 0., dt_5 = 0., dt_6 = 0., dt_7 = 0.;
  double dt_det_pdf = 0.;
  double dt_rdist = 0., dt_pdf = 0., dt_sort = 0.;
  double time_MC_1 = 0., time_MC_2 = 0., time_MC_3 = 0., time_MC_4 = 0., time_MC_5 = 0., time_MC_6 = 0., time_MC_7 = 0., time_MC_8 = 0.;

  printf("GENERATING RANDOM VECTOR BOOST ... ");
  RNG_BOOST RNG(given_condition);
  
  INDEX_MC *IDX_ARR = (INDEX_MC*) mkl_malloc(N_THREADS_BD*sizeof(INDEX_MC), BIT);
  for(MKL_LONG i=0; i<N_THREADS_SS; i++)
    {
      IDX_ARR[i].initial();
      IDX_ARR[i].set_initial_variables();
    }
  
  printf("DONE\n");
  printf("START SIMULATION\n");
  // MKL_LONG N_associations = 0;

  double N_diff = 0.;
  MKL_LONG N_tot_associable_chain = TRAJ.Np*atoi(given_condition("N_chains_per_particle").c_str());
  
  for(MKL_LONG t = 0; t<Nt-1; t++)
    {
      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
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
          // to get all pairs of particles in cell
          R_boost.compute_RDIST_particle(index_particle, TRAJ, index_t_now);
        } // index_particle
      dt_rdist += dsecnd() - time_st_rdist;
      double time_st_MC = dsecnd();
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
              double time_st_pdf = dsecnd();
              
              /*
                The following will check only for the existing brdige.
                Since the bridges will only occurred withing neighbouring cell-list, this bridge chain will not violate the condition for cell-list.
              */
#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, R_boost, vec_boost_Nd_parallel) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
              for(MKL_LONG i=0; i<TRAJ.Np; i++)
                {
                  CONNECT.update_CHAIN_SUGGESTION_MAP_particle(i, POTs, R_boost);
                }
            }  // else for MC_renewal check
          // #pragma omp for
#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, R_boost, N_THREADS_BD, dt_pdf, dt_sort) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
          for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
            {
              CONNECT.update_ASSOCIATION_MAP_particle(index_particle, POTs, R_boost);
            } // index_particle, local parallel
#pragma omp parallel for default(none) shared(given_condition, DATA, TRAJ, POTs, CONNECT, CHAIN, LOCKER, IDX_ARR, index_t_now, vec_boost_Nd_parallel, R_boost, dt_rdist, dt_pdf, dt_sort, cnt_arr, cnt_add, cnt_del, cnt_mov, cnt_cancel, cnt_lock, cnt_SS, N_steps_block, RNG, cnt, N_THREADS_SS, N_tot_associable_chain) private(time_MC_1, time_MC_2, time_MC_3, time_MC_4, time_MC_5, time_MC_6, time_MC_7, time_MC_8) num_threads(N_THREADS_SS) if(N_THREADS_SS > 1) reduction(+:dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7)
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
              index_itself = RANDOM::return_LONG_INT_rand_boost(RNG.BOOST_SS[it], TRAJ.Np);
              // choice for selected chain end
              double rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(RNG.BOOST_SS[it]);
              time_MC_2 = dsecnd();
              index_hash_attached_bead = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF); 
              index_attached_bead = CONNECT.HASH[index_itself](index_hash_attached_bead); 
              time_MC_3 = dsecnd();
              // choice for behaviour of selected chain end
              double rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(RNG.BOOST_SS[it]);
              // the PDF is already computed in the previous map
              time_MC_4 = dsecnd();
              
              MKL_LONG k = SEARCHING::backtrace_cell_list(CONNECT.dCDF_ASSOCIATION[index_itself], CONNECT.TOKEN_ASSOCIATION[index_itself], rolling_dCDF_U, index_itself, R_boost);
              index_new_attached_bead = CONNECT.INDEX_ASSOCIATION[index_itself](k);
              time_MC_5 = dsecnd();
              MKL_LONG IDENTIFIER_ACTION = TRUE; // it can be 1 (IDX.ADD) but just true value
              MKL_LONG IDENTIFIER_LOCKING = FALSE;
              if(N_THREADS_SS > 1)
                {
#pragma omp critical(LOCKING)  // LOCKING is the name for this critical blocks
                  {
                    cnt_SS ++;
                    /*
                      On the omp critical region, the block will work only one thread.
                      If the other thread reaching this reason while there is one thread already working on this block, then the reached thread will wait until finishing the job of the other thread.
                      This benefits to identify the working beads index on this case, since the 
                    */
                    // CHECKING

                    for(MKL_LONG I_BEADS = 0; I_BEADS < 3; I_BEADS++)
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
                        // for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++) // N_THREADS_SS > 1 is not necessary
                        for(MKL_LONG I_BEADS = 0; I_BEADS < 3; I_BEADS++) 
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
                      double rolling_transition = RANDOM::return_double_rand_SUP1_boost(RNG.BOOST_SS[it]);
                      if (rolling_transition < tpa)
                        {
                          IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX_ARR[it]);
                        }
                      else
                        IDENTIFIER_ACTION = IDX_ARR[it].CANCEL;
                    }

                  ACTION::ACT(index_t_now, POTs, CONNECT, IDX_ARR[it], R_boost.Rsca, IDENTIFIER_ACTION);
                  if(CHAIN.INITIALIZATION)
                    {
                      // note that this is affected by LOCKING scheme for parallelism of stochastic simulation part
                      // hence the tracking individual chain is not affected by the SS parallelisation scheme.
                      CHAIN.TRACKING_ACTION(CONNECT, IDENTIFIER_ACTION, IDX_ARR[it]); // it will track individual chain information
                    }
                  
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
#pragma omp atomic // it provide lower overhead compared with critical directive
                  cnt_arr[IDENTIFIER_ACTION] ++;

                  if (given_condition("MC_LOG") == "TRUE")
                    {
#pragma omp critical(MC_LOG)
                      {
                        MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
			  
                        // MKL_LONG count_N_associagtions = cnt_add - cnt_del;
                        {
                          DATA.MC_LOG << cnt << '\t' << index_itself << '\t' << setprecision(7) << rolling_dCDF<< '\t'  << index_attached_bead << '\t'  << index_new_attached_bead<< '\t'  << setprecision(7) << rolling_dCDF_U<< '\t'  << k<< '\t'  << index_new_attached_bead << '\t'  << CONNECT.TOKEN[index_itself]<< '\t'<< CONNECT.N_CONNECTED_ENDS(index_itself) << '\t' << CONNECT.weight[index_itself](0) <<'\t' <<  total_bonds << '\t'  << cnt_add<< '\t'  << cnt_mov<< '\t'  << cnt_del<< '\t'  << cnt_cancel << '\t' << cnt_lock << endl;
                        }
                        // FILE_LOG << boost::format("%10d\t%4d\t")
                      } // MC_LOG
                    }
                } // if(!IDENTIFIER_LOCKING)
            } // for loop (ASSOCIATION)
        } // region for topological update
      double time_end_MC = dsecnd();

#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next, R_boost, vec_boost_Nd_parallel, force_spring, force_repulsion, force_random, RNG, N_THREADS_BD, given_condition) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      for (MKL_LONG i=0; i<TRAJ.Np; i++)
        {
          MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
          
          force_spring[i].set_value(0);
          force_repulsion[i].set_value(0);
          force_random[i].set_value(0);

          INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(POTs, CONNECT, force_spring[i], i, R_boost.Rvec, R_boost.Rsca);
          INTEGRATOR::EULER::cal_repulsion_force_R_boost(POTs, force_repulsion[i], i, R_boost);
          INTEGRATOR::EULER::cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]); 
          
          for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
            {
              TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k)) + sqrt(TRAJ.dt)*force_random[i](k);
            }
        }
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC
      MKL_LONG N_associations = cnt_add - cnt_del;
      double time_end_LV = dsecnd();
      double time_end_AN = time_end_LV;
      if(t%N_skip==0)
        {
          time_end_LV = dsecnd();
          energy(0) = TRAJ(index_t_now);
          energy(4) = (double)N_associations; // number of associations
          energy(5) = dsecnd() - time_st_simulation;
          ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY_R_boost(POTs, CONNECT, energy, (TRAJ.c_t - 1.)*TRAJ.dt, vec_boost_Nd_parallel[0], R_boost);
          time_end_AN = dsecnd();
          double total_dt = dt_1 + dt_2 + dt_3 + dt_4 + dt_5 + dt_6 + dt_7;
          double total_dt_pdf = dt_rdist + dt_pdf + dt_sort;
          double total_time = time_MC + time_LV + time_AN + time_file + total_dt_pdf;
          double dt_pdf_all = dt_pdf + dt_sort;
          
          printf("##### STEPS = %ld\tTIME = %8.6e tau_0\tENERGY = %6.3e (computing time = %4.3e)\n", TRAJ.c_t, TRAJ(index_t_now)/atof(given_condition("Rt").c_str()), energy(1), energy(5));
          printf("time consuming: MC = %3.2e (%3.1f), LV = %3.2e (%3.1f), AN = %3.2e (%3.1f), FILE = %3.2e (%3.1f), DIST = %3.2e (%3.1f)\n", time_MC, time_MC*100/total_time, time_LV, time_LV*100/total_time, time_AN, time_AN*100/total_time, time_file, time_file*100/total_time, total_dt_pdf, total_dt_pdf*100/total_time);
          printf("MC: all pdf = %3.2e (%3.1f), basic_random = %3.2e (%3.1f), getting_hash = %3.2e (%3.1f), det_jump = %3.2e (%3.1f), new_end = %3.2e (%3.1f), LOCKING = %3.2e (%3.1f), action = %3.2e (%3.1f), update = %3.2e (%3.1f)\n", dt_pdf_all, dt_pdf_all*100./total_dt, dt_1, dt_1*100/total_dt, dt_2, dt_2*100./total_dt, dt_3, dt_3*100./total_dt, dt_4, dt_4*100./total_dt, dt_5, dt_5*100./total_dt, dt_6, dt_6*100./total_dt, dt_7, dt_7*100./total_dt);
          printf("DIST: rdist = %6.3e (%3.1f), computing pdf = %6.3e (%3.1f), sorting pdf = %6.3e (%3.1f)\n", dt_rdist, 100.*dt_rdist/total_dt_pdf, dt_pdf, 100.*dt_pdf/total_dt_pdf, dt_sort, dt_sort*100./total_dt_pdf);
          printf("CHECK LAST STATISTICS: N_tot_asso = %ld, NAS = %ld, fraction=%4.3f, cnt_lock/cnt_normal=%3.2e ####\n\n", N_tot_associable_chain, N_associations, N_associations/(double)N_tot_associable_chain, cnt_lock/(double)cnt_SS);

          TRAJ.fprint_row(DATA.traj, index_t_now);
          energy.fprint_row(DATA.ener, 0);
          for(MKL_LONG ip=0; ip<TRAJ.Np; ip++)
            {
              CONNECT.HASH[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.hash, 1, CONNECT.TOKEN[ip]);
              CONNECT.weight[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.weight, 1, CONNECT.TOKEN[ip]);
            }
          if(CHAIN.INITIALIZATION)
            CHAIN.write(DATA.chain);
        }
      

      double time_end_save = dsecnd();

      time_MC += time_end_MC - time_st_MC;
      time_LV += time_end_LV - time_end_MC;
      time_AN += time_end_AN - time_end_LV;
      time_file += time_end_save - time_end_AN;
    }

  double time_simulation = dsecnd() - time_st_simulation;
  printf("Total simulation time = %6.3e\n", time_simulation);
  mkl_free(vec_boost_Nd_parallel);
  mkl_free(force_spring);
  mkl_free(force_repulsion);
  mkl_free(force_random);
  mkl_free(tmp_index_vec);
  mkl_free(IDX_ARR);
  return 0;
}

