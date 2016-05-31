#include "stochastic_HEUR_flowers.h"

MKL_LONG stochastic_simulation_HEUR_flowers(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, RECORD_DATA& DATA, COND& given_condition)
{

  using namespace std;
  double time_st_simulation = dsecnd();
  
  printf("### STOCHASTIC SIMULATION FOR HEUR FLOWER ###\n");
  printf("SETTING simulation environment ...");
  TEMPORAL_VARIABLE_HEUR VAR(given_condition, TRAJ.rows);
  mkl_set_num_threads(VAR.N_THREADS_BD);
  LOCK LOCKER(TRAJ.Np);
  RDIST R_boost(given_condition);
  
  MKL_LONG &cnt_cancel = VAR.cnt_arr[INDEX_MC::CANCEL], &cnt_add = VAR.cnt_arr[INDEX_MC::ADD], &cnt_del = VAR.cnt_arr[INDEX_MC::OPP_DEL], &cnt_mov = VAR.cnt_arr[INDEX_MC::MOV], &cnt_lock = VAR.cnt_arr[INDEX_MC::LOCK], &cnt_SS = VAR.cnt_arr[5];
  cnt_add = CONNECT.N_ASSOCIATION;
  
  if(given_condition("MC_LOG") == "TRUE")
    {
      DATA.MC_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dCDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << '\t' << "16_cnt_lock" << endl;
    }
  printf("SET SIMULATION PARAMETERS ...");

  MATRIX energy(1, 6, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]

  VAR.time_AN +=
    ANALYSIS::CAL_ENERGY_R_boost(POTs, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);

  RNG_BOOST RNG(given_condition);
  
  // INDEX_MC *IDX_ARR = (INDEX_MC*) mkl_malloc(VAR.N_THREADS_SS*sizeof(INDEX_MC), BIT);
  INDEX_MC *IDX_ARR = new INDEX_MC [VAR.N_THREADS_SS];
  for(MKL_LONG i=0; i<VAR.N_THREADS_SS; i++)
    {
      IDX_ARR[i].initial();
      IDX_ARR[i].set_initial_variables();
    }

  printf("DONE\nSTART SIMULATION\n\n");
  
  for(MKL_LONG t = 0; t<VAR.Nt-1; t++)
    {
      MKL_LONG index_t_now = t % VAR.N_basic;
      MKL_LONG index_t_next = (t+1) % VAR.N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
      MKL_LONG cnt = 1;
      VAR.time_DIST +=
        REPULSIVE_BROWNIAN::OMP_compute_RDIST(TRAJ, index_t_now, R_boost, VAR.tmp_index_vec, VAR.N_THREADS_BD);
      VAR.time_SS +=
        OMP_SS_topological_time_evolution(t, CONNECT, CHAIN, POTs, R_boost, RNG, IDX_ARR, DATA, given_condition, LOCKER, VAR);
      VAR.time_LV +=
        OMP_time_evolution_Euler(TRAJ, index_t_now, index_t_next, CONNECT, POTs, R_boost, VAR.vec_boost_Nd_parallel, VAR.force_repulsion, VAR.force_random, VAR.force_spring, RNG, VAR.N_THREADS_BD, given_condition, VAR);
      
      VAR.time_LV +=
        GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC

      VAR.N_associations = cnt_add - cnt_del;
      
      if(t%VAR.N_skip==0)
        {
          VAR.time_AN +=
            ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY_R_boost(POTs, CONNECT, energy, (TRAJ.c_t - 1.)*TRAJ.dt, VAR.vec_boost_Nd_parallel[0], R_boost);

          energy(4) = (double)VAR.N_associations; // number of associations
          energy(5) = dsecnd() - time_st_simulation; // record computation time

          VAR.time_file +=
            record_simulation_data(DATA, TRAJ, CONNECT, CHAIN, energy, index_t_now); // neeed review

          VAR.time_RECORDED +=
            report_simulation_info(TRAJ, energy, VAR);
        }
    }

  double time_simulation = dsecnd() - time_st_simulation;
  printf("Total simulation time = %6.3e, recorded time = %6.3e (record/total = %3.1f)\n", time_simulation, VAR.time_RECORDED, VAR.time_RECORDED/time_simulation);
  // mkl_free(IDX_ARR);
  delete[] IDX_ARR;
  return 0;
}

double record_simulation_data(RECORD_DATA& DATA, TRAJECTORY& TRAJ, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, MATRIX& energy, const MKL_LONG index_t_now) 
{
  double time_st = dsecnd();
  REPULSIVE_BROWNIAN::record_simulation_data(DATA, TRAJ, energy, index_t_now);

  for(MKL_LONG ip=0; ip<TRAJ.Np; ip++)
    {
      CONNECT.HASH[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.hash, 1, CONNECT.TOKEN[ip]);
      CONNECT.weight[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.weight, 1, CONNECT.TOKEN[ip]);
    }

  if(CHAIN.INITIALIZATION)
    CHAIN.write(DATA.chain);
  return dsecnd() - time_st;
}

double report_simulation_info(TRAJECTORY& TRAJ, MATRIX& energy, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double total_time = VAR.time_SS + VAR.time_LV + VAR.time_AN + VAR.time_file + VAR.time_DIST;
  double sum_time_LV = VAR.time_LV_init + VAR.time_LV_force + VAR.time_LV_update;
  printf("##### STEPS = %ld\tTIME = %8.6e tau_B\tENERGY = %6.3e (computing time: %4.3e)\n", TRAJ.c_t, VAR.simulation_time, energy(1), energy(5));
  printf("time consuming: LV = %3.2e (%3.1f%%), SS = %3.2e (%3.1f%%), AN = %3.2e (%3.1f%%), FILE = %3.2e (%3.1f%%), DIST = %3.2e (%3.1f%%)\n", VAR.time_LV, VAR.time_LV*100/total_time, VAR.time_SS, VAR.time_SS*100/total_time, VAR.time_AN, VAR.time_AN*100/total_time, VAR.time_file, VAR.time_file*100/total_time, VAR.time_DIST, VAR.time_DIST*100/total_time);
  printf("LV: init = %3.1f%%, force = %3.1f%% (repulsion = %3.1f%%, random = %3.1f%%, connector = %3.1f%%), update = %3.1f%%\n", 100.*VAR.time_LV_init/sum_time_LV, 100.*VAR.time_LV_force/sum_time_LV, 100.*VAR.time_LV_force_repulsion/VAR.time_LV_force, 100.*VAR.time_LV_force_random/VAR.time_LV_force, 100.*VAR.time_LV_force_connector/VAR.time_LV_force, 100.*VAR.time_LV_update/sum_time_LV);
  printf("SS: index process = %3.1f%%, LOCKING = %3.1f%%, check_dissociation = %3.1f%%, transition = %3.1f%%, update info. = %3.1f%%, pre-ASSOCIATION = %3.1f%%, pre-SUGGESTION = %3.1f%%)  \n", 100.*VAR.time_SS_index/VAR.time_SS, 100.*VAR.time_SS_LOCK/VAR.time_SS, 100.*VAR.time_SS_check/VAR.time_SS, 100.*VAR.time_SS_transition/VAR.time_SS, 100.*VAR.time_SS_update_info/VAR.time_SS, 100.*VAR.time_SS_update_ASSOCIATION_MAP/VAR.time_SS, 100.*VAR.time_SS_update_CHAIN_SUGGESTION/VAR.time_SS);
  // 	 100.*VAR.time_SS_update/VAR.time_SS, VAR.time_SS_update_ASSOCIATION_MAP, VAR.time_SS_update_ASSOCIATION_MAP/VAR.time_SS_update, VAR.time_SS_update_CHAIN_SUGGESTION, VAR.time_SS_update_CHAIN_SUGGESTION/VAR.time_SS_update);
  printf("CHECK LAST STATISTICS: N_tot_asso = %ld, NAS = %ld, fraction=%4.3f ####\n\n", VAR.N_tot_associable_chain, VAR.N_associations, VAR.N_associations/(double)VAR.N_tot_associable_chain);
  return total_time;
}


double OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next, ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, MATRIX* vec_boost_Nd_parallel, MATRIX* force_repulsion, MATRIX* force_random, MATRIX* force_spring, RNG_BOOST& RNG, const MKL_LONG N_THREADS_BD, COND& given_condition, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double time_st = dsecnd();
  double time_LV_init = 0., time_LV_force = 0., time_LV_update = 0.;
  double time_LV_force_repulsion = 0., time_LV_force_random = 0., time_LV_force_connector = 0.;

#pragma omp parallel for default(none) shared(TRAJ, CONNECT, POTs, index_t_now, index_t_next, R_boost, vec_boost_Nd_parallel, force_repulsion, force_random, force_spring, RNG, given_condition) num_threads(N_THREADS_BD) reduction(+: time_LV_init, time_LV_force, time_LV_update, time_LV_force_repulsion, time_LV_force_random, time_LV_force_connector) if(N_THREADS_BD > 1)
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
      double time_st_init = dsecnd();
      force_repulsion[i].set_value(0);
      force_random[i].set_value(0);
      force_spring[i].set_value(0);
      double time_st_force = dsecnd();
      time_LV_init += time_st_force - time_st_init;
      time_LV_force_connector +=
        INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(POTs, CONNECT, force_spring[i], i, R_boost.Rvec, R_boost.Rsca);
      time_LV_force_repulsion +=
        INTEGRATOR::EULER::cal_repulsion_force_R_boost(POTs, force_repulsion[i], i, R_boost);
      time_LV_force_random +=
        // INTEGRATOR::EULER::cal_random_force_boost_simplified(POTs, force_random[i], RNG.BOOST_BD[it]);
        INTEGRATOR::EULER::cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]);
      
      double time_st_update = dsecnd();
      time_LV_force += time_st_update - time_st_force;
      
      for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k)) + sqrt(TRAJ.dt)*force_random[i](k);
        }
      
      time_LV_update += dsecnd() - time_st_update;
    }
  VAR.time_LV_init += time_LV_init;
  VAR.time_LV_force += time_LV_force;
  VAR.time_LV_update += time_LV_update;

  VAR.time_LV_force_repulsion += time_LV_force_repulsion;
  VAR.time_LV_force_random += time_LV_force_random;
  VAR.time_LV_force_connector += time_LV_force_connector;
  return dsecnd() - time_st;
}

double write_MC_LOG_if_TRUE(bool flag_MC_LOG, RECORD_DATA& DATA, ASSOCIATION& CONNECT, const INDEX_MC& IDX, const MKL_LONG cnt, const MKL_LONG* cnt_arr, const double rolling_dCDF, const double rolling_dCDF_U) 
// const MKL_LONG index_itself, const double rolling_dCDF, const MKL_LONG index_attached_bead, const MKL_LONG index_new_attached_bead, const double rolling_dCDF_U, const MKL_LONG k, MKL_LONG* cnt_arr)
{
  double time_st = dsecnd();
  if(flag_MC_LOG)
    {
      MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
			  
      // MKL_LONG count_N_associagtions = cnt_add - cnt_del;
      {
        DATA.MC_LOG << cnt << '\t' << IDX.beads[CONNECT.flag_itself] << '\t' << setprecision(7) << rolling_dCDF<< '\t'  << IDX.beads[CONNECT.flag_hash_other] << '\t'  << IDX.beads[CONNECT.flag_other] << '\t'  << setprecision(7) << rolling_dCDF_U<< '\t'  << IDX.beads[CONNECT.flag_hash_backtrace] << '\t'  << IDX.beads[CONNECT.flag_new] << '\t'  << CONNECT.TOKEN[IDX.beads[CONNECT.flag_itself]]<< '\t'<< CONNECT.N_CONNECTED_ENDS(IDX.beads[CONNECT.flag_itself]) << '\t' << CONNECT.weight[IDX.beads[CONNECT.flag_itself]](0) <<'\t' <<  total_bonds << '\t'  << cnt_arr[INDEX_MC::ADD]<< '\t'  << cnt_arr[INDEX_MC::MOV]<< '\t'  << cnt_arr[INDEX_MC::OPP_DEL]<< '\t'  << cnt_arr[INDEX_MC::CANCEL] << '\t' << cnt_arr[INDEX_MC::LOCK] << endl;
      }
      // FILE_LOG << boost::format("%10d\t%4d\t")
    } // MC_LOG
  return dsecnd() - time_st;
}

double transition_single_chain_end(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, CHAIN_HANDLE& CHAIN, RDIST& R_boost, INDEX_MC& IDX, const MKL_LONG IDENTIFIER_ACTION)
{
  double time_st = dsecnd();
  ACTION::ACT(POTs, CONNECT, IDX, R_boost.Rsca, IDENTIFIER_ACTION);
  if(CHAIN.INITIALIZATION)
    {
      // note that this is affected by LOCKING scheme for parallelism of stochastic simulation part
      // hence the tracking individual chain is not affected by the SS parallelisation scheme.
      CHAIN.TRACKING_ACTION(CONNECT, IDENTIFIER_ACTION, IDX); // it will track individual chain information
    }
  return dsecnd() - time_st;
}

double check_dissociation_probability(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, INDEX_MC& IDX, gsl_rng* RNG_BOOST_SS_IT, MKL_LONG& IDENTIFIER_ACTION)
{
  double time_st = dsecnd();
  MKL_LONG index_itself = IDX.beads[CONNECT.flag_itself];
  MKL_LONG index_attached_bead = IDX.beads[CONNECT.flag_other];
  
  double distance_exist_bridge = R_boost.Rsca[index_itself](index_attached_bead);
  double tpa = POTs.transition(distance_exist_bridge, POTs.f_connector(distance_exist_bridge, POTs.force_variables), POTs.force_variables);
  if (tpa == 1.0)
    {
      IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX);
    }
  else
    {
      double rolling_transition = RANDOM::return_double_rand_SUP1_boost(RNG_BOOST_SS_IT);
      if (rolling_transition < tpa)
        {
          IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX);
        }
      else
        IDENTIFIER_ACTION = IDX.CANCEL;
    }
  return dsecnd() - time_st;
}

double LOCKING_PARALLEL(LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR, const INDEX_MC& IDX, MKL_LONG& IDENTIFIER_ACTION, MKL_LONG& IDENTIFIER_LOCKING)
{
  double time_st = dsecnd();
  /*
    On the omp critical region, the block will work only one thread.
    If the other thread reaching this reason while there is one thread already working on this block, then the reached thread will wait until finishing the job of the other thread.
    This benefits to identify the working beads index on this case, since the 
  */
  // CHECKING

  for(MKL_LONG I_BEADS = 0; I_BEADS < 3; I_BEADS++)
    {
      if(LOCKER(IDX.beads[I_BEADS]))
        {
          IDENTIFIER_ACTION = IDX.CANCEL;
          IDENTIFIER_LOCKING = TRUE;
          break;
        }
    }

  // this is LOCKING procedure
  if(!IDENTIFIER_LOCKING)
    {
      VAR.cnt_SS ++;
      // cnt++;  // preventing LOCKING affect to the IDENTIFICATION of stochastic balance
      // for(MKL_LONG I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++) // N_THREADS_SS > 1 is not necessary
      for(MKL_LONG I_BEADS = 0; I_BEADS < 3; I_BEADS++) 
        {
          LOCKER(IDX.beads[I_BEADS]) = TRUE;
        }
    }
  else
    {
      VAR.cnt_arr[INDEX_MC::LOCK] ++;
    }
  return dsecnd() - time_st;
}

double release_LOCKING(LOCK& LOCKER, INDEX_MC& IDX)
{
  double time_st = dsecnd();
  for(MKL_LONG I_BEADS = 0; I_BEADS < 3; I_BEADS++)
    {
      LOCKER(IDX.beads[I_BEADS]) = FALSE;
    }
  return dsecnd() - time_st;
}

// double micelle_selection(ASSOCIATION& CONNECT, RNG_BOOST& RNG, MKL_LONG& index_itself, MKL_LONG& index_hash_attached_bead, MKL_LONG& index_attached_bead, MKL_LONG& index_new_attached_bead, const MKL_LONG index_thread, RDIST& R_boost, TEMPORAL_VARIABLE_HEUR& VAR)
double micelle_selection(ASSOCIATION& CONNECT, gsl_rng* RNG_BOOST_SS_IT, INDEX_MC& IDX, RDIST& R_boost, TEMPORAL_VARIABLE_HEUR& VAR, double& rolling_dCDF, double& rolling_dCDF_U)
{
  double time_st = dsecnd();
  // index_itself = RANDOM::return_LONG_INT_rand_boost(RNG.BOOST_SS[index_thread], VAR.Np);
  // double rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(RNG.BOOST_SS[index_thread]);
  // index_hash_attached_bead = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF);
  // index_attached_bead = CONNECT.HASH[index_itself](index_hash_attached_bead);
  // double rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(RNG.BOOST_SS[index_thread]);
  // MKL_LONG k = SEARCHING::backtrace_cell_list(CONNECT.dCDF_ASSOCIATION[index_itself], CONNECT.TOKEN_ASSOCIATION[index_itself], rolling_dCDF_U, index_itself, R_boost);
  // index_new_attached_bead = CONNECT.INDEX_ASSOCIATION[index_itself](k);
  
  IDX.beads[CONNECT.flag_itself] = RANDOM::return_LONG_INT_rand_boost(RNG_BOOST_SS_IT, VAR.Np);
  // choice for selected chain end
  rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(RNG_BOOST_SS_IT);
  IDX.beads[CONNECT.flag_hash_other] = CONNECT.GET_INDEX_HASH_FROM_ROLL(IDX.beads[CONNECT.flag_itself], rolling_dCDF); 
  IDX.beads[CONNECT.flag_other] = CONNECT.HASH[IDX.beads[CONNECT.flag_itself]](IDX.beads[CONNECT.flag_hash_other]); 
  // choice for behaviour of selected chain end
  rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(RNG_BOOST_SS_IT);
  // the PDF is already computed in the previous map
  IDX.beads[CONNECT.flag_hash_backtrace] = SEARCHING::backtrace_cell_list(CONNECT.dCDF_ASSOCIATION[IDX.beads[CONNECT.flag_itself]], CONNECT.TOKEN_ASSOCIATION[IDX.beads[CONNECT.flag_itself]], rolling_dCDF_U, IDX.beads[CONNECT.flag_itself], R_boost);
  IDX.beads[CONNECT.flag_new] = CONNECT.INDEX_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](IDX.beads[CONNECT.flag_hash_backtrace]);
  // IDX.beads[CONNECT.flag_new] = CONNECT.INDEX_ASSOCIATION[IDX.beads[CONNECT.flag_itself]][k];
  
  return dsecnd() - time_st;
}

double OMP_SS_update_topology(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, CHAIN_HANDLE& CHAIN, RNG_BOOST& RNG, RECORD_DATA& DATA, INDEX_MC* IDX_ARR, LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double time_st = dsecnd();
  double time_SS_index = 0., time_SS_LOCK = 0., time_SS_check = 0., time_SS_transition = 0., time_SS_update_info = 0.;
#pragma omp parallel for default(none) shared(DATA, POTs, CONNECT, CHAIN, LOCKER, IDX_ARR, VAR, R_boost, RNG) num_threads(VAR.N_THREADS_SS) reduction(+: time_SS_index, time_SS_LOCK, time_SS_check, time_SS_transition, time_SS_update_info) if(VAR.N_THREADS_SS > 1) 
  // reduction(+:dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7)
  for(MKL_LONG tp=0; tp<VAR.N_tot_associable_chain; tp++)
    {
      /*
        'it' have the identity number for current thread. Then, the reference variables IDX and r_boost will be used in order to usability and readability. In this case, IDX, r_boost is just reference of existing one, but IDX and r_boost itself is local reference variables which will varied thread to thread
      */
      MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects

      // name definition for convenience
      MKL_LONG &index_itself = IDX_ARR[it].beads[CONNECT.flag_itself];
      MKL_LONG &index_attached_bead = IDX_ARR[it].beads[CONNECT.flag_other];
      MKL_LONG &index_new_attached_bead = IDX_ARR[it].beads[CONNECT.flag_new];
      MKL_LONG &index_hash_attached_bead = IDX_ARR[it].beads[CONNECT.flag_hash_other];

      double rolling_dCDF = 0., rolling_dCDF_U = 0.;
      
      time_SS_index +=
        micelle_selection(CONNECT, RNG.BOOST_SS[it], IDX_ARR[it], R_boost, VAR, rolling_dCDF, rolling_dCDF_U);
      // micelle_selection(CONNECT, RNG, index_itself, index_hash_attached_bead, index_attached_bead, index_new_attached_bead, it, R_boost, VAR);
          
      MKL_LONG IDENTIFIER_ACTION = TRUE; // it can be 1 (IDX.ADD) but just true value
      MKL_LONG IDENTIFIER_LOCKING = FALSE;
      // locking
      if(VAR.N_THREADS_SS > 1)
        {
#pragma omp critical (LOCKING)  // LOCKING is the name for this critical blocks
          {
          time_SS_LOCK += // this is differ from time_SS_LOCK since it is inside critical directive
            LOCKING_PARALLEL(LOCKER, VAR, IDX_ARR[it], IDENTIFIER_ACTION, IDENTIFIER_LOCKING);
        }
    }
      // Note that the critical region only applicable with single thread while the others will be used in parallel regime.
      // In addition, the gap for passing the critical region will tune further gaps, then the computation speed for passing critical region will not be real critical issue.
      if(!IDENTIFIER_LOCKING) 
        {
          // This block only compute when the thread is NOT LOCKED
          time_SS_check +=
            check_dissociation_probability(CONNECT, POTs, R_boost, IDX_ARR[it], RNG.BOOST_SS[it], IDENTIFIER_ACTION);

          time_SS_transition +=
            transition_single_chain_end(CONNECT, POTs, CHAIN, R_boost, IDX_ARR[it], IDENTIFIER_ACTION);

          time_SS_update_info +=
            ACTION::UPDATE_INFORMATION(CONNECT, IDX_ARR[it], VAR.cnt_arr, IDENTIFIER_ACTION);

          // UNLOCKING
          // The critical directive is no more necessarly since only one thread visited each beads
          if(VAR.N_THREADS_SS > 1)
            time_SS_LOCK +=
              release_LOCKING(LOCKER, IDX_ARR[it]);
          
#pragma omp critical (COUNTING)
          {
            /*
              critical(COUNTING) blocks:
              This is counting the action information that will be used for the future.
              Notice that the writing MC_LOG file is inside of this COUNTING critical directive, since all the information should be the same for writing (temporal)
            */
            
            VAR.cnt_arr[IDENTIFIER_ACTION] ++;

            // write_MC_LOG_if_TRUE(VAR.MC_LOG, DATA, CONNECT, VAR.cnt_SS, IDX_ARR[it], index_itself, rolling_dCDF, index_attached_bead, index_new_attached_bead, rolling_dCDF_U, VAR.cnt_arr);
            write_MC_LOG_if_TRUE(VAR.MC_LOG, DATA, CONNECT, IDX_ARR[it], VAR.cnt_SS + 1, VAR.cnt_arr, rolling_dCDF, rolling_dCDF_U);
          }

        } // if(!IDENTIFIER_LOCKING)
}

  // adding the measured time to the variable structure
  VAR.time_SS_index += time_SS_index;
VAR.time_SS_LOCK += time_SS_LOCK;
VAR.time_SS_check += time_SS_check;
VAR.time_SS_transition += time_SS_transition;
VAR.time_SS_update_info += time_SS_update_info;

return dsecnd() - time_st;
}

double OMP_SS_update_STATISTICS(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double time_st = dsecnd();
  MKL_LONG Np = VAR.Np;
  if (VAR.MC_renewal)
    {
      CONNECT.set_initial_condition();
      for(MKL_LONG i=0; i<4; i++)
        VAR.cnt_arr[i] = 0;
    }

  double time_update_ASSOCIATION_MAP = 0;
  double time_update_CHAIN_SUGGESTION = 0;
  
#pragma omp parallel for default(none) shared(Np, POTs, CONNECT, R_boost) num_threads(VAR.N_THREADS_BD) reduction(+: time_update_ASSOCIATION_MAP, time_update_CHAIN_SUGGESTION) if(VAR.N_THREADS_BD > 1)
  for(MKL_LONG index_particle = 0; index_particle < Np; index_particle++)
    {
      // update suggestion map after time evolution of Langevin equation
      time_update_CHAIN_SUGGESTION +=
        CONNECT.update_CHAIN_SUGGESTION_MAP_particle(index_particle, POTs, R_boost);
      // update association map after time evolution of Langevin equation
      
      time_update_ASSOCIATION_MAP +=
        CONNECT.update_ASSOCIATION_MAP_particle(index_particle, POTs, R_boost);
    }

  VAR.time_SS_update_ASSOCIATION_MAP += time_update_ASSOCIATION_MAP;
  VAR.time_SS_update_CHAIN_SUGGESTION += time_update_CHAIN_SUGGESTION;
  return dsecnd() - time_st;
}


double OMP_SS_topological_time_evolution(const MKL_LONG time_step_LV, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, POTENTIAL_SET& POTs, RDIST& R_boost, RNG_BOOST& RNG, INDEX_MC* IDX_ARR, RECORD_DATA& DATA, COND& given_condition, LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR)
// note that the TRJECTORY class dependency is deleted because of RDIST class.
{
  double time_st = dsecnd();
  if(time_step_LV%VAR.N_steps_block == 0) // the equilibration functionality is disabled at this moment. (it will be seperated for future works)
    {
      // this is rearranged in order to use the previously updated information when MC_renewal is not turned on.
      // #pragma omp for
      // initialization_topological_update(CONNECT, POTs, given_condition, VAR);
      
      VAR.time_SS_update +=
        OMP_SS_update_STATISTICS(CONNECT, POTs, R_boost, VAR);

      VAR.time_SS_CORE +=
        OMP_SS_update_topology(CONNECT, POTs, R_boost, CHAIN, RNG, DATA, IDX_ARR, LOCKER, VAR);
    }
  return dsecnd() - time_st;
}




