#include "repulsive_brownian.h"


MKL_LONG main_EQUILIBRATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition)
{
  using namespace std;

  
  printf("STARTING_MAIN_ROOT\n"); 
  double time_st_simulation = dsecnd(); 
  MKL_LONG N_steps_block = atol(given_condition("N_steps_block").c_str());
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_base") + ".traj").c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_base") + ".ener").c_str();

  MKL_LONG N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
  printf("THREAD_SETTING: %ld ... ", N_THREADS_BD); 
  mkl_set_num_threads(N_THREADS_BD);
  double time_MC = 0.;
  double time_LV = 0.;
  double time_AN = 0.;
  double time_file = 0.;
  printf("DONE\n");
  printf("GENERATING BOOSTING VECTORS\n");

  MKL_LONG *tmp_index_vec = (MKL_LONG*) mkl_malloc(TRAJ.N_dimension*sizeof(MKL_LONG), BIT);
  MATRIX *vec_boost_Nd_parallel = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT); 

  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(TRAJ.N_dimension, 1, 0.);
    }
  RDIST R_boost(given_condition);
  
  printf("DONE\n"); 
  printf("FORCE VECTOR GENERATING ... "); 


  MATRIX *force_repulsion = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  MATRIX *force_random = (MATRIX*) mkl_malloc(TRAJ.Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      force_repulsion[i].initial(TRAJ.N_dimension, 1, 0.);
      force_random[i].initial(TRAJ.N_dimension, 1, 0.);
    }
  printf("DONE\n"); 
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
  
  printf("DONE\n");
  printf("START SIMULATION\n");
  // MKL_LONG N_associations = 0;

  double max_try_ASSOC = tolerance_association;
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
          R_boost.compute_RDIST_particle(index_particle, TRAJ, index_t_now);
        } // index_particle
      dt_rdist += dsecnd() - time_st_rdist;
      double time_st_MC = dsecnd();
      double time_end_MC = dsecnd();

#pragma omp parallel for default(none) shared(TRAJ, POTs, index_t_now, index_t_next, R_boost, vec_boost_Nd_parallel, force_repulsion, force_random, RNG, N_THREADS_BD, given_condition) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
      for (MKL_LONG i=0; i<TRAJ.Np; i++)
        {
          MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
          
          force_repulsion[i].set_value(0);
          force_random[i].set_value(0);

          INTEGRATOR::EULER::cal_repulsion_force_R_boost(POTs, force_repulsion[i], i, R_boost);
          INTEGRATOR::EULER::cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]); 
          
          for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
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
          // energy(4) = (double)N_associations;
          energy(4) = 0;
          energy(5) = dsecnd() - time_st_simulation;
          // ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_now, vec_boost_Nd_parallel[0]);
          time_end_AN = dsecnd();
          double total_dt = dt_1 + dt_2 + dt_3 + dt_4 + dt_5 + dt_6 + dt_7;
          double total_dt_pdf = dt_rdist + dt_pdf + dt_sort;
          double total_time = time_MC + time_LV + time_AN + time_file + total_dt_pdf;
          double dt_pdf_all = dt_pdf + dt_sort;

          printf("##### STEPS = %ld\tTIME = %8.6e tau_B\tENERGY = %6.3e (computing time: %4.3e)\n", TRAJ.c_t, TRAJ(index_t_now)/atof(given_condition("repulsion_coefficient").c_str()), energy(1), energy(5));
          printf("time consuming: MC = %3.2e (%3.1f), LV = %3.2e (%3.1f), AN = %3.2e (%3.1f), FILE = %3.2e (%3.1f), DIST = %3.2e (%3.1f)\n", time_MC, time_MC*100/total_time, time_LV, time_LV*100/total_time, time_AN, time_AN*100/total_time, time_file, time_file*100/total_time, total_dt_pdf, total_dt_pdf*100/total_time);
          printf("DIST: rdist = %6.3e (%3.1f), computing pdf = %6.3e (%3.1f), sorting pdf = %6.3e (%3.1f)\n\n", dt_rdist, 100.*dt_rdist/total_dt_pdf, dt_pdf, 100.*dt_pdf/total_dt_pdf, dt_sort, dt_sort*100./total_dt_pdf);
          TRAJ.fprint_row(DATA.traj, index_t_now);
          energy.fprint_row(DATA.ener, 0);
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
  mkl_free(force_repulsion);
  mkl_free(force_random);
  mkl_free(tmp_index_vec);
  return 0;
}
