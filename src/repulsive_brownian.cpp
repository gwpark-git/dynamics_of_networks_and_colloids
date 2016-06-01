#include "repulsive_brownian.h"

using namespace REPULSIVE_BROWNIAN;


double REPULSIVE_BROWNIAN::OMP_compute_RDIST(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, RDIST& R_boost, MKL_LONG* tmp_index_vec, const MKL_LONG N_THREADS_BD)
{
  double time_st_rdist = dsecnd();
  R_boost.allocate_cells_from_positions(TRAJ, index_t_now, tmp_index_vec);
  
#pragma omp parallel for default(none) shared(TRAJ, index_t_now, R_boost) num_threads(N_THREADS_BD) if(N_THREADS_BD > 1)
  /*
    Originally, this parallel regime is designed to use N_cells and TOKEN.
    Because the chunk size is not easily specified, it is used to parallel with index of particles instead of cell-based.
    On this regards, the cell_index array is used which will return the cell index for the subjected particle.
  */
  for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
    {
      R_boost.compute_RDIST_particle(index_particle, TRAJ, index_t_now);
    } // index_particle
  // dt_rdist += dsecnd() - time_st_rdist;
  return dsecnd() - time_st_rdist;
}


double REPULSIVE_BROWNIAN::OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next, POTENTIAL_SET& POTs, RDIST& R_boost, MATRIX* vec_boost_Nd_parallel, MATRIX* force_repulsion, MATRIX* force_random, RNG_BOOST& RNG, const MKL_LONG N_THREADS_BD, COND& given_condition)
{
  double time_st = dsecnd();
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
  return dsecnd() - time_st;
}


MKL_LONG REPULSIVE_BROWNIAN::main_EQUILIBRATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition)
{
  using namespace std;
  
  double time_st_simulation = dsecnd(); 

  printf("### SIMULATION: REPULSIVE BROWNIAN ###\n");
  printf("SETTING simulation environment ...");
  
  TEMPORAL_VARIABLE VAR(given_condition, TRAJ.rows);
  // initial_VAR(VAR, given_condition, TRAJ.rows);
  mkl_set_num_threads(VAR.N_THREADS_BD);
  
  RDIST R_boost(given_condition);
  RNG_BOOST RNG(given_condition);
  MATRIX energy(1, 6, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  // MATRIX energy(1, 12, 0.);
  
  printf("DONE\nSTART SIMULATION\n\n");
  
  VAR.time_AN += // this part related with the initial analysis from the given (or generated) positions of micelle
    ANALYSIS::CAL_ENERGY_R_boost(POTs, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);
  
  for(MKL_LONG t = 0; t<VAR.Nt-1; t++)
    {
      MKL_LONG index_t_now = t % VAR.N_basic;
      MKL_LONG index_t_next = (t+1) % VAR.N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;

      VAR.time_DIST +=         // compute RDIST with cell_list advantage
        REPULSIVE_BROWNIAN::OMP_compute_RDIST(TRAJ, index_t_now, R_boost, VAR.tmp_index_vec, VAR.N_THREADS_BD);

      VAR.time_LV +=           // update Langevin equation using Euler integrator
        REPULSIVE_BROWNIAN::OMP_time_evolution_Euler(TRAJ, index_t_now, index_t_next, POTs, R_boost, VAR.vec_boost_Nd_parallel, VAR.force_repulsion, VAR.force_random, RNG, VAR.N_THREADS_BD, given_condition); // check arguments

      VAR.time_LV +=           // keep periodic box condition
        GEOMETRY::minimum_image_convention(TRAJ, index_t_next); // applying minimum image convention for PBC
      
      if(t%VAR.N_skip==0)
        {
          VAR.time_AN += // measuring energy of system
            ANALYSIS::CAL_ENERGY_R_boost(POTs, energy, TRAJ(index_t_now), R_boost);
          
          energy(4) = 0;        // information related with number of association
          energy(5) = dsecnd() - time_st_simulation; // computation time for simulation

          VAR.time_file += // write simulation data file
            record_simulation_data(DATA, TRAJ, energy, index_t_now);

          VAR.simulation_time = TRAJ(index_t_now)/atof(given_condition("repulsion_coefficient").c_str());
          VAR.time_RECORDED += // print simulation information for users
            report_simulation_info(TRAJ, energy, VAR);
            // report_simulation_info(TRAJ, VAR.time_LV, VAR.time_AN, VAR.time_file, VAR.time_DIST);
        }
    }

  double time_simulation = dsecnd() - time_st_simulation;
  printf("Total simulation time = %6.3e, recorded time = %6.3e (record/total = %3.1f)\n", time_simulation, VAR.time_RECORDED, VAR.time_RECORDED/time_simulation);
  return 0;
}

double REPULSIVE_BROWNIAN::record_simulation_data(RECORD_DATA& DATA, TRAJECTORY& TRAJ, MATRIX& energy, const MKL_LONG index_t_now)
{
  double time_st = dsecnd();
  TRAJ.fprint_row(DATA.traj, index_t_now);
  energy.fprint_row(DATA.ener, 0);
  return dsecnd() - time_st;
}

double REPULSIVE_BROWNIAN::report_simulation_info(TRAJECTORY& TRAJ, MATRIX& energy, TEMPORAL_VARIABLE& VAR)
{
  double total_time = VAR.time_LV + VAR.time_AN + VAR.time_file + VAR.time_DIST;
  printf("##### STEPS = %ld\tTIME = %8.6e tau_B\tENERGY = %6.3e (computing time: %4.3e)\n", TRAJ.c_t, VAR.simulation_time, energy(1), energy(5));
  printf("time consuming: LV = %3.2e (%3.1f), AN = %3.2e (%3.1f), FILE = %3.2e (%3.1f), DIST = %3.2e (%3.1f) ###\n\n", VAR.time_LV, VAR.time_LV*100/total_time, VAR.time_AN, VAR.time_AN*100/total_time, VAR.time_file, VAR.time_file*100/total_time, VAR.time_DIST, VAR.time_DIST*100/total_time);
  return total_time;
}

REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::TEMPORAL_VARIABLE(COND& given_condition, MKL_LONG given_N_basic)
{
  Np = atoi(given_condition("Np").c_str());
  MKL_LONG N_dimension = atoi(given_condition("N_dimension").c_str());
  // MKL_LONG N_dimension = Np;
  N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
  // tmp_index_vec = (MKL_LONG*) mkl_malloc(N_dimension*sizeof(MKL_LONG), BIT);
  // vec_boost_Nd_parallel = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT); 
  tmp_index_vec = new MKL_LONG [N_dimension];
  vec_boost_Nd_parallel = new MATRIX [Np];
  for(MKL_LONG i=0; i<Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(N_dimension, 1, 0.);
    }

  // force_repulsion = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  // force_random = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  force_repulsion = new MATRIX [Np];
  force_random = new MATRIX [Np];
  for(MKL_LONG i=0; i<Np; i++)
    {
      force_repulsion[i].initial(N_dimension, 1, 0.);
      force_random[i].initial(N_dimension, 1, 0.);
    }

  N_skip = atol(given_condition("N_skip").c_str());

  Nt = atol(given_condition("Nt").c_str());
  N_basic = given_N_basic;

  time_LV = 0.; time_DIST = 0.; time_file = 0.; time_AN = 0.; time_RECORDED = 0.;
  time_LV_init = 0.; time_LV_force = 0.; time_LV_update = 0.;
  simulation_time = 0.;
  INITIALIZATION = TRUE;
};

// double destruct_VAR(TEMPORAL_VARIABLES& VAR)
REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::~TEMPORAL_VARIABLE()
{
  if(INITIALIZATION)
    {
      // mkl_free(vec_boost_Nd_parallel);
      // mkl_free(force_repulsion);
      // mkl_free(force_random);
      // mkl_free(tmp_index_vec);
      delete[] vec_boost_Nd_parallel;
      delete[] force_repulsion;
      delete[] force_random;
      delete[] tmp_index_vec;
    }
}
