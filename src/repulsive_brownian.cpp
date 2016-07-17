#include "repulsive_brownian.h"

using namespace REPULSIVE_BROWNIAN;

double
REPULSIVE_BROWNIAN::
OMP_compute_RDIST(TRAJECTORY& TRAJ, const MKL_LONG index_t_now,
                  RDIST& R_boost, MKL_LONG* tmp_index_vec,
                  const MKL_LONG N_THREADS_BD)
{
  double time_st_rdist = dsecnd();
  R_boost.allocate_cells_from_positions(TRAJ, index_t_now, tmp_index_vec);
  
#pragma omp parallel for default(none) if(N_THREADS_BD > 1)	\
  shared(TRAJ, index_t_now, R_boost)                        \
  num_threads(N_THREADS_BD)                                 
    
  /*
    Originally, this parallel regime is designed to use N_cells and TOKEN.
    Because the chunk size is not easily specified, it is used to parallel with index of particles instead of cell-based.
    On this regards, the cell_index array is used which will return the cell index for the subjected particle.
  */
  for(MKL_LONG index_particle=0; index_particle<TRAJ.Np; index_particle++)
    {
      GEOMETRY::
        compute_RDIST_particle(R_boost, index_particle,
                               TRAJ, index_t_now);
    } // index_particle
  // dt_rdist += dsecnd() - time_st_rdist;
  return dsecnd() - time_st_rdist;
}

double
REPULSIVE_BROWNIAN::
OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next,
                         POTENTIAL_SET& POTs, MATRIX* force_repulsion, MATRIX* force_random,
                         RDIST& R_boost, MATRIX* vec_boost_Nd_parallel,
                         RNG_BOOST& RNG,
                         const MKL_LONG N_THREADS_BD,
                         COND& given_condition, TEMPORAL_VARIABLE& VAR)
{
  double time_st = dsecnd();
  double RF_random_xx = 0., RF_random_yy = 0., RF_random_zz = 0.;
  double RF_random_xy = 0., RF_random_xz = 0., RF_random_yz = 0.;
  double RF_repulsion_xx = 0., RF_repulsion_yy = 0., RF_repulsion_zz = 0.;
  double RF_repulsion_xy = 0., RF_repulsion_xz = 0., RF_repulsion_yz = 0.;
  double energy_repulsive_potential = 0.;
  
#pragma omp parallel for default(none) if(N_THREADS_BD > 1)     \
  shared(TRAJ, index_t_now, index_t_next,                       \
         POTs,force_repulsion, force_random,                    \
         R_boost, vec_boost_Nd_parallel,                        \
         RNG, N_THREADS_BD, given_condition, VAR)               \
  num_threads(N_THREADS_BD)                                     \
  reduction(+:RF_random_xx, RF_random_yy, RF_random_zz,         \
            RF_random_xy, RF_random_xz, RF_random_yz,           \
            RF_repulsion_xx, RF_repulsion_yy, RF_repulsion_zz,	\
            RF_repulsion_xy, RF_repulsion_xz, RF_repulsion_yz,  \
	    energy_repulsive_potential)
  
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
          
      force_repulsion[i].set_value(0);
      force_random[i].set_value(0);

      // INTEGRATOR::EULER::
      //   cal_repulsion_force_R_boost(POTs, force_repulsion[i], i, R_boost);

      INTEGRATOR::EULER::
	cal_repulsion_force_R_boost_with_RF(POTs, force_repulsion[i], i, R_boost,
					    RF_repulsion_xx, RF_repulsion_yy, RF_repulsion_zz,
					    RF_repulsion_xy, RF_repulsion_xz, RF_repulsion_yz,
					    energy_repulsive_potential);
      
      INTEGRATOR::EULER::
        cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]); 
          
      for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k)
            + TRAJ.dt*(force_repulsion[i](k))
            + sqrt(TRAJ.dt)*force_random[i](k);
        }
      if(VAR.SIMPLE_SHEAR)
        TRAJ(index_t_next, i, VAR.shear_axis) += TRAJ.dt*VAR.Wi_tau_R*TRAJ(index_t_now, i, VAR.shear_grad_axis);

      RF_random_xx += TRAJ(index_t_now, i, 0)*force_random[i](0)/sqrt(TRAJ.dt);
      RF_random_yy += TRAJ(index_t_now, i, 1)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_zz += TRAJ(index_t_now, i, 2)*force_random[i](2)/sqrt(TRAJ.dt);

      RF_random_xy += TRAJ(index_t_now, i, 0)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_xz += TRAJ(index_t_now, i, 0)*force_random[i](2)/sqrt(TRAJ.dt);
      RF_random_yz += TRAJ(index_t_now, i, 1)*force_random[i](2)/sqrt(TRAJ.dt);

      
    }
  // allocationc omputed RF values into VAR
  VAR.RF_random_xx = RF_random_xx; VAR.RF_random_yy = RF_random_yy; VAR.RF_random_zz = RF_random_zz;
  VAR.RF_random_xy = RF_random_xy; VAR.RF_random_xz = RF_random_xz; VAR.RF_random_yz = RF_random_yz;

  VAR.RF_repulsion_xx = RF_repulsion_xx/2.; VAR.RF_repulsion_yy = RF_repulsion_yy/2.; VAR.RF_repulsion_zz = RF_repulsion_zz/2.;
  VAR.RF_repulsion_xy = RF_repulsion_xy/2.; VAR.RF_repulsion_xz = RF_repulsion_xz/2.; VAR.RF_repulsion_yz = RF_repulsion_yz/2.;

  VAR.energy_repulsive_potential = energy_repulsive_potential/2.;
  
  return dsecnd() - time_st;
}


MKL_LONG
REPULSIVE_BROWNIAN::
main_EQUILIBRATION(TRAJECTORY& TRAJ,
                   POTENTIAL_SET& POTs,
                   RECORD_DATA& DATA,
                   COND& given_condition)
{
  using namespace std;
  
  double time_st_simulation = dsecnd(); 

  printf("### SIMULATION: REPULSIVE BROWNIAN ###\n");
  printf("SETTING simulation environment ...");
  
  TEMPORAL_VARIABLE VAR(given_condition, TRAJ.rows);
  mkl_set_num_threads(VAR.N_THREADS_BD);
  
  RDIST R_boost(given_condition);
  RNG_BOOST RNG(given_condition);
  MATRIX energy(1, VAR.N_components_energy, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  
  printf("DONE\nSTART SIMULATION\n\n");

  VAR.time_DIST +=         // compute RDIST with cell_list advantage
    // note that even if there is shear flow implementation,
    // the time zero is not affected by implemented shear (which means the shear flow simulation cannot be inheritance from the previous shear flow implementation at this moment)
    REPULSIVE_BROWNIAN::
    OMP_compute_RDIST(TRAJ, 0,
                      R_boost, VAR.tmp_index_vec,
                      VAR.N_THREADS_BD);

  
  VAR.time_AN += // this part related with the initial analysis from the given (or generated) positions of micelle
    ANALYSIS::
    CAL_ENERGY_R_boost(POTs, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);
  
  for(MKL_LONG t = 0; t<VAR.Nt-1; t++)
    {
      MKL_LONG index_t_now = t % VAR.N_basic;
      MKL_LONG index_t_next = (t+1) % VAR.N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
      VAR.virial_initial();
      
      if(VAR.SIMPLE_SHEAR)
        {
          double time_div_tau_R = t*TRAJ.dt;
          // printf("info PBC shift: ");
          // VAR.shear_PBC_shift = (VAR.Wi_tau_R*TRAJ.box_dimension[VAR.shear_grad_axis]*time_div_tau_R);
          VAR.shear_PBC_shift = fmod(VAR.Wi_tau_R*TRAJ.box_dimension[VAR.shear_grad_axis]*time_div_tau_R, TRAJ.box_dimension[VAR.shear_axis]); // it will remap the shear_PBC_shift into modulo scheme
          R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
          // printf("S0: %3.2f, ", R_boost.map_to_central_box_image);
          MKL_LONG central_standard = (MKL_LONG)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
          R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
          // printf("Z(S0/(L/2)): %d, M0: %3.2f\n", central_standard, R_boost.map_to_central_box_image);
        }
      VAR.time_DIST +=         // compute RDIST with cell_list advantage
        REPULSIVE_BROWNIAN::
        OMP_compute_RDIST(TRAJ, index_t_now,
                          R_boost, VAR.tmp_index_vec,
                          VAR.N_THREADS_BD);

      VAR.time_LV +=           // update Langevin equation using Euler integrator
        REPULSIVE_BROWNIAN::
        OMP_time_evolution_Euler(TRAJ, index_t_now, index_t_next,
                                 POTs, VAR.force_repulsion, VAR.force_random,
                                 R_boost, VAR.vec_boost_Nd_parallel,
                                 RNG,
                                 VAR.N_THREADS_BD,
                                 given_condition, VAR); // check arguments

      if(VAR.SIMPLE_SHEAR)
        {
          VAR.time_LV +=
            GEOMETRY::
            apply_shear_boundary_condition(TRAJ, index_t_next, VAR.shear_axis, VAR.shear_grad_axis, VAR.shear_PBC_shift);
        }
      
      VAR.time_LV +=           // keep periodic box condition
        GEOMETRY::
        minimum_image_convention_loop(TRAJ, index_t_next); // applying minimum image convention for PBC
      
      if(t%VAR.N_skip_ener==0 || t%VAR.N_skip_file==0)
        {
          // VAR.time_AN += // measuring energy of system
          //   ANALYSIS::
          //   CAL_ENERGY_R_boost(POTs, energy, TRAJ(index_t_now), R_boost);
	  energy(0) = TRAJ(index_t_now);

          VAR.time_AN +=
            VAR.record_virial_into_energy_array(energy);
	  
          VAR.time_AN +=
            REPULSIVE_BROWNIAN::
            sum_virial_components(energy);
	  
          energy(4) = 0;        // information related with number of association
          energy(5) = dsecnd() - time_st_simulation; // computation time for simulation

          VAR.time_file +=
            energy.fprint_row(DATA.ener, 0);
	  
          if(t%VAR.N_skip_file==0)
            {
              VAR.time_file += // write simulation data file
                TRAJ.fprint_row(DATA.traj, index_t_now);

              VAR.simulation_time = TRAJ(index_t_now)/atof(given_condition("repulsion_coefficient").c_str());
              VAR.time_RECORDED += // print simulation information for users
                report_simulation_info(TRAJ, energy, VAR);
            }
        }
    }

  double time_simulation = dsecnd() - time_st_simulation;
  printf("Total simulation time = %6.3e, recorded time = %6.3e (record/total = %3.1f)\n", time_simulation, VAR.time_RECORDED, VAR.time_RECORDED/time_simulation);
  return 0;
}

double
REPULSIVE_BROWNIAN::
record_simulation_data(RECORD_DATA& DATA,
                       TRAJECTORY& TRAJ, const MKL_LONG index_t_now,
                       MATRIX& energy)
{
  double time_st = dsecnd();
  TRAJ.fprint_row(DATA.traj, index_t_now);
  energy.fprint_row(DATA.ener, 0);
  return dsecnd() - time_st;
}

double
REPULSIVE_BROWNIAN::
report_simulation_info(TRAJECTORY& TRAJ,
                       MATRIX& energy,
                       TEMPORAL_VARIABLE& VAR)
{
  double total_time = VAR.time_LV + VAR.time_AN + VAR.time_file + VAR.time_DIST;
  printf("##### STEPS = %ld\tTIME = %8.6e tau_B\tENERGY = %6.3e (computing time: %4.3e)\n", TRAJ.c_t, VAR.simulation_time, energy(1), energy(5));
  printf("time consuming: LV = %3.2e (%3.1f), AN = %3.2e (%3.1f), FILE = %3.2e (%3.1f), DIST = %3.2e (%3.1f) ###\n\n", VAR.time_LV, VAR.time_LV*100/total_time, VAR.time_AN, VAR.time_AN*100/total_time, VAR.time_file, VAR.time_file*100/total_time, VAR.time_DIST, VAR.time_DIST*100/total_time);
  return total_time;
}

REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::
TEMPORAL_VARIABLE(COND& given_condition, MKL_LONG given_N_basic)
  : BROWNIAN::BROWNIAN_VARIABLE(given_condition, given_N_basic)
{
  MKL_LONG N_dimension = atoi(given_condition("N_dimension").c_str());
  vec_boost_Nd_parallel = new MATRIX [Np];
  for(MKL_LONG i=0; i<Np; i++)
    {
      vec_boost_Nd_parallel[i].initial(N_dimension, 1, 0.);
    }

  force_repulsion = new MATRIX [Np];
  for(MKL_LONG i=0; i<Np; i++)
    {
      force_repulsion[i].initial(N_dimension, 1, 0.);
    }

  virial_initial();
  N_components_energy = 24;

  if(given_condition("SIMPLE_SHEAR") == "TRUE")
    {
      Wi_tau_R = atof(given_condition("Wi_tau_C").c_str());
      Wi_tau_B = Wi_tau_R/atof(given_condition("repulsion_coefficient").c_str()); 
    }
  else
    {
      Wi_tau_R = 0.;
    }
  
  INITIALIZATION = TRUE;
};

REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::
~TEMPORAL_VARIABLE()
{
  if(INITIALIZATION)
    {
      delete[] vec_boost_Nd_parallel;
      delete[] force_repulsion;
    }
}
