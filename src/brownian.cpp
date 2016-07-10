#include "brownian.h"

using namespace BROWNIAN;

double
BROWNIAN::
OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next,
			 POTENTIAL_SET& POTs, MATRIX* force_random,
			 RNG_BOOST& RNG,
			 const MKL_LONG N_THREADS_BD,
			 COND& given_condition, BROWNIAN_VARIABLE& VAR)
{
  double time_st = dsecnd();
  double RF_random_xx = 0., RF_random_yy = 0., RF_random_zz = 0.;
  double RF_random_xy = 0., RF_random_xz = 0., RF_random_yz = 0.;
  // note that reduction in OpenMP cannot take any of member variable and array type variables
#pragma omp parallel for default(none) if(N_THREADS_BD > 1)		\
  shared(TRAJ, index_t_now, index_t_next,				\
	 POTs, force_random,						\
	 RNG, N_THREADS_BD, given_condition, VAR)			\
  num_threads(N_THREADS_BD)						\
  reduction(+:RF_random_xx, RF_random_yy, RF_random_zz, RF_random_xy, RF_random_xz, RF_random_yz)
  
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
      
      force_random[i].set_value(0);

      INTEGRATOR::EULER::cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]); 
      
      for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) // inheritance the current positions
            + sqrt(TRAJ.dt)*force_random[i](k);              // apply Wiener process
        }

      // compute for virials
      // note that weird sequencial sum instead of array index type
      // is for compatibility with reduction operator in OpenMP
      RF_random_xx += TRAJ(index_t_now, i, 0)*force_random[i](0)/sqrt(TRAJ.dt);
      RF_random_yy += TRAJ(index_t_now, i, 1)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_zz += TRAJ(index_t_now, i, 2)*force_random[i](2)/sqrt(TRAJ.dt);

      RF_random_xy += TRAJ(index_t_now, i, 0)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_xz += TRAJ(index_t_now, i, 0)*force_random[i](2)/sqrt(TRAJ.dt);
      RF_random_yz += TRAJ(index_t_now, i, 1)*force_random[i](2)/sqrt(TRAJ.dt);
      
      
      // apply simple shear into the time evolution
      if(VAR.SIMPLE_SHEAR)
        TRAJ(index_t_next, i, VAR.shear_axis) += TRAJ.dt*VAR.Wi_tau_B*TRAJ(index_t_now, i, VAR.shear_grad_axis);

    }

  // allocationc omputed RF values into VAR
  VAR.RF_random_xx = RF_random_xx; VAR.RF_random_yy = RF_random_yy; VAR.RF_random_zz = RF_random_zz;
  VAR.RF_random_xy = RF_random_xy; VAR.RF_random_xz = RF_random_xz; VAR.RF_random_yz = RF_random_yz;
  return dsecnd() - time_st;
}

MKL_LONG
BROWNIAN::
main_PURE_BROWNIAN
(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition)
{
  using namespace std;
  
  double time_st_simulation = dsecnd(); 

  printf("### SIMULATION: PURE BROWNIAN MOTION ###\n");
  printf("SETTING simulation environment ...");
  
  BROWNIAN_VARIABLE VAR(given_condition, TRAJ.rows);
  mkl_set_num_threads(VAR.N_THREADS_BD);
  
  RNG_BOOST RNG(given_condition);
  MATRIX energy(1, VAR.N_components_energy, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // from 6 to 11: total virial stress
  // from 12 to 17: virial stress come from random collision
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  // MATRIX energy(1, 12, 0.);
  
  printf("DONE\nSTART SIMULATION\n\n");

  VAR.time_AN += // this part related with the initial analysis from the given (or generated) positions of micelle
    ANALYSIS::CAL_ENERGY_BROWNIAN(POTs, energy, TRAJ(0));
  
  for(MKL_LONG t = 0; t<VAR.Nt-1; t++)
    {
      MKL_LONG index_t_now = t % VAR.N_basic;
      MKL_LONG index_t_next = (t+1) % VAR.N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;

      VAR.virial_initial();

      VAR.time_LV +=           // update Langevin equation using Euler integrator
        BROWNIAN::OMP_time_evolution_Euler(TRAJ, index_t_now, index_t_next, POTs, VAR.force_random, RNG, VAR.N_THREADS_BD, given_condition, VAR); // check arguments

      if(VAR.SIMPLE_SHEAR)
        {
          double time_div_tau_B = t*TRAJ.dt; // note that TRAJ.dt == dt/tau_B.

          VAR.shear_PBC_shift = fmod(VAR.Wi_tau_B*TRAJ.box_dimension[VAR.shear_grad_axis]*time_div_tau_B, TRAJ.box_dimension[VAR.shear_axis]);
	  
          // the modulo for float type, fmod, is applied in order to reduce potential overhead for minimum_image_convention function, since shift_factor is proportional to time.
          // note that the original one, shifted by shift_factor without modulo, is tested with loop-type minimum image convention without any changes of modulo scheme with single if-phrase.
          // to be on the safe side, the loop style will be used from now on
          // however, here the modulo scheme will be used instead of full shift factor in order to reduce overhead to apply minimum_image_convention
          
          VAR.time_LV +=
            GEOMETRY::apply_shear_boundary_condition(TRAJ, index_t_next, VAR.shear_axis, VAR.shear_grad_axis, VAR.shear_PBC_shift);
	  
        }
      
      VAR.time_LV += // boundary condition through shear flow is extracted into the previous function, apply_shear_boundary_condition.
        GEOMETRY::minimum_image_convention_loop(TRAJ, index_t_next);

      if(t%VAR.N_skip_ener==0 || t%VAR.N_skip_file==0)
        {
          VAR.time_AN += // measuring energy of system
            ANALYSIS::CAL_ENERGY_BROWNIAN(POTs, energy, TRAJ(index_t_now));

	  VAR.time_AN +=
	    VAR.record_virial_into_energy_array(energy);

	  VAR.time_AN +=
	    BROWNIAN::sum_virial_components(energy);
	  
          energy(4) = 0;        // information related with number of association
          energy(5) = dsecnd() - time_st_simulation; // computation time for simulation
	  VAR.time_file +=
	    energy.fprint_row(DATA.ener, 0);

	  if(t%VAR.N_skip_file==0)
	    {
	      // the writing file is individual from writing energy
	      VAR.time_file += // write simulation data file
		TRAJ.fprint_row(DATA.traj, index_t_now);
	      VAR.simulation_time = TRAJ(index_t_now); // the repulsion coefficient is not necessary for pure Brownian motion
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
BROWNIAN::
report_simulation_info
(TRAJECTORY& TRAJ, MATRIX& energy, BROWNIAN_VARIABLE& VAR)
{
  double total_time = VAR.time_LV + VAR.time_AN + VAR.time_file + VAR.time_DIST;
  printf("##### STEPS = %ld\tTIME = %8.6e tau_B\tENERGY = %6.3e (computing time: %4.3e)\n", TRAJ.c_t, VAR.simulation_time, energy(1), energy(5));
  printf("time consuming: LV = %3.2e (%3.1f), AN = %3.2e (%3.1f), FILE = %3.2e (%3.1f), DIST = %3.2e (%3.1f) ###\n\n", VAR.time_LV, VAR.time_LV*100/total_time, VAR.time_AN, VAR.time_AN*100/total_time, VAR.time_file, VAR.time_file*100/total_time, VAR.time_DIST, VAR.time_DIST*100/total_time);
  return total_time;
}

BROWNIAN::BROWNIAN_VARIABLE::
BROWNIAN_VARIABLE
(COND& given_condition, MKL_LONG given_N_basic)
{
  Np = atoi(given_condition("Np").c_str());
  MKL_LONG N_dimension = atoi(given_condition("N_dimension").c_str());
  volume_PBC_box = pow(atof(given_condition("box_dimension").c_str()), N_dimension); // note that the definition should be changed when the box dimension have dependency with the direction
  N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
  tmp_index_vec = new MKL_LONG [N_dimension];
  force_random = new MATRIX [Np];
  for(MKL_LONG i=0; i<Np; i++)
    {
      force_random[i].initial(N_dimension, 1, 0.);
    }

  N_skip_ener = atol(given_condition("N_skip_ener").c_str());
  N_skip_file = atol(given_condition("N_skip_file").c_str());  

  Nt = atol(given_condition("Nt").c_str());
  N_basic = given_N_basic;

  time_LV = 0.; time_DIST = 0.; time_file = 0.; time_AN = 0.; time_RECORDED = 0.;
  time_LV_init = 0.; time_LV_force = 0.; time_LV_update = 0.;
  time_LV_force_random = 0.;
  simulation_time = 0.;

  virial_initial();
  N_components_energy = 18;
  
  if(given_condition("SIMPLE_SHEAR")=="TRUE")
    {
      SIMPLE_SHEAR = TRUE;
      Wi_tau_B = atof(given_condition("Wi_tau_C").c_str()); // the tau_C in the pure Brownian come with tau_B
      shear_axis = atoi(given_condition("shear_axis").c_str()); // 0 will be set as default. (x-axis)
      shear_grad_axis = atoi(given_condition("shear_grad_axis").c_str()); // 1 will be set as default. (y-axis)
      shear_PBC_shift = 0.; // initially, it is zero
    }
  else
    {
      SIMPLE_SHEAR = FALSE;
      Wi_tau_B = 0.;
      // the following are set to be on the safe side (sequence effct)
      shear_axis = 0;
      shear_grad_axis = 1;
      shear_PBC_shift = 0;
    }

  INITIALIZATION = TRUE;
};

BROWNIAN::BROWNIAN_VARIABLE::
~BROWNIAN_VARIABLE()
{
  if(INITIALIZATION)
    {
      delete[] force_random;
      delete[] tmp_index_vec;
    }
}
