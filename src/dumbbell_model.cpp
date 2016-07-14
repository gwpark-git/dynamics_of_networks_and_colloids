#include "dumbbell_model.h"

using namespace DUMBBELL;

double
DUMBBELL::
OMP_compute_RDIST(TRAJECTORY& TRAJ, const MKL_LONG index_t_now,
                  RDIST& R_boost, MKL_LONG* tmp_index_vec,
                  CONNECTIVITY& CONNECT,
                  const MKL_LONG N_THREADS_BD)
{
  double time_st_rdist = dsecnd();

#pragma omp parallel for default(none) if(N_THREADS_BD > 1) \
  shared(TRAJ, index_t_now, R_boost, CONNECT) \
  num_threads(N_THREADS_BD)

  for(MKL_LONG index_particle = 0; index_particle < TRAJ.Np; index_particle++)
    {
      // itself
      MKL_LONG index_target = CONNECT.HASH[index_particle](1); // note that there is the only one connectivity
      // note that BEYOND_BOX will not be used since cell_list will be off for dumbbell model
      double distance = R_boost.measure_minimum_distance(R_boost, TRAJ, index_t_now, index_particle, index_target, R_boost.BEYOND_BOX[0][0]);
      R_boost.Rsca[index_particle](index_target) = distance;
    }
  
  return dsecnd() - time_st_rdist;
}

MKL_LONG
DUMBBELL::
generate_dumbbell_connectivity
(CONNECTIVITY& CONNECT)
{
  if (!(CONNECT.Np % 2 == 0))
    {
      printf("ERR: Non-even number of particle is not possible to simulate dumbbell model\n");
      return -1;
    }
  MKL_LONG half_Np = CONNECT.Np/2;
  MKL_LONG hash_index_target = 1; // dumbbell only have 1 connectivity
  // for(MKL_LONG i=0; i<half_Np; i++)
  //   {
  //     CONNECT.HASH[i](hash_index_target) = i + half_Np; // it means index 0 connected with 0 + half_Np
  //     CONNECT.TOKEN[i] += 1;
  //   }

  // for(MKL_LONG i=half_Np; i<CONNECT.Np; i++)
  //   {
  //     CONNECT.HASH[i](hash_index_target) = i - half_Np;
  //     CONNECT.TOKEN[i] += 1;
  //   }

  for(MKL_LONG i=0; i<CONNECT.Np; i++)
    {
      CONNECT.HASH[i](hash_index_target) = (i + half_Np)%CONNECT.Np;
      CONNECT.TOKEN[i] = 2; // TOKEN is fixed

    }
  
  return 0;
}

double
DUMBBELL::
OMP_time_evolution_Euler
(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next, CONNECTIVITY& CONNECT, POTENTIAL_SET& POTs, MATRIX* force_random, MATRIX* force_spring, RNG_BOOST& RNG, RDIST& R_boost, const MKL_LONG N_THREADS_BD, COND& given_condition, DUMBBELL_VARIABLE& VAR)
{
  double time_st = dsecnd();

  double RF_random_xx = 0., RF_random_yy = 0., RF_random_zz = 0.;
  double RF_random_xy = 0., RF_random_xz = 0., RF_random_yz = 0.;
  double RF_connector_xx = 0., RF_connector_yy = 0., RF_connector_zz = 0.;
  double RF_connector_xy = 0., RF_connector_xz = 0., RF_connector_yz = 0.;

  
#pragma omp parallel for default(none) if(N_THREADS_BD > 1)    \
  shared(TRAJ,index_t_now, index_t_next,                       \
         CONNECT, POTs, force_random, force_spring,            \
         R_boost, RNG, N_THREADS_BD, given_condition, VAR)     \
  num_threads(N_THREADS_BD)                                    \
  reduction(+:RF_random_xx, RF_random_yy, RF_random_zz,		   \
            RF_random_xy, RF_random_xz, RF_random_yz,		   \
            RF_connector_xx, RF_connector_yy, RF_connector_zz, \
            RF_connector_xy, RF_connector_xz, RF_connector_yz)
  
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG it = omp_get_thread_num(); // get thread number for shared array objects
      
      force_random[i].set_value(0);
      force_spring[i].set_value(0);

      VAR.time_LV_force_connector +=
        INTEGRATOR::EULER::cal_connector_force_boost(POTs, CONNECT, force_spring[i], i, R_boost.Rvec, R_boost.Rsca);

      VAR.time_LV_force_random +=
        INTEGRATOR::EULER::cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]); 
      
      for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k)	// inheritance the current positions
            + TRAJ.dt*force_spring[i](k)			// connector contribution
            + sqrt(TRAJ.dt)*force_random[i](k);			// apply Wiener process
        }

      if(VAR.SIMPLE_SHEAR)					// apply simple shear into the time evolution
        TRAJ(index_t_next, i, VAR.shear_axis) += TRAJ.dt*VAR.Wi_tau_B*TRAJ(index_t_now, i, VAR.shear_grad_axis);

      RF_random_xx += TRAJ(index_t_now, i, 0)*force_random[i](0)/sqrt(TRAJ.dt);
      RF_random_yy += TRAJ(index_t_now, i, 1)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_zz += TRAJ(index_t_now, i, 2)*force_random[i](2)/sqrt(TRAJ.dt);

      RF_random_xy += TRAJ(index_t_now, i, 0)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_xz += TRAJ(index_t_now, i, 0)*force_random[i](2)/sqrt(TRAJ.dt);
      RF_random_yz += TRAJ(index_t_now, i, 1)*force_random[i](2)/sqrt(TRAJ.dt);

      RF_connector_xx += TRAJ(index_t_now, i, 0)*force_spring[i](0);
      RF_connector_yy += TRAJ(index_t_now, i, 1)*force_spring[i](1);
      RF_connector_zz += TRAJ(index_t_now, i, 2)*force_spring[i](2);

      RF_connector_xy += TRAJ(index_t_now, i, 0)*force_spring[i](1);
      RF_connector_xz += TRAJ(index_t_now, i, 0)*force_spring[i](2);
      RF_connector_yz += TRAJ(index_t_now, i, 1)*force_spring[i](2);

      
    }

  VAR.RF_random_xx = RF_random_xx; VAR.RF_random_yy = RF_random_yy; VAR.RF_random_zz = RF_random_zz;
  VAR.RF_random_xy = RF_random_xy; VAR.RF_random_xz = RF_random_xz; VAR.RF_random_yz = RF_random_yz;
  VAR.RF_connector_xx = RF_connector_xx; VAR.RF_connector_yy = RF_connector_yy; VAR.RF_connector_zz = RF_connector_zz;
  VAR.RF_connector_xy = RF_connector_xy; VAR.RF_connector_xz = RF_connector_xz; VAR.RF_connector_yz = RF_connector_yz;

  
  return dsecnd() - time_st;
}

MKL_LONG
DUMBBELL::
main_DUMBBELL
(TRAJECTORY& TRAJ, CONNECTIVITY& CONNECT, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition)
{
  using namespace std;
  
  double time_st_simulation = dsecnd(); 

  printf("### SIMULATION: PURE BROWNIAN MOTION ###\n");
  printf("SETTING simulation environment ...");
  
  DUMBBELL_VARIABLE VAR(given_condition, TRAJ.rows);
  RDIST R_boost(given_condition);
  mkl_set_num_threads(VAR.N_THREADS_BD);
  
  RNG_BOOST RNG(given_condition);
  MATRIX energy(1, VAR.N_components_energy, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]
  // MATRIX energy(1, 12, 0.);
  
  printf("DONE\nSTART SIMULATION\n\n");

  
  VAR.time_DIST +=         // compute RDIST with cell_list advantage
    DUMBBELL::
    OMP_compute_RDIST(TRAJ, 0, R_boost, VAR.tmp_index_vec, CONNECT, VAR.N_THREADS_BD);

  VAR.time_AN += // this part related with the initial analysis from the given (or generated) positions of micelle
    ANALYSIS::DUMBBELL::CAL_ENERGY_R_boost(POTs, CONNECT, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);
  
  for(MKL_LONG t = 0; t<VAR.Nt-1; t++)
    {
      MKL_LONG index_t_now = t % VAR.N_basic;
      MKL_LONG index_t_next = (t+1) % VAR.N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;

      VAR.virial_initial();

      if(VAR.SIMPLE_SHEAR)
        {
          double time_div_tau_B = t*TRAJ.dt; // note that TRAJ.dt == dt/tau_B.
          VAR.shear_PBC_shift = fmod(VAR.Wi_tau_B*TRAJ.box_dimension[VAR.shear_grad_axis]*time_div_tau_B, TRAJ.box_dimension[VAR.shear_axis]);
          R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
          MKL_LONG central_standard = (MKL_LONG)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
          R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
        }

      VAR.time_DIST +=
        DUMBBELL::OMP_compute_RDIST(TRAJ, index_t_now, R_boost, VAR.tmp_index_vec, CONNECT, VAR.N_THREADS_BD);
      
      VAR.time_LV +=           // update Langevin equation using Euler integrator
        DUMBBELL::OMP_time_evolution_Euler(TRAJ, index_t_now, index_t_next, CONNECT, POTs, VAR.force_random, VAR.force_spring, RNG, R_boost, VAR.N_THREADS_BD, given_condition, VAR); // check arguments

      if(VAR.SIMPLE_SHEAR)
        {
          // the modulo for float type, fmod, is applied in order to reduce potential overhead for minimum_image_convention function, since shift_factor is proportional to time.
          // note that the original one, shifted by shift_factor without modulo, is tested with loop-type minimum image convention without any changes of modulo scheme with single if-phrase.
          // to be on the safe side, the loop style will be used from now on
          // however, here the modulo scheme will be used instead of full shift factor in order to reduce overhead to apply minimum_image_convention
          VAR.time_LV +=
            GEOMETRY::apply_shear_boundary_condition(TRAJ, index_t_next, VAR.shear_axis, VAR.shear_grad_axis, VAR.shear_PBC_shift);
        }
      
      VAR.time_LV += // boundary condition through shear flow is extracted into the previous function, apply_shear_boundary_condition.
        GEOMETRY::minimum_image_convention_loop(TRAJ, index_t_next);

      
      if(t%VAR.N_skip_ener==0 || t%VAR.N_skip_file == 0)
        {
          VAR.time_AN += // measuring energy of system
            ANALYSIS::DUMBBELL::CAL_ENERGY_R_boost(POTs, CONNECT, energy, TRAJ(index_t_now), R_boost);
          VAR.time_AN +=
            VAR.record_virial_into_energy_array(energy);
          
          energy(4) = 0;        // information related with number of association
          energy(5) = dsecnd() - time_st_simulation; // computation time for simulation
          VAR.time_file += // write simulation data file
            energy.fprint_row(DATA.ener, 0);

          if (t%VAR.N_skip_file == 0)
            {
              VAR.time_file +=
                TRAJ.fprint_row(DATA.traj, index_t_now);

              for(MKL_LONG ip=0; ip<TRAJ.Np; ip++)
                {
                  CONNECT.HASH[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.hash, 1, CONNECT.TOKEN[ip]);
                }
              
              // VAR.time_file +=
              //   CONNECT.HASH.fprint_row(DATA.hash, index_t_now);
              
              VAR.simulation_time = TRAJ(index_t_now); // it follows the brownian time exactly
              VAR.time_RECORDED += // print simulation information for users
                BROWNIAN::report_simulation_info(TRAJ, energy, VAR);
            }
        }
    }

  double time_simulation = dsecnd() - time_st_simulation;
  printf("Total simulation time = %6.3e, recorded time = %6.3e (record/total = %3.1f)\n", time_simulation, VAR.time_RECORDED, VAR.time_RECORDED/time_simulation);
  return 0;
}

