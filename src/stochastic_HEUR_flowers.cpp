#include "stochastic_HEUR_flowers.h"

using namespace HEUR;

double
HEUR::
record_RDIST_PARTICLE(ofstream &file_DIST,
                      RDIST& R_boost)
{
  double time_st = dsecnd();

  for(MKL_LONG index_particle=0; index_particle<R_boost.Np/2; index_particle++)
    {
      MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
      for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
        {
          MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
          for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
            {
              MKL_LONG index_target = R_boost(cell_index_neighbor, p);
              double distance = R_boost.Rsca[index_particle](index_target);
              if(distance <= R_boost.cell_length && index_particle != index_target)
                file_DIST << std::scientific << distance << '\t';
            }
        }
    }
  // file_DIST << '\n';
  return dsecnd() - time_st;
}

double
HEUR::
record_RDIST_BRIDGE(ofstream &file_DIST,
                    RDIST& R_boost, ASSOCIATION& CONNECT)
{
  double time_st = dsecnd();
  for(MKL_LONG i=0; i<CONNECT.Np/2; i++)
    {
      for(MKL_LONG j=1; j<CONNECT.TOKEN[i]; j++)
        {
          double distance = R_boost.Rsca[i](CONNECT.HASH[i](j));
          file_DIST << std::scientific << distance << '\t';
        }
    }
  // file_DIST << '\n';
  return dsecnd() - time_st;
}

MKL_LONG
HEUR::
stochastic_simulation_HEUR_flowers(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, RECORD_DATA& DATA, COND& given_condition)
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

  MATRIX energy(1, VAR.N_components_energy, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]

  VAR.time_DIST +=
    REPULSIVE_BROWNIAN::OMP_compute_RDIST(TRAJ, 0, R_boost, VAR.tmp_index_vec, VAR.N_THREADS_BD);
  
  // VAR.time_AN +=
  //   ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY_R_boost(POTs, CONNECT, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);

  RNG_BOOST RNG(given_condition);
  
  INDEX_MC *IDX_ARR = new INDEX_MC [VAR.N_THREADS_SS]; // it will call default constructor

  printf("DONE\nSTART SIMULATION\n\n");
  if(VAR.STEP_SHEAR)
    {
      MKL_LONG time_init = 0;
      
      VAR.shear_PBC_shift = fmod(VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis], TRAJ.box_dimension[VAR.shear_axis]);
      // VAR.shear_PBC_shift = VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis];
      R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
      MKL_LONG central_standard = (MKL_LONG)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
      R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
      
      // VAR.shear_PBC_shift = VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis];
      // // VAR.shear_PBC_shift = VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis];
      // R_boost.map_to_central_box_image = VAR.shear_PBC_shift;
      // MKL_LONG central_standard = 0;
      // // MKL_LONG central_standard = (MKL_LONG)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
      // // R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;

      
      GEOMETRY::
        apply_step_shear(TRAJ, time_init,
                         VAR.shear_axis, VAR.shear_grad_axis,
                         VAR.gamma_0, TRAJ.box_dimension[VAR.shear_grad_axis]);
      printf("RECORD_STEP_DEFORMATION: shear_PBC_shift = %3.1f, central_standard = %d, map_to_central_box_image = %3.1f\n", VAR.shear_PBC_shift, central_standard, R_boost.map_to_central_box_image);
    }
  if(VAR.SIMPLE_SHEAR)
    {
      printf("APPLY SIMPLE_SHEAR\n");
    }
  for(MKL_LONG t = 0; t<VAR.Nt-1; t++)
    {
      MKL_LONG index_t_now = t % VAR.N_basic;
      MKL_LONG index_t_next = (t+1) % VAR.N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
      MKL_LONG cnt = 1;
      VAR.virial_initial();

      // printf("test1\n");
      if(VAR.SIMPLE_SHEAR)
        {
          double time_div_tau_R = t*TRAJ.dt;
          // check the basic shift factor. Note that relative distance in the upper and lower box in shear gradient direction becomes Wi_R*L_D*t
          VAR.shear_PBC_shift = fmod(VAR.Wi_tau_R*TRAJ.box_dimension[VAR.shear_grad_axis]*time_div_tau_R, TRAJ.box_dimension[VAR.shear_axis]);
          // the modulo scheme will works for efficiency of getting the directly connected boundary of PBC box
          R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
          
          // central_standard to identify the middle-point shift factor since something beyond middle point represent the longer upper (or lower) boxes.
          // it is of importance to identify minimum image convention
          MKL_LONG central_standard = (MKL_LONG)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
          // this is the last parts. The corrected re-shift the relative box distance based on middle point identification
          R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
          // note that map_to_central_box_image will be used when we try to measure minimum distance between two micelles
        }
      // printf("test1\n");      
      VAR.time_DIST +=
        REPULSIVE_BROWNIAN::OMP_compute_RDIST(TRAJ, index_t_now, R_boost, VAR.tmp_index_vec, VAR.N_THREADS_BD);
      // printf("test2\n");      
      VAR.time_SS +=
        OMP_SS_topological_time_evolution(t, CONNECT, CHAIN, POTs, R_boost, RNG, IDX_ARR, DATA, given_condition, LOCKER, VAR);
      // printf("test3\n");      
      VAR.time_LV +=
        OMP_time_evolution_Euler(TRAJ, index_t_now, index_t_next,
                                 CONNECT,
                                 POTs, VAR.force_repulsion, VAR.force_random, VAR.force_spring,
                                 R_boost, VAR.vec_boost_Nd_parallel,
                                 RNG,
                                 VAR.N_THREADS_BD,
                                 given_condition, VAR);
      // printf("test4\n");
      if(VAR.SIMPLE_SHEAR || VAR.STEP_SHEAR)
        {
          // when the particles beyond boundary of PBC box in shear gradient direction, the minimum image convention should be changed
          // this functionality will remap in this case
          VAR.time_LV += 
            GEOMETRY::apply_shear_boundary_condition(TRAJ, index_t_next, VAR.shear_axis, VAR.shear_grad_axis, VAR.shear_PBC_shift);
        }
      VAR.time_LV +=
        GEOMETRY::minimum_image_convention_loop(TRAJ, index_t_next); // applying minimum image convention for PBC
      // note that minimum_image_convention must be with loop when shear flow is implemented
      // this is due to the fact that the previous minimum_image_convention cannot cope the shift higher than two times of box dimension

      VAR.N_associations = cnt_add - cnt_del;

      if(t%VAR.N_skip_ener==0 || t%VAR.N_skip_file==0)
        {
          // VAR.time_AN +=
          //   ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY_R_boost(POTs, CONNECT, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);
          energy(0) = TRAJ(index_t_now);
	  
          VAR.time_AN +=
            VAR.record_virial_into_energy_array(energy);
	  
          VAR.time_AN +=
            HEUR::sum_virial_components(energy);

          
          energy(4) = (double)VAR.N_associations; // number of associations
          energy(5) = dsecnd() - time_st_simulation; // record computation time
          energy(30) = VAR.N_diff_associations; // it record number of different associations
          
          VAR.time_file +=
            energy.fprint_row(DATA.ener, 0);

          if(t%VAR.N_skip_file==0)
            {
              VAR.time_file +=
                record_simulation_data(DATA, TRAJ, CONNECT, CHAIN, index_t_now); // neeed review

              VAR.simulation_time = TRAJ(index_t_now)/(atof(given_condition("repulsion_coefficient").c_str())*atof(given_condition("Rt").c_str()));
          
              VAR.time_RECORDED +=
                report_simulation_info(TRAJ, energy, VAR);
              if(DATA.r_dist_bridge && t%VAR.N_skip_rdist==0)
                {
                  VAR.time_RECORDED +=
                    record_RDIST_BRIDGE(DATA.r_dist_bridge, R_boost, CONNECT);
                  VAR.time_RECORDED +=
                    record_RDIST_PARTICLE(DATA.r_dist_particle, R_boost);
                }
            }

        }
    }

  double time_simulation = dsecnd() - time_st_simulation;
  printf("Total simulation time = %6.3e, recorded time = %6.3e (record/total = %3.1f)\n", time_simulation, VAR.time_RECORDED, VAR.time_RECORDED/time_simulation);
  // mkl_free(IDX_ARR);
  delete[] IDX_ARR;
  return 0;
}

double
HEUR::
record_simulation_data(RECORD_DATA& DATA, TRAJECTORY& TRAJ, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, const MKL_LONG index_t_now) 
{
  double time_st = dsecnd();
  TRAJ.fprint_row(DATA.traj, index_t_now);
  
  for(MKL_LONG ip=0; ip<TRAJ.Np; ip++)
    {
      CONNECT.HASH[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.hash, 1, CONNECT.TOKEN[ip]);
      CONNECT.weight[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.weight, 1, CONNECT.TOKEN[ip]);
    }

  if(CHAIN.INITIALIZATION)
    CHAIN.write(DATA.chain);
  return dsecnd() - time_st;
}

double
HEUR::
report_simulation_info(TRAJECTORY& TRAJ, MATRIX& energy, TEMPORAL_VARIABLE_HEUR& VAR)
{
  printf(" PBC VOLUME = %3.2e\n", VAR.volume_PBC_box);
  double total_time = VAR.time_SS + VAR.time_LV + VAR.time_AN + VAR.time_file + VAR.time_DIST;
  double sum_time_LV = VAR.time_LV_init + VAR.time_LV_force + VAR.time_LV_update;
  printf("##### STEPS = %ld\tTIME = %8.6e tau_B\tENERGY = %6.3e (computing time: %4.3e)\n", TRAJ.c_t, VAR.simulation_time, energy(1), energy(5));
  printf("time consuming: LV = %3.2e (%3.1f%%), SS = %3.2e (%3.1f%%), AN = %3.2e (%3.1f%%), FILE = %3.2e (%3.1f%%), DIST = %3.2e (%3.1f%%)\n", VAR.time_LV, VAR.time_LV*100/total_time, VAR.time_SS, VAR.time_SS*100/total_time, VAR.time_AN, VAR.time_AN*100/total_time, VAR.time_file, VAR.time_file*100/total_time, VAR.time_DIST, VAR.time_DIST*100/total_time);
  printf("LV: init = %3.1f%%, force = %3.1f%% (repulsion = %3.1f%%, random = %3.1f%%, connector = %3.1f%%), update = %3.1f%%\n", 100.*VAR.time_LV_init/sum_time_LV, 100.*VAR.time_LV_force/sum_time_LV, 100.*VAR.time_LV_force_repulsion/VAR.time_LV_force, 100.*VAR.time_LV_force_random/VAR.time_LV_force, 100.*VAR.time_LV_force_connector/VAR.time_LV_force, 100.*VAR.time_LV_update/sum_time_LV);
  printf("SS: index process = %3.1f%%, LOCKING = %3.1f%%, check_dissociation = %3.1f%%, transition = %3.1f%%, update info. = %3.1f%%, pre-ASSOCIATION = %3.1f%%, pre-SUGGESTION = %3.1f%%)  \n", 100.*VAR.time_SS_index/VAR.time_SS, 100.*VAR.time_SS_LOCK/VAR.time_SS, 100.*VAR.time_SS_check/VAR.time_SS, 100.*VAR.time_SS_transition/VAR.time_SS, 100.*VAR.time_SS_update_info/VAR.time_SS, 100.*VAR.time_SS_update_ASSOCIATION_MAP/VAR.time_SS, 100.*VAR.time_SS_update_CHAIN_SUGGESTION/VAR.time_SS);
  // 	 100.*VAR.time_SS_update/VAR.time_SS, VAR.time_SS_update_ASSOCIATION_MAP, VAR.time_SS_update_ASSOCIATION_MAP/VAR.time_SS_update, VAR.time_SS_update_CHAIN_SUGGESTION, VAR.time_SS_update_CHAIN_SUGGESTION/VAR.time_SS_update);
  printf("CHECK LAST STATISTICS: N_tot_asso = %ld, NAS = %ld, fraction=%4.3f, dNAS = %ld ####\n\n", VAR.N_tot_associable_chain, VAR.N_associations, VAR.N_associations/(double)VAR.N_tot_associable_chain, VAR.N_diff_associations);
  return total_time;
}


double
HEUR::
OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next,
                         ASSOCIATION& CONNECT,
                         POTENTIAL_SET& POTs, MATRIX* force_repulsion, MATRIX* force_random, MATRIX* force_spring,
                         RDIST& R_boost, MATRIX* vec_boost_Nd_parallel,
                         RNG_BOOST& RNG,
                         const MKL_LONG N_THREADS_BD,
                         COND& given_condition, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double time_st = dsecnd();
  double time_LV_init = 0., time_LV_force = 0., time_LV_update = 0.;
  double time_LV_force_repulsion = 0., time_LV_force_random = 0., time_LV_force_connector = 0.;

  double RF_random_xx = 0., RF_random_yy = 0., RF_random_zz = 0.;
  double RF_random_xy = 0., RF_random_xz = 0., RF_random_yz = 0.;
  double RF_repulsion_xx = 0., RF_repulsion_yy = 0., RF_repulsion_zz = 0.;
  double RF_repulsion_xy = 0., RF_repulsion_xz = 0., RF_repulsion_yz = 0.;
  double RF_connector_xx = 0., RF_connector_yy = 0., RF_connector_zz = 0.;
  double RF_connector_xy = 0., RF_connector_xz = 0., RF_connector_yz = 0.;

  double energy_repulsive_potential = 0., energy_elastic_potential = 0.;

  MKL_LONG N_diff_associations = 0;
  
#pragma omp parallel for default(none) if(N_THREADS_BD > 1)             \
  shared(TRAJ, index_t_now, index_t_next,                               \
         CONNECT, POTs, force_repulsion, force_random, force_spring,	\
         R_boost, vec_boost_Nd_parallel,  RNG, VAR, given_condition)	\
  num_threads(N_THREADS_BD)                                             \
  reduction(+: time_LV_init, time_LV_force, time_LV_update, time_LV_force_repulsion, time_LV_force_random, time_LV_force_connector, \
            RF_random_xx, RF_random_yy, RF_random_zz,                   \
            RF_random_xy, RF_random_xz, RF_random_yz,                   \
            RF_repulsion_xx, RF_repulsion_yy, RF_repulsion_zz,          \
            RF_repulsion_xy, RF_repulsion_xz, RF_repulsion_yz,          \
            RF_connector_xx, RF_connector_yy, RF_connector_zz,          \
            RF_connector_xy, RF_connector_xz, RF_connector_yz,          \
            N_diff_associations, energy_repulsive_potential, energy_elastic_potential)
	
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
        INTEGRATOR::EULER_ASSOCIATION::
        cal_connector_force_boost_with_RF(POTs, CONNECT, force_spring[i], i, R_boost.Rvec, R_boost.Rsca,
                                          RF_connector_xx, RF_connector_yy, RF_connector_zz,
                                          RF_connector_xy, RF_connector_xz, RF_connector_yz,
                                          energy_elastic_potential);
	
      // INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(POTs, CONNECT, force_spring[i], i, R_boost.Rvec, R_boost.Rsca);
	
      time_LV_force_repulsion +=
        INTEGRATOR::EULER::
        cal_repulsion_force_R_boost_with_RF(POTs, force_repulsion[i], i, R_boost,
                                            RF_repulsion_xx, RF_repulsion_yy, RF_repulsion_zz,
                                            RF_repulsion_xy, RF_repulsion_xz, RF_repulsion_yz,
                                            energy_repulsive_potential);
      // 	INTEGRATOR::EULER::cal_repulsion_force_R_boost(POTs, force_repulsion[i], i, R_boost);
      // time_LV_force_random +=
      // INTEGRATOR::EULER::cal_random_force_boost_simplified(POTs, force_random[i], RNG.BOOST_BD[it]);
      INTEGRATOR::EULER::cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]);
      
      double time_st_update = dsecnd();
      time_LV_force += time_st_update - time_st_force;
      
      for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k)
            + TRAJ.dt*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k))
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

      N_diff_associations += CONNECT.TOKEN[i] - 1; // 
        
      time_LV_update += dsecnd() - time_st_update;
    }
  VAR.time_LV_init += time_LV_init;
  VAR.time_LV_force += time_LV_force;
  VAR.time_LV_update += time_LV_update;

  VAR.time_LV_force_repulsion += time_LV_force_repulsion;
  VAR.time_LV_force_random += time_LV_force_random;
  VAR.time_LV_force_connector += time_LV_force_connector;

  /*
    The following virial stress measurements are related with the time scale.
    In this simulation code, the repulsive time scale is used in order to non-dimensionalization of Langevin equation.
    Hence, the repulsive potential is computed without repulsion coefficient C and the other potential is affected by this repulsion coefficient C. 
    For instance, bridge interaction force is divided by C and the random force is divided by sqrt(C).
   */

  // allocationc omputed RF values into VAR
  VAR.RF_random_xx = RF_random_xx/VAR.volume_PBC_box; VAR.RF_random_yy = RF_random_yy/VAR.volume_PBC_box; VAR.RF_random_zz = RF_random_zz/VAR.volume_PBC_box;
  VAR.RF_random_xy = RF_random_xy/VAR.volume_PBC_box; VAR.RF_random_xz = RF_random_xz/VAR.volume_PBC_box; VAR.RF_random_yz = RF_random_yz/VAR.volume_PBC_box;

  // the following are changed in order tune the time scale in rightful way

  // VAR.RF_repulsion_xx = RF_repulsion_xx*POTs.force_variables[0]; VAR.RF_repulsion_yy = RF_repulsion_yy*POTs.force_variables[0]; VAR.RF_repulsion_zz = RF_repulsion_zz*POTs.force_variables[0];
  // VAR.RF_repulsion_xy = RF_repulsion_xy*POTs.force_variables[0]; VAR.RF_repulsion_xz = RF_repulsion_xz*POTs.force_variables[0]; VAR.RF_repulsion_yz = RF_repulsion_yz*POTs.force_variables[0];
  double duplication_divisor = 2.;
  VAR.RF_repulsion_xx = RF_repulsion_xx/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_yy = RF_repulsion_yy/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_zz = RF_repulsion_zz/(duplication_divisor*VAR.volume_PBC_box);
  VAR.RF_repulsion_xy = RF_repulsion_xy/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_xz = RF_repulsion_xz/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_yz = RF_repulsion_yz/(duplication_divisor*VAR.volume_PBC_box);

  VAR.RF_connector_xx = RF_connector_xx/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_yy = RF_connector_yy/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_zz = RF_connector_zz/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box);
  VAR.RF_connector_xy = RF_connector_xy/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_xz = RF_connector_xz/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_yz = RF_connector_yz/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box);

  VAR.energy_elastic_potential = energy_elastic_potential/(POTs.force_variables[0]*duplication_divisor);
  // the divisor C_rep (POTs.force_variables[0]) is used in order to tune time scale as tau_R

  VAR.energy_repulsive_potential = energy_repulsive_potential/duplication_divisor;

  VAR.N_diff_associations = (MKL_LONG)N_diff_associations/2; // divisor 2 for removing duplicate count
  
  return dsecnd() - time_st;
}

double
HEUR::
write_MC_LOG_if_TRUE(bool flag_MC_LOG,
                     RECORD_DATA& DATA,
                     ASSOCIATION& CONNECT,
                     const INDEX_MC& IDX,
                     const MKL_LONG cnt, const MKL_LONG* cnt_arr,
                     const double rolling_dCDF, const double rolling_dCDF_U) 
{
  double time_st = dsecnd();
  if(flag_MC_LOG)
    {
      MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
			  
      {
        DATA.MC_LOG << cnt << '\t' << IDX.beads[CONNECT.flag_itself] << '\t' << setprecision(7) << rolling_dCDF<< '\t'  << IDX.beads[CONNECT.flag_hash_other] << '\t'  << IDX.beads[CONNECT.flag_other] << '\t'  << setprecision(7) << rolling_dCDF_U<< '\t'  << IDX.beads[CONNECT.flag_hash_backtrace] << '\t'  << IDX.beads[CONNECT.flag_new] << '\t'  << CONNECT.TOKEN[IDX.beads[CONNECT.flag_itself]]<< '\t'<< CONNECT.N_CONNECTED_ENDS(IDX.beads[CONNECT.flag_itself]) << '\t' << CONNECT.weight[IDX.beads[CONNECT.flag_itself]](0) <<'\t' <<  total_bonds << '\t'  << cnt_arr[INDEX_MC::ADD]<< '\t'  << cnt_arr[INDEX_MC::MOV]<< '\t'  << cnt_arr[INDEX_MC::OPP_DEL]<< '\t'  << cnt_arr[INDEX_MC::CANCEL] << '\t' << cnt_arr[INDEX_MC::LOCK] << endl;
      }
    } // MC_LOG
  return dsecnd() - time_st;
}

double
HEUR::
transition_single_chain_end(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, CHAIN_HANDLE& CHAIN, RDIST& R_boost, INDEX_MC& IDX, const MKL_LONG IDENTIFIER_ACTION)
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

double
HEUR::
check_dissociation_probability(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, INDEX_MC& IDX, gsl_rng* RNG_BOOST_SS_IT, MKL_LONG& IDENTIFIER_ACTION)
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

double
HEUR::
LOCKING_PARALLEL(LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR, const INDEX_MC& IDX, MKL_LONG& IDENTIFIER_ACTION, MKL_LONG& IDENTIFIER_LOCKING)
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

double
HEUR::
release_LOCKING(LOCK& LOCKER, INDEX_MC& IDX)
{
  double time_st = dsecnd();
  for(MKL_LONG I_BEADS = 0; I_BEADS < 3; I_BEADS++)
    {
      LOCKER(IDX.beads[I_BEADS]) = FALSE;
    }
  return dsecnd() - time_st;
}

double
HEUR::
micelle_selection(ASSOCIATION& CONNECT, gsl_rng* RNG_BOOST_SS_IT, INDEX_MC& IDX, RDIST& R_boost, TEMPORAL_VARIABLE_HEUR& VAR, double& rolling_dCDF, double& rolling_dCDF_U)
{
  double time_st = dsecnd();
  
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
  
  return dsecnd() - time_st;
}

double
HEUR::
OMP_SS_update_topology(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, CHAIN_HANDLE& CHAIN, RNG_BOOST& RNG, RECORD_DATA& DATA, INDEX_MC* IDX_ARR, LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double time_st = dsecnd();
  double time_SS_index = 0., time_SS_LOCK = 0., time_SS_check = 0., time_SS_transition = 0., time_SS_update_info = 0.;
#pragma omp parallel for default(none) if(VAR.N_THREADS_SS > 1)         \
  shared(DATA, POTs, CONNECT, CHAIN, LOCKER, IDX_ARR, VAR, R_boost, RNG) \
  num_threads(VAR.N_THREADS_SS)                                         \
  reduction(+: time_SS_index, time_SS_LOCK, time_SS_check, time_SS_transition, time_SS_update_info)    
  
  
  // reduction(+:dt_1, dt_2, dt_3, dt_4, dt_5, dt_6, dt_7)
  for(MKL_LONG tp=0; tp<VAR.N_tot_associable_chain; tp++)
    {
      MKL_LONG IDENTIFIER_ACTION = TRUE; // it can be 1 (IDX.ADD) but just true value
      MKL_LONG IDENTIFIER_LOCKING = FALSE;
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
        HEUR::micelle_selection(CONNECT, RNG.BOOST_SS[it], IDX_ARR[it], R_boost, VAR, rolling_dCDF, rolling_dCDF_U);
      if(CONNECT.TOKEN[IDX_ARR[it].beads[CONNECT.flag_itself]] == 0) // if there is no chain end attached for selected bead
        IDENTIFIER_ACTION = FALSE;                                   // there are not action take
      
          
      // locking
      if(VAR.N_THREADS_SS > 1)
        {
#pragma omp critical (LOCKING) // LOCKING is the name for this critical blocks
          {
            time_SS_LOCK +=      // this is differ from time_SS_LOCK since it is inside critical directive
              HEUR::LOCKING_PARALLEL(LOCKER, VAR, IDX_ARR[it], IDENTIFIER_ACTION, IDENTIFIER_LOCKING);
          }
        }
      // Note that the critical region only applicable with single thread while the others will be used in parallel regime.
      // In addition, the gap for passing the critical region will tune further gaps, then the computation speed for passing critical region will not be real critical issue.
      if(!IDENTIFIER_LOCKING) 
        {
          // This block only compute when the thread is NOT LOCKED
          time_SS_check +=
            HEUR::check_dissociation_probability(CONNECT, POTs, R_boost, IDX_ARR[it], RNG.BOOST_SS[it], IDENTIFIER_ACTION);

          time_SS_transition +=
            HEUR::transition_single_chain_end(CONNECT, POTs, CHAIN, R_boost, IDX_ARR[it], IDENTIFIER_ACTION);

          time_SS_update_info +=
            ACTION::UPDATE_INFORMATION(CONNECT, IDX_ARR[it], VAR.cnt_arr, IDENTIFIER_ACTION);

          // UNLOCKING
          // The critical directive is no more necessarly since only one thread visited each beads
          if(VAR.N_THREADS_SS > 1)
            time_SS_LOCK +=
              HEUR::release_LOCKING(LOCKER, IDX_ARR[it]);
          
#pragma omp critical (COUNTING)
          {
            /*
              critical(COUNTING) blocks:
              This is counting the action information that will be used for the future.
              Notice that the writing MC_LOG file is inside of this COUNTING critical directive, since all the information should be the same for writing (temporal)
            */
            
            VAR.cnt_arr[IDENTIFIER_ACTION] ++;

            HEUR::write_MC_LOG_if_TRUE(VAR.MC_LOG, DATA, CONNECT, IDX_ARR[it], VAR.cnt_SS + 1, VAR.cnt_arr, rolling_dCDF, rolling_dCDF_U);
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

double
HEUR::
OMP_SS_update_STATISTICS(ASSOCIATION& CONNECT,
                         POTENTIAL_SET& POTs,
                         RDIST& R_boost,
                         TEMPORAL_VARIABLE_HEUR& VAR)
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
  
#pragma omp parallel for default(none) if(VAR.N_THREADS_BD > 1)         \
  shared(Np, POTs, CONNECT, R_boost)                                    \
  num_threads(VAR.N_THREADS_BD)                                         \
  reduction(+: time_update_ASSOCIATION_MAP, time_update_CHAIN_SUGGESTION)	
  
  
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


double
HEUR::
OMP_SS_topological_time_evolution(const MKL_LONG time_step_LV,
                                  ASSOCIATION& CONNECT,
                                  CHAIN_HANDLE& CHAIN,
                                  POTENTIAL_SET& POTs,
                                  RDIST& R_boost,
                                  RNG_BOOST& RNG,
                                  INDEX_MC* IDX_ARR,
                                  RECORD_DATA& DATA,
                                  COND& given_condition,
                                  LOCK& LOCKER,
                                  TEMPORAL_VARIABLE_HEUR& VAR)
// note that the TRJECTORY class dependency is deleted because of RDIST class.
{
  double time_st = dsecnd();
  if(time_step_LV%VAR.N_steps_block == 0) // the equilibration functionality is disabled at this moment. (it will be seperated for future works)
    {
      // this is rearranged in order to use the previously updated information when MC_renewal is not turned on.
      // #pragma omp for
      // initialization_topological_update(CONNECT, POTs, given_condition, VAR);
      
      VAR.time_SS_update +=
        HEUR::OMP_SS_update_STATISTICS(CONNECT, POTs, R_boost, VAR);

      VAR.time_SS_CORE +=
        HEUR::OMP_SS_update_topology(CONNECT, POTs, R_boost, CHAIN, RNG, DATA, IDX_ARR, LOCKER, VAR);
    }
  return dsecnd() - time_st;
}

HEUR::TEMPORAL_VARIABLE_HEUR::
TEMPORAL_VARIABLE_HEUR(COND& given_condition, MKL_LONG given_N_basic)
  : REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE(given_condition, given_N_basic)
{

  MKL_LONG N_dimension = atoi(given_condition("N_dimension").c_str());
  N_steps_block = atol(given_condition("N_steps_block").c_str());
  N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());

  time_SS = 0.;
  time_SS_CORE = 0.;
  time_SS_index = 0.;
  time_SS_LOCK = 0.;
  time_SS_check = 0.;
  time_SS_transition = 0.;
  time_SS_update_info = 0.;

  time_SS_update = 0.;
  time_SS_update_ASSOCIATION_MAP = 0.;
  time_SS_update_CHAIN_SUGGESTION = 0.;

  time_LV_force_repulsion = 0.;
  time_LV_force_random = 0.;
  time_LV_force_connector = 0.;
      
  for(MKL_LONG i=0; i<6; i++)
    cnt_arr[i] = 0;

  cnt_SS = 0;
      
  force_spring = new MATRIX [Np];
  for(MKL_LONG i=0; i<Np; i++)
    force_spring[i].initial(N_dimension, 1, 0.);

  N_tot_associable_chain = Np*atoi(given_condition("N_chains_per_particle").c_str());

  MC_renewal = FALSE;
  if(given_condition("MC_renewal")=="TRUE")
    MC_renewal = TRUE;

  MC_LOG = FALSE;
  if(given_condition("MC_LOG")=="TRUE")
    MC_LOG = TRUE;

  virial_initial();
  N_components_energy = 31; // last 30 index for N_diff_associations
      
}


