#include <iostream>
#include "lib_ed_cpp_BD/matrix_ed.h"
#include "lib_ed_cpp_BD/lib_traj.h"
#include "lib_ed_cpp_BD/lib_evolution.h"
#include "lib_ed_cpp_BD/lib_association.h"
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
MKL_LONG main_NAPLE_ASSOCIATION_TEST(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, COND& given_condition);

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
      // MATRIX energy_info(Np*((MKL_LONG)Np/2 + 1), 2*TRAJ.dimension + 2, 0.);
      GENERATOR::random_position_generator(TRAJ);

      POTENTIAL_SET POTs;
      if(given_condition("Method") == "NAPLE_REPULSION")
        {
          FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_set(POTs, given_condition);
          main_NAPLE_REPULSION(TRAJ, POTs, given_condition);
        }
      else if(given_condition("Method") == "NAPLE_ASSOCIATION_TEST")
        {
          MKL_LONG Np = atol(given_condition("Np").c_str());
          MKL_LONG N_chains_per_particle = atol(given_condition("N_chains_per_particle").c_str());
          MKL_LONG TOL_connection = atol(given_condition("tolerance_allowing_connections").c_str());
          bool allowing_multiple_connections = FALSE;
          if (given_condition("allowing_multiple_connections") == "TRUE")
            allowing_multiple_connections = TRUE;
          ASSOCIATION CONNECT(Np, N_chains_per_particle, TOL_connection, allowing_multiple_connections);
          FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition);
          main_NAPLE_ASSOCIATION_TEST(TRAJ, POTs, CONNECT, given_condition);
        }
      else if(given_condition("Method") == "NAPLE_ASSOCIATION")
        {
          MKL_LONG Np = atol(given_condition("Np").c_str());
          MKL_LONG N_chains_per_particle = atol(given_condition("N_chains_per_particle").c_str());
          MKL_LONG TOL_connection = atol(given_condition("tolerance_allowing_connections").c_str());
          bool allowing_multiple_connections = FALSE;
          if (given_condition("allowing_multiple_connections") == "TRUE")
            allowing_multiple_connections = TRUE;
          ASSOCIATION CONNECT(Np, N_chains_per_particle, TOL_connection, allowing_multiple_connections);
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
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_trajectory")).c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_energy")).c_str();
  string filename_energy_info  = (given_condition("output_path") + '/' + given_condition("filename_energy_info")).c_str();
  string filename_HASH = (given_condition("output_path") + '/' + given_condition("filename_HASH")).c_str();
  string filename_CASE = (given_condition("output_path") + '/' + given_condition("filename_CASE")).c_str();
  string filename_weight = (given_condition("output_path") + '/' + given_condition("filename_weight")).c_str();
  string filename_MC_LOG = (given_condition("output_path") + '/' + given_condition("filename_MC_LOG")).c_str();
  ofstream FILE_LOG;

  if (given_condition("MC_LOG") == "TRUE")
    {
      FILE_LOG.open(filename_MC_LOG, std::ios_base::out);
      FILE_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dPDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << endl;
    }
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 4, 0.);
  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  TRAJ.fprint_row(filename_trajectory.c_str(), 0);

  MKL_LONG N_basic = TRAJ.rows;

  double bMSE_tolerance = atof(given_condition("tolerance_MSE_association").c_str());
  MKL_LONG N_max_steps = atol(given_condition("N_max_steps").c_str());

  MATRIX dPDF_U(TRAJ.Np, 1, 0.);
  MATRIX INDEX_dPDF_U(TRAJ.Np, 1, 0);

  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);

  for(MKL_LONG t=0; t<Nt-1; t++)
    {
      tmp_vec.set_value(0.);

      MKL_LONG index_t_now = t % N_basic;
      MKL_LONG index_t_next = (t+1) % N_basic;
      TRAJ(index_t_next) = (++TRAJ.c_t) * TRAJ.dt;

      // TRAJ(index_t_next)++;
      
      // while(fabs(block_MSE_next - block_MSE_now) > bMSE_tolerance || ++cnt < N_max_steps)
      MKL_LONG IDENTIFIER_MSE = TRUE;
      double block_MSE_now = 0., block_MSE_next = 10000.;
      double block_MSE_square_mean = 0., block_MSE_mean_square = 0.;
      MKL_LONG cnt = 1;

      dPDF_U.set_value(0.);
      INDEX_dPDF_U.set_value(0.);
      MKL_LONG cnt_del = 0;
      MKL_LONG cnt_add = 0;
      MKL_LONG cnt_mov = 0;
      MKL_LONG cnt_cancel = 0;

      CONNECT.initial();
      for(MKL_LONG i=0; i<CONNECT.Np; i++)
        CONNECT.TOKEN(i) = 1;
      while(IDENTIFIER_MSE && ++cnt < N_max_steps)
        {
          // choice for visiting bead          
          MKL_LONG index_itself = RANDOM::return_LONG_INT_rand(TRAJ.Np);

          // choice for selected chain end
          double rolling_dCDF = RANDOM::return_double_rand_SUP1();

          // MKL_LONG HASH_INDEX = CONNECT.GET_HASH_FROM_ROLL(index_itself, rolling_P);
          MKL_LONG index_hash_selected_chain = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF);
          MKL_LONG index_other_end_of_selected_chain = CONNECT.HASH(index_itself, index_hash_selected_chain);

          // choice for behaviour of selected chain end
          double rolling_dPDF_U = RANDOM::return_double_rand_SUP1();
          // ANALYSIS::GET_ORDERED_PDF_POTENTIAL(TRAJ, index_t_now, POTs, index_itself, INDEX_dPDF_U, dPDF_U);
          ANALYSIS::GET_ORDERED_PDF_POTENTIAL(TRAJ, index_t_now, POTs, index_other_end_of_selected_chain, INDEX_dPDF_U, dPDF_U);

          MKL_LONG k=CONNECT.Np-1;
          for(k=CONNECT.Np-1; k>=0; k--)
            {
              if(dPDF_U(k) < rolling_dPDF_U)
                {
                  k++;
                  break;
                }
            }
          MKL_LONG index_new_end_of_selected_chain = INDEX_dPDF_U(k);
          // for(MKL_LONG i=CONNECT.Np-10; i<CONNECT.Np; i++)
          //   {
          //     printf("(i, I(i), dPDF(i)) = (%ld, %ld, %6.3f), (k, roll) = (%ld, %6.3f)\n", i, (MKL_LONG) INDEX_dPDF_U(i), dPDF_U(i), k, rolling_dPDF_U);
          //   }

          // printf("(cnt, i_V, i_H, i_T) = (%ld, %ld, %ld, %ld),\t (k, i_k, N_c, TOKEN, CASE) = (%ld, %ld, %ld, %ld, %6.3e)\t(TOL=%lf)\n", cnt, index_itself, index_hash_selected_chain, index_other_end_of_selected_chain, k, index_new_end_of_selected_chain, (MKL_LONG)CONNECT.CASE(index_itself, 0), CONNECT.TOKEN(index_itself), CONNECT.CASE(index_itself, CONNECT.TOKEN(index_itself)-1), fabs(block_MSE_next-block_MSE_now));
          
          if(CONNECT.DEL_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
            {
              CONNECT.del_association_hash(index_itself, index_hash_selected_chain);
              cnt_del ++;
            }
          else if (CONNECT.NEW_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
            {
              CONNECT.add_association_INFO(POTs, index_itself, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_itself, index_new_end_of_selected_chain, tmp_vec));
              cnt_add ++;
            }
          else if (CONNECT.MOV_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
            {
              CONNECT.del_association_hash(index_itself, index_hash_selected_chain);
              CONNECT.add_association_INFO(POTs, index_itself, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_itself, index_new_end_of_selected_chain, tmp_vec));
              cnt_mov ++;
            }
          else
            {
              cnt_cancel ++;
            }
          // MKL_LONG k=0;
          MKL_LONG index_set[3] = {index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain};
          // for(MKL_LONG k=0; k<CONNECT.Np; k++)
          for(MKL_LONG i=0; i<3; i++)
            {
              CONNECTIVITY_update_Z_particle(CONNECT, index_set[i]);
              CONNECTIVITY_update_dPDF_particle(CONNECT, index_set[i]);
              CONNECTIVITY_update_dCDF_particle(CONNECT, index_set[i]);
            }
          MKL_LONG N_associations = CONNECT.N_TOTAL_ASSOCIATION()/2.;
          block_MSE_square_mean += pow(N_associations,2.);
          block_MSE_mean_square += N_associations/2.;
          if(cnt%1000 == 0 && cnt != 0)
            {
              block_MSE_now = block_MSE_next;
              block_MSE_next = (block_MSE_square_mean - pow(block_MSE_mean_square,2.))/1000.;
              if (fabs(block_MSE_next - block_MSE_now)/1000. < bMSE_tolerance)
                {
                  // printf("t, cnt, block_MSE_now, block_MSE_next, iden: %ld, %ld, %lf, %lf, %lf\n", t, cnt, block_MSE_now, block_MSE_next, fabs(block_MSE_next - block_MSE_now)/1000.);
                  // printf("\tN_asso, cnt_add, cnt_mov, cnt_del, cnt_cancel = %ld, %ld, %ld, %ld, %ld\n", CONNECT.N_TOTAL_ASSOCIATION()/2, cnt_add, cnt_mov, cnt_del, cnt_cancel);
              
                  IDENTIFIER_MSE = FALSE;
                }
              block_MSE_square_mean = 0.;
              block_MSE_mean_square = 0.;
            }
          // if(cnt%N_energy_frequency==0)
          //   {
          //     CONNECT.HASH.fprint_row(filename_HASH.c_str(), index_itself);
          //     CONNECT.CASE.fprint_row(filename_CASE.c_str(), index_itself);
          //     CONNECT.weight.fprint_row(filename_weight.c_str(), index_itself);
          if (given_condition("MC_LOG") == "TRUE")
            {
              MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
              FILE_LOG << cnt << '\t' << index_itself << '\t' << rolling_dCDF<< '\t'  << index_hash_selected_chain<< '\t'  << index_other_end_of_selected_chain<< '\t'  << rolling_dPDF_U<< '\t'  << k<< '\t'  << index_new_end_of_selected_chain<< '\t'  << CONNECT.TOKEN(index_itself)<< '\t'<<CONNECT.N_CONNECTED_ENDS(index_itself) << '\t' << CONNECT.weight(index_itself, 0) <<'\t' <<  total_bonds << '\t'  << cnt_add<< '\t'  << cnt_mov<< '\t'  << cnt_del<< '\t'  << cnt_cancel << endl;
}

          //   }
        } // while
      // IDENTIFIER_MSE = TRUE;
      // cnt = 1;
      
      // MKL_LONG index_t_next = (index_t_now + 1)%TRAJ.Nt;
      // std::cout << TRAJ(index_t_next) << std::endl;
      // TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;
      MATRIX force_spring(TRAJ.dimension, 1, 0.);
      MATRIX force_repulsion(TRAJ.dimension, 1, 0.);
      MATRIX force_random(TRAJ.dimension, 1, 0.);

#pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next) firstprivate(force_spring, force_repulsion, force_random) // firstprivate called copy-constructor while private called default constructor
      for (MKL_LONG i=0; i<TRAJ.Np; i++)
        {

          // from its print functionality, MATRIX objects are cleary working properly.
          // The reason is that with parallization, all the matrix object showed different address
          INTEGRATOR::EULER_ASSOCIATION::cal_connector_force(TRAJ, POTs, CONNECT, force_spring, index_t_now, i);
          // // the followings are for testing help conditional
          // printf("index_particle = %ld, (Fx, Fy) = (%6.3f, %6.3f), |F| = %6.3f\n", i, force_spring(0), force_spring(1), force_spring.norm());
          // if(isnan(force_spring.norm()))
          //   {
          //     printf("HASH: ");
          //     for(MKL_LONG j=0; j<CONNECT.TOKEN(i); j++)
          //       {
          //         printf("%ld, ", CONNECT.HASH(i,j));
          //       }
          //     printf("\nweight: ");
          //     for(MKL_LONG j=0; j<CONNECT.TOKEN(i); j++)
          //       {
          //         printf("%ld, ", CONNECT.weight(i,j));
          //       }
          //     printf("\n");
          //     printf("positions\n");
          //     for(MKL_LONG j=0; j<CONNECT.TOKEN(i); j++)
          //       {
          //         printf("(x, y)[P_%ld] = (%6.3f, %6.3f)\n", CONNECT.HASH(i,j), TRAJ(index_t_now, CONNECT.HASH(i,j), 0), TRAJ(index_t_now, CONNECT.HASH(i,j), 1));
          //       }
          //   }
          // // cout << i << '\t' << force_spring(0) << '\t' << force_springforce_spring.norm() << endl;
          
          INTEGRATOR::EULER::cal_repulsion_force(TRAJ, POTs, force_repulsion, index_t_now, i);
          INTEGRATOR::EULER::cal_random_force(TRAJ, POTs, force_random, index_t_now);
          for (MKL_LONG k=0; k<TRAJ.dimension; k++)
            {
              // printf("fc = %6.3e, fr = %6.3e, fR = %6.3e\n", force_spring(k), force_repulsion(k), force_random(k));
              TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*((1./POTs.force_variables[0])*force_spring(k) + force_repulsion(k)) + sqrt(TRAJ.dt)*force_random(k);
            }
        }


      // for(MKL_LONG i=0; i<CONNECT.Np; i++)
      //   {
      //     for(MKL_LONG j=0; j<CONNECT.TOKEN(i); j++)
      //       {
      //         MKL_LONG index_target = CONNECT.HASH(i,j);
      //         CONNECTIVITY_update_CASE_particle_target(CONNECT, POTs, i, index_target, GEOMETRY::get_minimum_distance(TRAJ, index_t_next, i, index_target, tmp_vec));
      //       }
      //     CONNECTIVITY_update_Z_particle(CONNECT, i);
      //     CONNECTIVITY_update_dPDF_particle(CONNECT, i);
      //     CONNECTIVITY_update_dCDF_particle(CONNECT, i);
      //   }

      // INTEGRATOR::EULER::simple_Euler(TRAJ, POTs, index_t_now); 
      GEOMETRY::minimum_image_convention(TRAJ, index_t_next);
      // ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, index_t_next);
      ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_next);
      if(t%N_skip==0)
        {
          printf("STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_next), energy(1));
          TRAJ.fprint_row(filename_trajectory.c_str(), index_t_next);
          energy.fprint(filename_energy.c_str());
          CONNECT.HASH.fprint(filename_HASH.c_str());
          CONNECT.CASE.fprint(filename_CASE.c_str());
          CONNECT.weight.fprint(filename_weight.c_str());
        }
      
      if(t%N_energy_frequency==0)
        {
          ANALYSIS::cal_detail_repulsion(TRAJ, POTs, filename_energy_info.c_str(), index_t_next);
        }
    }
  if (given_condition("MC_LOG") == "TRUE")
    FILE_LOG.close();
 
  return 0;
}

MKL_LONG main_NAPLE_ASSOCIATION_TEST(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, COND& given_condition)
{
  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_trajectory")).c_str();
  string filename_energy = (given_condition("output_path") + '/' + given_condition("filename_energy")).c_str();
  string filename_energy_info  = (given_condition("output_path") + '/' + given_condition("filename_energy_info")).c_str();
  string filename_HASH = (given_condition("output_path") + '/' + given_condition("filename_HASH")).c_str();
  string filename_CASE = (given_condition("output_path") + '/' + given_condition("filename_CASE")).c_str();
  string filename_weight = (given_condition("output_path") + '/' + given_condition("filename_weight")).c_str();
  string filename_MC_LOG = (given_condition("output_path") + '/' + given_condition("filename_MC_LOG")).c_str();
  ofstream FILE_LOG;
  FILE_LOG.open(filename_MC_LOG, std::ios_base::out);
  FILE_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dPDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << endl;
  MKL_LONG N_skip = atol(given_condition("N_skip").c_str());
  MKL_LONG N_energy_frequency = atol(given_condition("N_energy_frequency").c_str()); 

  MATRIX energy(1, 4, 0.);
  ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, 0);

  MKL_LONG Nt = atol(given_condition("Nt").c_str());
  TRAJ.fprint_row(filename_trajectory.c_str(), 0);
  TRAJ.fprint_row(filename_trajectory.c_str(), 0);

  MKL_LONG N_basic = TRAJ.rows;

  double bMSE_tolerance = atof(given_condition("tolerance_MSE_association").c_str());
  MKL_LONG N_max_steps = atol(given_condition("N_max_steps").c_str());
  MKL_LONG t=0;
  MATRIX tmp_vec(TRAJ.dimension, 1, 0.);
      
  MKL_LONG index_t_now = t % N_basic;
  MKL_LONG index_t_next = (t+1) % N_basic;
  TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;

  // TRAJ(index_t_next)++;
      
  // while(fabs(block_MSE_next - block_MSE_now) > bMSE_tolerance || ++cnt < N_max_steps)
  MKL_LONG IDENTIFIER_MSE = TRUE;
  double block_MSE_now = 0., block_MSE_next = 10000.;
  double block_MSE_square_mean = 0., block_MSE_mean_square = 0.;
  MKL_LONG cnt = 1;

  MATRIX dPDF_U(TRAJ.Np, 1, 0.);
  MATRIX INDEX_dPDF_U(TRAJ.Np, 1, 0);
  MKL_LONG cnt_del = 0;
  MKL_LONG cnt_add = 0;
  MKL_LONG cnt_mov = 0;
  MKL_LONG cnt_cancel = 0;

  while(IDENTIFIER_MSE && ++cnt < N_max_steps)
    {
      // choice for visiting bead          
      MKL_LONG index_itself = RANDOM::return_LONG_INT_rand(TRAJ.Np);

      // choice for selected chain end
      double rolling_dCDF = RANDOM::return_double_rand_SUP1();

      // MKL_LONG HASH_INDEX = CONNECT.GET_HASH_FROM_ROLL(index_itself, rolling_P);
      MKL_LONG index_hash_selected_chain = CONNECT.GET_INDEX_HASH_FROM_ROLL(index_itself, rolling_dCDF);
      MKL_LONG index_other_end_of_selected_chain = CONNECT.HASH(index_itself, index_hash_selected_chain);

      // choice for behaviour of selected chain end
      double rolling_dPDF_U = RANDOM::return_double_rand_SUP1();
      ANALYSIS::GET_ORDERED_PDF_POTENTIAL(TRAJ, index_t_now, POTs, index_itself, INDEX_dPDF_U, dPDF_U);
      MKL_LONG k=CONNECT.Np-1;
      for(k=CONNECT.Np-1; k>=0; k--)
        {
          if(dPDF_U(k) < rolling_dPDF_U)
            {
              k++;
              break;
            }
        }
      MKL_LONG index_new_end_of_selected_chain = INDEX_dPDF_U(k);

      // printf("(cnt, i_V, i_H, i_T) = (%ld, %ld, %ld, %ld),\t (k, i_k, N_c, TOKEN, CASE) = (%ld, %ld, %ld, %ld, %6.3e)\t(TOL=%lf)\n", cnt, index_itself, index_hash_selected_chain, index_other_end_of_selected_chain, k, index_new_end_of_selected_chain, (MKL_LONG)CONNECT.CASE(index_itself, 0), CONNECT.TOKEN(index_itself), CONNECT.CASE(index_itself, CONNECT.TOKEN(index_itself)-1), fabs(block_MSE_next-block_MSE_now));
          
      if(CONNECT.DEL_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
        {
          CONNECT.del_association_hash(index_itself, index_hash_selected_chain);
          cnt_del ++;
        }
      else if (CONNECT.NEW_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
        {
          CONNECT.add_association_INFO(POTs, index_itself, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_itself, index_new_end_of_selected_chain, tmp_vec));
          cnt_add ++;
        }
      else if (CONNECT.MOV_ASSOCIATION(CONNECT, index_itself, index_other_end_of_selected_chain, index_new_end_of_selected_chain))
        {
          CONNECT.del_association_hash(index_itself, index_hash_selected_chain);
          CONNECT.add_association_INFO(POTs, index_itself, index_new_end_of_selected_chain, GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_itself, index_new_end_of_selected_chain, tmp_vec));
          cnt_mov ++;
        }
      else
        {
          cnt_cancel ++;
        }
      for(MKL_LONG k=0; k<TRAJ.Np; k++)
        {
          CONNECTIVITY_update_Z_particle(CONNECT, k);
          CONNECTIVITY_update_dPDF_particle(CONNECT, k);
          CONNECTIVITY_update_dCDF_particle(CONNECT, k);
        }
      MKL_LONG N_associations = CONNECT.N_TOTAL_ASSOCIATION()/2.;
      block_MSE_square_mean += pow(N_associations,2.);
      block_MSE_mean_square += N_associations/2.;
      if(cnt%1000 == 0)
        {
          block_MSE_now = block_MSE_next;
          block_MSE_next = (block_MSE_square_mean - pow(block_MSE_mean_square,2.))/1000.;
          if (fabs(block_MSE_next - block_MSE_now)/1000. < bMSE_tolerance)
            {
              // printf("t, cnt, block_MSE_now, block_MSE_next, iden: %ld, %ld, %lf, %lf, %lf\n", t, cnt, block_MSE_now, block_MSE_next, fabs(block_MSE_next - block_MSE_now)/1000.);
              
              IDENTIFIER_MSE = FALSE;
            }
          block_MSE_square_mean = 0.;
          block_MSE_mean_square = 0.;
        }
      MKL_LONG total_bonds = CONNECT.N_TOTAL_ASSOCIATION();

      FILE_LOG << cnt << '\t' << index_itself << '\t' << rolling_dCDF<< '\t'  << index_hash_selected_chain<< '\t'  << index_other_end_of_selected_chain<< '\t'  << rolling_dPDF_U<< '\t'  << k<< '\t'  << index_new_end_of_selected_chain<< '\t'  << CONNECT.TOKEN(index_itself)<< '\t'<<CONNECT.N_CONNECTED_ENDS(index_itself) << '\t' << CONNECT.weight(index_itself, 0) <<'\t' <<  total_bonds << '\t'  << cnt_add<< '\t'  << cnt_mov<< '\t'  << cnt_del<< '\t'  << cnt_cancel << endl;

      if(!IDENTIFIER_MSE)
        {
          CONNECT.HASH.fprint(filename_HASH.c_str());
          CONNECT.CASE.fprint(filename_CASE.c_str());
          CONNECT.weight.fprint(filename_weight.c_str());
        }
          
    }

  // MKL_LONG index_t_next = (index_t_now + 1)%TRAJ.Nt;
  // std::cout << TRAJ(index_t_next) << std::endl;
  // TRAJ(index_t_next) = (++TRAJ.c_t)*TRAJ.dt;
  // MATRIX force_spring(TRAJ.dimension, 1, 0.);
  // MATRIX force_repulsion(TRAJ.dimension, 1, 0.);
  // MATRIX force_random(TRAJ.dimension, 1, 0.);

  // MKL_LONG i=0;
// #pragma omp parallel for default(none) shared(TRAJ, POTs, CONNECT, index_t_now, index_t_next) firstprivate(force_spring, force_repulsion, force_random) // firstprivate called copy-constructor while private called default constructor
//   for (i=0; i<TRAJ.Np; i++)
//     {

//       // from its print functionality, MATRIX objects are cleary working properly.
//       // The reason is that with parallization, all the matrix object showed different address
//       INTEGRATOR::EULER_ASSOCIATION::cal_connector_force(TRAJ, POTs, CONNECT, force_spring, index_t_now, i);
//       cout << force_spring.norm() << endl;
//       INTEGRATOR::EULER::cal_repulsion_force(TRAJ, POTs, force_repulsion, index_t_now, i);
//       INTEGRATOR::EULER::cal_random_force(TRAJ, POTs, force_random, index_t_now);
//       for (MKL_LONG k=0; k<TRAJ.dimension; k++)
//         {
//           // printf("fc = %6.3e, fr = %6.3e, fR = %6.3e\n", force_spring(k), force_repulsion(k), force_random(k));
//           TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k) + TRAJ.dt*(force_spring(k) + force_repulsion(k)) + sqrt(TRAJ.dt)*force_random(k);
//         }
//     }
      
  // INTEGRATOR::EULER::simple_Euler(TRAJ, POTs, index_t_now); 
  // GEOMETRY::minimum_image_convention(TRAJ, index_t_next);
  // ANALYSIS::CAL_ENERGY(TRAJ, POTs, energy, index_t_next);
  // ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY(TRAJ, POTs, CONNECT, energy, index_t_next);
  // if(t%N_skip==0)
  //   {
  //     printf("STEPS = %ld\tTIME_WR = %8.6e\tENERGY = %6.3e\n", TRAJ.c_t, TRAJ(index_t_next), energy(1));
  //     TRAJ.fprint_row(filename_trajectory.c_str(), index_t_next);
  //     energy.fprint(filename_energy.c_str());
  //   }
      
  // if(t%N_energy_frequency==0)
  //   {
  //     ANALYSIS::cal_detail_repulsion(TRAJ, POTs, filename_energy_info.c_str(), index_t_next);
  //   }

// FILE_LOG.close();
  
  return 0;
}



// -ipo option for interprocedure optimization for all files

/*
 * Local variables:
 * compile-command: "icpc -ipo -openmp -Wall -mkl -o Brownian_simulation lib_ed_cpp_BD/lib_traj.cpp lib_ed_cpp_BD/read_file_condition.cpp lib_ed_cpp_BD/lib_evolution.cpp lib_ed_cpp_BD/matrix_ed.cpp lib_ed_cpp_BD/matrix_long_ed.cpp lib_ed_cpp_BD/lib_potential.cpp lib_ed_cpp_BD/lib_connectivity.cpp lib_ed_cpp_BD/lib_association.cpp lib_ed_cpp_BD/lib_geometry.cpp lib_ed_cpp_BD/lib_random.cpp Brownian_simulation.cpp -lgsl -lm"
 * End:
 */



