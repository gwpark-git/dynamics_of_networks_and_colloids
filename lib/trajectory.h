
// Trajectory class definition for first aim
// As test purpose, it include 2-dimensional case, then extend to the 3-dimensional geometry.
#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include "matrix.h"
#include "random.h"
#include <omp.h>

/* #include "read_file_condition.h" */
#include "file_IO.h"
class TRAJECTORY : public MATRIX
{

 public:
  enum {x = 0, y = 1, z = 2}; // this is for convenience

  MKL_LONG init_ref();
  enum {N_t_test = 10, N_p_test = 3};
  // as mentioned in previously, the dimensionality on this class is set to 2 for testing purpose. In order to make generalized code, it will changed and implemented into initiator of this class, later.

  string Method;
  MKL_LONG N_energy_frequency;
  
  MKL_LONG Nt;
  MKL_LONG Np;
  // geometry. Note that all the proposed scheme should be defined in dimensionless scheme.
  MKL_LONG N_dimension;
  double *box_dimension;

  // time step
  MKL_LONG c_t;
  double dt;


  // member function
  MKL_LONG Rt_i(MKL_LONG t, MKL_LONG i_bead, MKL_LONG direction)
  {
    return index(t, i_bead*2*N_dimension + 1 + direction);
  }
  
  MKL_LONG Vt_i(MKL_LONG t, MKL_LONG i_bead, MKL_LONG direction)
  {
    return index(t, i_bead*2*N_dimension + 1 + N_dimension + direction);
  }
  
  MKL_LONG traj_read(char fn[]);
  MKL_LONG traj_write(char fn[])
  {
    fprint_out(fn);
    return 0;
  }

  // initiator
  MKL_LONG initialization(MKL_LONG N_time, MKL_LONG N_particle, double given_dt);
  MKL_LONG initialization_COND(COND& given_condition);

 TRAJECTORY() ;
 TRAJECTORY(const TRAJECTORY& TRAJ); // its for copy-constructor. It must be classified for openmp variable
 TRAJECTORY(COND& given_condition, MKL_LONG N_basic);
 // operator overloading
 // simple operator for time
 double& operator()(MKL_LONG time_t);
 // simple operator for position
 double& operator()(MKL_LONG time_t, MKL_LONG bead_i, MKL_LONG dimension_k);
 // simple operator for position and velocity. i_RV = 0 for position, i_RV = 1 for velocity
 double& operator()(MKL_LONG i_RV, MKL_LONG time_t, MKL_LONG bead_i, MKL_LONG dimension_k);
 MKL_LONG read_exist_traj(const char* fn_given_traj);
 
 virtual ~TRAJECTORY()
    {
      mkl_free(box_dimension);
      // The MATRIX and REF_MATRIX objects are automatically killed-by virtual destructor option of the library. This fact is tested with packages.
    }

 // data loader
};



MKL_LONG traj_count_line(char fn[]);

namespace GENERATOR
{
  MKL_LONG random_position_generator(TRAJECTORY& TRAJ);
  /* MKL_LONG random_vector_generator(MATRIX& R_VEC_TRANS); */
  MKL_LONG random_position_generator_REF(TRAJECTORY& TRAJ, MATRIX& R_VEC_TRANS);
}



class TRAJECTORY_HDF5
{
 public:
  double ***data;
  MKL_LONG Np;
  MKL_LONG N_dimension;
  double *box_dimension;

  MKL_LONG c_t;
  double dt;
  
  bool INITIALIZATION;
  MKL_LONG N_basic;
  
  TRAJECTORY_HDF5()
    {
      std::cout << "ERR: TRAJECTORY_HDF5 class is not support basic constructor\n";
    }
  TRAJECTORY_HDF5(COND& given_condition, MKL_LONG given_N_basic)
    {
      Np = atoi(given_condition("Np").c_str());
      N_dimension = atoi(given_condition("N_dimension").c_str());
      N_basic = given_N_basic; // it report how many time step should be stored in memory
      c_t = 0;
      dt = atof(given_condition("dt").c_str());
      box_dimension = (double*)mkl_malloc(N_dimension*sizeof(double), BIT);
      for(MKL_LONG k=0; k<N_dimension; k++)
        {
          box_dimension[k] = atof(given_condition("box_dimension").c_str());
        }
      data = (double***)mkl_malloc(N_basic*sizeof(double**), BIT);
      for(MKL_LONG t=0; t<N_basic; t++)
        {
          data[t] = (double**)mkl_malloc(Np*sizeof(double*), BIT);
          for(MKL_LONG i=0; i<Np; i++)
            {
              data[t][i] = (double*)mkl_malloc(N_dimension*sizeof(double), BIT);
            }
        }
      INITIALIZATION = TRUE;
    }
  ~TRAJECTORY_HDF5()
    {
      if(INITIALIZATION)
        {
          for(MKL_LONG t=0; t<N_basic; t++)
            {
              for(MKL_LONG i=0; i<Np; i++)
                {
                  delete data[t][i];
                }
              delete data[t];
            }
          delete data;
          delete box_dimension;
          INITIALIZATION = FALSE;
        }
    }
  
};

#endif

