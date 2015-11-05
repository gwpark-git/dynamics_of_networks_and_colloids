

// Trajectory class definition for first aim
// As test purpose, it include 2-dimensional case, then extend to the 3-dimensional geometry.
#ifndef LIB_TRAJ_H
#define LIB_TRAJ_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include "matrix_ed.h"
#include "lib_random.h"
/* #include "lib_ed_cpp_BD/matrix_mem_ref_ed.h" */
/* #include "lib_evolution.h" */
#include <omp.h>

#include "read_file_condition.h"


class TRAJECTORY : public MATRIX
{
 private:
  enum {x = 0, y = 1, z = 2}; // this is for convenience

 public:

  long init_ref();
  enum {N_t_test = 10, N_p_test = 3};
  // as mentioned in previously, the dimensionality on this class is set to 2 for testing purpose. In order to make generalized code, it will changed and implemented into initiator of this class, later.

  string Method;
  long N_energy_frequency;
  
  long Nt;
  long Np;
  // geometry. Note that all the proposed scheme should be defined in dimensionless scheme.
  long dimension;
  double *box_dimension;

  // time step
  long c_t;
  double dt;


  // member function
  long Rt_i(long t, long i_bead, long direction)
  {
    return index(t, i_bead*2*dimension + 1 + direction);
  }
  
  long Vt_i(long t, long i_bead, long direction)
  {
    return index(t, i_bead*2*dimension + 1 + dimension + direction);
  }
  
  long traj_read(char fn[]);
  long traj_write(char fn[])
  {
    fprint_out(fn);
    return 0;
  }

  // initiator
  long initialization(long N_time, long N_particle, double given_dt);
  long initialization_COND(COND& given_condition);

 TRAJECTORY() ;
 TRAJECTORY(const TRAJECTORY& TRAJ); // its for copy-constructor. It must be classified for openmp variable
 TRAJECTORY(COND& given_condition, long N_basic);
 // operator overloading
 // simple operator for time
 double& operator()(long time_t);
 // simple operator for position
 double& operator()(long time_t, long bead_i, long dimension_k);
 // simple operator for position and velocity. i_RV = 0 for position, i_RV = 1 for velocity
 double& operator()(long i_RV, long time_t, long bead_i, long dimension_k);
 long read_exist_traj(const char* fn_given_traj);
 
 virtual ~TRAJECTORY()
    {
      /* if(force_variables) */
      /*   delete force_variables; */
      mkl_free(box_dimension);
      /* if(FILE_ENERGY_INFO) */
      /*   { */
      /*     FILE_ENERGY_INFO.close(); */
      /*   } */
      /* if(FILE_TRAJECTORY) */
      /*   { */
      /*     FILE_TRAJECTORY.close(); */
      /*   } */
      // The MATRIX and REF_MATRIX objects are automatically killed-by virtual destructor option of the library. This fact is tested with packages.
    }

 // data loader
 /* long CONTINUE_TRAJ(char fn[]); */
};


long traj_count_line(char fn[]);

namespace GENERATOR
{
  long random_position_generator(TRAJECTORY& TRAJ);
  /* long random_vector_generator(MATRIX& R_VEC_TRANS); */
  long random_position_generator_REF(TRAJECTORY& TRAJ, MATRIX& R_VEC_TRANS);
}


// namespace GENERATOR include test function for generating position and velocity for making handiful code for testing

#endif

