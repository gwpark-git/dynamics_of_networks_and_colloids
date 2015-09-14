
#ifndef LIB_CONNECTIVITY_H
#define LIB_CONNECTIVITY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <mkl.h>
#include "matrix_long_ed.h"
#include "read_file_condition.h"

class CONNECTIVITY
{
 public:
  MATRIX_LONG HASH;
  MATRIX_LONG TOKEN;
  
  MKL_LONG Np; // number of particles (particle == bead == micelle)
  MKL_LONG Mc; // maximally avaliable connections
  
  CONNECTIVITY()
    {
      std::cout << "ERR:: Class CONNECTIVITY must have argument\n";
    }
  CONNECTIVITY(MKL_LONG number_of_particles, MKL_LONG maximum_connections);
  CONNECTIVITY(COND& given_condition);
  virtual ~CONNECTIVITY(){}

  MKL_LONG read_exist_hash(const char* fn_hash);
  
  MKL_LONG check_valid(MKL_LONG index_particle, MKL_LONG index_target);
  /* MKL_LONG& operator()(MKL_LONG index_particle); // return TOKEN for particle */
  /* MKL_LONG& operator()(MKL_LONG index_particle, MKL_LONG index_target); // return target index */
};

  /* MKL_LONG get_connected_bead(TRAJECTORY& TRAJ, MKL_LONG given_index); */

#endif
