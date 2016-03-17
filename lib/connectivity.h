
#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <mkl.h>
#include "matrix.h"
#include "read_file_condition.h"

class CONNECTIVITY
{
 public:
  MATRIX *HASH;
  MKL_LONG *TOKEN;
  
  MKL_LONG Np; // number of particles (particle == bead == micelle)
  MKL_LONG Mc; // maximally avaliable connections
  
  CONNECTIVITY()
    {
      std::cout << "ERR:: Class CONNECTIVITY must have argument\n";
    }
  CONNECTIVITY(MKL_LONG number_of_particles, MKL_LONG maximum_connections);
  CONNECTIVITY(COND& given_condition);
  virtual ~CONNECTIVITY()
    {
      mkl_free(HASH);
      mkl_free(TOKEN);
    }

  MKL_LONG read_exist_hash(const char* fn_hash);
  
  MKL_LONG check_valid(MKL_LONG index_particle, MKL_LONG index_target);
};


#endif
