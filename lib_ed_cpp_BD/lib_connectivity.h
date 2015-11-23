
#ifndef LIB_CONNECTIVITY_H
#define LIB_CONNECTIVITY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <mkl.h>
#include "matrix_ed.h"
#include "read_file_condition.h"

class CONNECTIVITY
{
 public:
  MATRIX *HASH;
  long *TOKEN;
  
  long Np; // number of particles (particle == bead == micelle)
  long Mc; // maximally avaliable connections
  
  CONNECTIVITY()
    {
      std::cout << "ERR:: Class CONNECTIVITY must have argument\n";
    }
  CONNECTIVITY(long number_of_particles, long maximum_connections);
  CONNECTIVITY(COND& given_condition);
  virtual ~CONNECTIVITY()
    {
      mkl_free(HASH);
      mkl_free(TOKEN);
    }

  long read_exist_hash(const char* fn_hash);
  
  long check_valid(long index_particle, long index_target);
  /* long& operator()(long index_particle); // return TOKEN for particle */
  /* long& operator()(long index_particle, long index_target); // return target index */
};


#endif
