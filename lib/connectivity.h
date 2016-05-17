
#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <mkl.h>
#include "matrix.h"
/* #include "read_file_condition.h" */
#include "file_IO.h"

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

class CHAIN_NODE
{
  /* 
     This class is designed for suggesting node for connectivity information based on single chain subject.
     Basically, its data structure is given by linked-list, and the index itself is related with the array of object.
  */
 public:
  /* MKL_LONG &HEAD; */
  /* MKL_LONG &TAIL; */
  MKL_LONG HEAD, TAIL; // the design is changed according to the design interface.

  MKL_LONG& index(MKL_LONG flag_HEAD_TAIL)
    {
      /*
        It prevent the complicate interface design which contains duplicated definition.
      */
      if (flag_HEAD_TAIL == 1)
        return TAIL;
      return HEAD;
      // note that the reference variable can omit the & symbol when return its values.
    }
  CHAIN_NODE()
    {
    }
  virtual ~CHAIN_NODE(){
  }
};


#endif
