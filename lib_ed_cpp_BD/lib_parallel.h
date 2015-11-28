
#ifndef LIB_PARALLEL_H
#define LIB_PARALLEL_H

#define BIT 64
#define TRUE 1
#define FALSE 0

#include <iostream>
#include "lib_handle_association.h"

class LOCK
{
 private:
  MKL_LONG INITIALIZATION;
  
 public:
  bool *locker;
  MKL_LONG NL;

  bool& operator()(MKL_LONG const index_target);
  bool CHECKER(MKL_LONG const index_target);
  MKL_LONG RESET();
  LOCK();
  LOCK(MKL_LONG N_target);
  ~LOCK();
};


#endif
