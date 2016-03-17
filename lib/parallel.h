
#ifndef PARALLEL_H
#define PARALLEL_H

#define BIT 64
#define TRUE 1
#define FALSE 0

#include <iostream>
#include "handle_association.h"

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
