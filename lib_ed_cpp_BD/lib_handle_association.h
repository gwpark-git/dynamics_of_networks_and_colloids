
#ifndef LIB_HANDLE_ASSOCIATION_H
#define LIB_HANDLE_ASSOCIATION_H
#include "lib_association.h"
#include "lib_geometry.h"
#include "lib_potential.h"
#include "lib_traj.h"

// this library is designed to handle association in easiler way.
// the reason to seperate it from lib_association is for removing dependencies of library

namespace ACTION
{
  MKL_LONG const IDX_CANCEL = 0;
  MKL_LONG const IDX_ADD = 1;
  MKL_LONG const IDX_DEL = 2;
  MKL_LONG const IDX_MOV = 3; // given by lib_association
    
  MKL_LONG CANCEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec);
  MKL_LONG MOV(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec);
  MKL_LONG DEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec);
  MKL_LONG ADD(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec);
}


#endif
