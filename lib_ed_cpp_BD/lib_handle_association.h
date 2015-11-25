
#ifndef LIB_HANDLE_ASSOCIATION_H
#define LIB_HANDLE_ASSOCIATION_H
#include "lib_association.h"
#include "lib_geometry.h"
#include "lib_potential.h"
#include "lib_traj.h"

// this library is designed to handle association in easiler way.
// the reason to seperate it from lib_association is for removing dependencies of library

class INDEX_MC
{
  // index set for MC steps
 public:
  /* MKL_LONG index_set[4]; */
  MKL_LONG beads[4];
  // note that following names are omitted with from itself
  MKL_LONG &itself, &attached_bead, &new_attached_bead, &hash_attached_bead;
  /* MKL_LONG &itself = beads[2]; */
  /* MKL_LONG &attached_bead = beads[0]; */
  /* MKL_LONG &new_attached_bead = beads[1]; */
  /* MKL_LONG &hash_attached_bead = beads[3]; */

  // action_arry setting for boosting up the boolean identifier
  MKL_LONG (*ACTION_ARR[4])(TRAJECTORY&, MKL_LONG, POTENTIAL_SET&, ASSOCIATION&, MKL_LONG[], MATRIX&);
  
  /*
    The static const wokring on this way:
    const: will work as natural const
    static: will be shared with all the object from this class
   */
  static const MKL_LONG CANCEL = 0;
  static const MKL_LONG ADD = 1;
  static const MKL_LONG DEL = 2;
  static const MKL_LONG MOV = 3;
  static const MKL_LONG N_BOOST_COUNT[]; // defined on sourcefile
 INDEX_MC();
  ~INDEX_MC()
    {
      /* mkl_free(beads); */
    }
};


namespace SEARCHING
{
  MKL_LONG backtrace(MATRIX& given_arr, double p);
  MKL_LONG bisection(MATRIX& given_arr, double p);
}

namespace ACTION
{
  /* MKL_LONG const IDX_CANCEL = 0; */
  /* MKL_LONG const IDX_ADD = 1; */
  /* MKL_LONG const IDX_DEL = 2; */
  /* MKL_LONG const IDX_MOV = 3; // given by lib_association */
  /* MKL_LONG N_BOOST_COUNT[] = {0, 2, 2, 3}; */
  /* N_BOOST_COUNT[IDX_CANCEL] = 0; */
  /* N_BOOST_COUNT[IDX_ADD] = 2; */
  /* N_BOOST_COUNT[IDX_DEL] = 2; */
  /* N_BOOST_COUNT[IDX_MOV] = 3; */
  MKL_LONG ACT(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX& tmp_vec, MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, INDEX_MC& IDX);  
  MKL_LONG UPDATE_INFORMATION(ASSOCIATION& CONNECT, INDEX_MC& IDX, MKL_LONG cnt_arr[], MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG CANCEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec);
  MKL_LONG MOV(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec);
  MKL_LONG DEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec);
  MKL_LONG ADD(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec);
}


#endif
