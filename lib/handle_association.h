
#ifndef HANDLE_ASSOCIATION_H
#define HANDLE_ASSOCIATION_H

#include "association.h"
#include "geometry.h"
#include "potential.h"
#include "trajectory.h"

// this library is designed to handle association in easiler way.
// the reason to seperate it from lib_association is for removing dependencies of library

class INDEX_MC
{
  // index set for MC steps
 public:
  // INDEX_MC::beads information include the all related bead indice.
  // The index value is related with the flag in the lib_association, which is not syncronized on this sourcecode at the moment.
  // Shortly, index_itself refers beads[2] this is because boolean optimization and handling indce easiler ways, and beads[0] is refer index_other and beads[1] is index_new. beads[3] is added in order to handle hash for index_other.
  // For better understanding, the reference variables are tagged and initialized in the constructor. For details, see the constructor in the sourcecode.

  MKL_LONG beads[4];
  /* // Commented for copy constructor validity. the references will be defined outside class object
  MKL_LONG &itself, &attached_bead, &new_attached_bead, &hash_attached_bead;
  */
  // action_arry setting for boosting up the boolean identifier
  // ACTION_ARR[IDX.CANCEL] act CANCEL. The IDX.CANCEL can be changed by IDX.ADD, IDX.DEL, IDX.MOV, respectively.
  MKL_LONG (*ACTION_ARR[4])(TRAJECTORY&, MKL_LONG, POTENTIAL_SET&, ASSOCIATION&, MKL_LONG[], MATRIX*);
  // Note that the last MATRIX* is related with R_minimum_distance_boost that has whole information
  
  /*
    The static const wokring on this way:
    const: will work as natural const
    static: will be shared with all the object from this class
   */
  static const MKL_LONG CANCEL = 0;
  static const MKL_LONG ADD = 1;
  static const MKL_LONG DEL = 2;
  static const MKL_LONG MOV = 3;
  static const MKL_LONG LOCK = 4;
  static const MKL_LONG N_BOOST_COUNT[]; // defined on sourcefile
  INDEX_MC();
  ~INDEX_MC(){}

  /* 
     This copy constructor is of importance to prevent the potential problem.
     Remind rule of three.
   */ 
  MKL_LONG copy_INDEX(const INDEX_MC& given_IDX);
  INDEX_MC& operator=(const INDEX_MC& given_IDX);
  INDEX_MC(const INDEX_MC& given_IDX); // copy constructore
  
  MKL_LONG set_initial_variables();
  MKL_LONG initial();
  
};

namespace SEARCHING
{
  //the basic backward searching method. Note that the probability map is given by sorted manner (this is due to the fact that the proability map is computed once a each time steps.) Therefore, if the cut-off scheme for FENE connector is properly working, there is no problem for the backtrace.
  MKL_LONG backtrace(MATRIX& given_arr, double p);
  MKL_LONG backtrace_cell_list(MATRIX& given_arr, MKL_LONG TOKEN, double p, MKL_LONG index_particle, RDIST& R_boost);
  
  // For general purpose, the bisection method is good aspect for approaching. Since the domain of this bisection is index space, that is defined on Integer space, the steps to reach solution is log_2(Np). In this scheme, however, some additional step is needed in order to compensate the loss of odd number division. Therefore, the overall step is determined log_2(Np) + #compensation, which approximately 11 in the Np640 cases (log2(640) is 8).
  MKL_LONG bisection(MATRIX& given_arr, double p);
}

namespace ACTION
{
  MKL_LONG ACT(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX* R_minimum_distance_boost, MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, INDEX_MC& IDX);  
  MKL_LONG UPDATE_INFORMATION(ASSOCIATION& CONNECT, INDEX_MC& IDX, MKL_LONG cnt_arr[], MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG CANCEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG MOV(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG DEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG ADD(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
}


#endif
