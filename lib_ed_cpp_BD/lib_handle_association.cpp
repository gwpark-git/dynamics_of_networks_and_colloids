

#include "lib_handle_association.h"

const MKL_LONG INDEX_MC::N_BOOST_COUNT[] = {0, 2, 2, 3};
/*
  This is defined outside of class definition in the header file.
  The initialization for MKL_LONG CACEL and so on are fine with static constant definition.
  However, in the case of array, it should be defined on here.
  N_BOOST_COUNT[IDX.CANCEL] = 0;
  N_BOOST_COUNT[IDX.ADD] = 2;
  N_BOOST_COUNT[IDX.DEL] = 2;
  N_BOOST_COUNT[IDX.MOV] = 3;
 */

/*
  the reason for inheritance for reference variable is that the reference variables cannot be along without given variables. In this way, the reference variables automatically allocated with exiting variables with beads array.
*/
INDEX_MC::INDEX_MC() : itself(beads[2]), attached_bead(beads[0]), new_attached_bead(beads[1]), hash_attached_bead(beads[3])
{
  // index_set[4] = {0};
  // beads = (MKL_LONG*) mkl_malloc(4*sizeof(MKL_LONG), BIT);
  // itself = beads[2];
  // attached_bead = beads[0];
  // new_attached_bead = beads[1];
  // hash_attached_bead = beads[3];
      
  ACTION_ARR[CANCEL] = ACTION::CANCEL;
  ACTION_ARR[ADD]    = ACTION::ADD   ;
  ACTION_ARR[DEL]    = ACTION::DEL   ;
  ACTION_ARR[MOV]    = ACTION::MOV   ;
}



MKL_LONG ACTION::ACT(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX& tmp_vec, MKL_LONG const IDENTIFIER_ACTION)
{
  IDX.ACTION_ARR[IDENTIFIER_ACTION](TRAJ, index_t_now, POTs, CONNECT, IDX.beads, tmp_vec);
  return 0;
}


MKL_LONG ACTION::UPDATE_INFORMATION(ASSOCIATION& CONNECT, INDEX_MC& IDX, MKL_LONG cnt_arr[], MKL_LONG const IDENTIFIER_ACTION)
{
  cnt_arr[IDENTIFIER_ACTION]++;
  for(MKL_LONG i=0; i<IDX.N_BOOST_COUNT[IDENTIFIER_ACTION]; i++)
    {
      CONNECT.update_Z_particle(IDX.beads[i]);
      CONNECT.update_dPDF_particle(IDX.beads[i]);
      CONNECT.update_dCDF_particle(IDX.beads[i]);
    }
  return 0;
}

MKL_LONG ACTION::CANCEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec)
{
  return 0;
}

MKL_LONG ACTION::MOV(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec)
{
  CONNECT.del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_set[CONNECT.flag_other], index_set[CONNECT.flag_new], tmp_vec));
  return 0;
}

MKL_LONG ACTION::DEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec)
{
  CONNECT.del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  return 0;
}

MKL_LONG ACTION::ADD(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX& tmp_vec)
{
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], tmp_vec));
  return 0;
}


MKL_LONG ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, INDEX_MC& IDX)
{
  // MKL_LONG &index_itself = index_set[2], &index_other = index_set[0], &index_new = index_set[1];
  // check index_other == index_new
  if (CONNECT.CANCEL_ASSOCIATION(CONNECT, IDX.beads))
    {
      // this is the most probable case for canceling
      return IDX.CANCEL;
    }
  // check index_itself == index_new
  else if (CONNECT.ADD_ASSOCIATION(CONNECT, IDX.beads))
    {
      if (CONNECT.CHECK_N_ADD_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.ADD;
      else
        return IDX.CANCEL;
    }
  // check index_itself == index_other
  else if (CONNECT.DEL_ASSOCIATION(CONNECT, IDX.beads))
    {
      if (CONNECT.CHECK_N_DEL_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.DEL;
      else return IDX.CANCEL;
    }
  // with nested case, this is only the case of MOV
  else
    {
      if (CONNECT.CHECK_N_MOV_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.MOV;
      return IDX.CANCEL;
    }
  printf("ERR:TRUTH_MAP::NO SPECIFIED ACTION CASE\n");
  return IDX.CANCEL;
}

MKL_LONG SEARCHING::bisection(MATRIX& given_arr, double p)
{
  MKL_LONG N = given_arr.size;
  MKL_LONG k = N/2;
  // double dk = N/2.;
  MKL_LONG dk = N/2;

  do
    {
      // dk /= 2;
      // the conditional phrase (dk%2 != 0) will return 0 or 1 when the modulo is zero or not, respectively
      // this compensate the loss of approaching because of given domain is in integer rather than real number
      // because of it, the convergence rate is slower than the previous one (log2(Np))
      dk = dk/2 + (dk%2 != 0);
      if (given_arr(k) < p)
        k += dk;
      else if (given_arr(k) > p)
        k -= dk;
      else
        return k;
    } while (dk > 1);
    
  if (given_arr(k) < p)
    k += 1;
  return k;
}

MKL_LONG SEARCHING::backtrace(MATRIX& given_arr, double p)
{
  for(MKL_LONG k= given_arr.size - 1; k >=0; k--)
    {
      if(given_arr(k) < p)
        {
          return k+1;
        }
    }
  return 0;
}

