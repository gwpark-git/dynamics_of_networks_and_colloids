

#include "lib_handle_association.h"

MKL_LONG ACTION::CANCEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec)
{
  return 0;
}

MKL_LONG ACTION::MOV(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec)
{
  CONNECT.del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_set[CONNECT.flag_other], index_set[CONNECT.flag_new], tmp_vec));
  return 0;
}

MKL_LONG ACTION::DEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec)
{
  CONNECT.del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  return 0;
}

MKL_LONG ACTION::ADD(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_set[], MATRIX& tmp_vec)
{
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], GEOMETRY::get_minimum_distance(TRAJ, index_t_now, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], tmp_vec));
  return 0;
}

