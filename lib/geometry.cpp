
#include "geometry.h"

MKL_LONG RDIST::compute_RDIST_particle(const MKL_LONG index_particle, TRAJECTORY& TRAJ, MKL_LONG index_t)
{
  MKL_LONG cell_index_particle = cell_index[index_particle];
  for(MKL_LONG k=0; k<N_neighbor_cells; k++)
    {
      MKL_LONG cell_index_neighbor = NEIGHBOR_CELLS[cell_index_particle][k];
      for(MKL_LONG p=0; p<TOKEN[cell_index_neighbor]; p++)
        {
          MKL_LONG index_target = (*this)(cell_index_neighbor, p);
          double distance = GEOMETRY::get_minimum_distance_cell_list(TRAJ, index_t, index_particle, index_target, Rvec[index_particle][index_target], BEYOND_BOX[cell_index_particle][k]);
          Rsca[index_particle](index_target) = distance;
        } // p
    } // k
  return 0;
}


RDIST::RDIST(COND& given_condition) : CLIST(given_condition)
{
  Rvec = (MATRIX**)mkl_malloc(Np*sizeof(MATRIX*), BIT);
  Rsca = (MATRIX*)mkl_malloc(Np*sizeof(MATRIX), BIT);
      
  for(MKL_LONG i=0; i<Np; i++)
    {
      // since the space complexity is not the matter for our simulation (at this moment),
      // the Rvec have Np*Np array that is much larger when we used cell-list approaches
      Rvec[i] = (MATRIX*)mkl_malloc(Np*sizeof(MATRIX), BIT);
      for(MKL_LONG j=0; j<Np; j++)
        {
          Rvec[i][j].initial(N_dimension, 1, 0.);
        }
      Rsca[i].initial(Np, 1, 0.);
    }
}


double GEOMETRY::get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost)
{
  // this function generating minimum relative vector into R_minimum_vec_boost,
  // and distance of it into R_minimum_boost.
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      R_minimum_boost_particle(i) = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, i, R_minimum_vec_boost[index_particle][i]);
    }
  return 0;
}

// double GEOMETRY::get_minimum_distance_for_particle(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost)
// {
//   // this function generating minimum relative vector into R_minimum_vec_boost,
//   // and distance of it into R_minimum_boost.
//   for(MKL_LONG i=0; i<TRAJ.Np; i++)
//     {
//       R_minimum_boost_particle(i) = GEOMETRY::get_minimum_distance(TRAJ, index_t, index_particle, i, R_minimum_vec_boost[index_particle][i]);
//     }
//   return 0;
// }


double GEOMETRY::get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec)
{
  // this comentation is related with the assumes that
  // the given dimensionality is properly working
  // therefore, it is assumed that there is no mistake
  // and no internal check for the dimensionality.
  GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, index_j, given_vec);
  return given_vec.norm();
}

// inlining
// <<<<<<< HEAD
// // double GEOMETRY::get_minimum_distance(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec)
// // {
// //   // this comentation is related with the assumes that
// //   // the given dimensionality is properly working
// //   // therefore, it is assumed that there is no mistake
// //   // and no internal check for the dimensionality.
// //   GEOMETRY::get_minimum_distance_rel_vector(TRAJ, index_t, index_i, index_j, given_vec);
// //   return given_vec.norm();
// // }

// double GEOMETRY::get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check)
// {
//   GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJ, index_t, index_i, index_j, given_vec, beyond_box_check);
//   return given_vec.norm();
// }
// =======
// >>>>>>> 8692b3ce744ba0597384a10095736b1f24e127f7

// double GEOMETRY::get_minimum_distance_cell_list(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check)
// {
//   GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJ, index_t, index_i, index_j, given_vec, beyond_box_check);
//   return given_vec.norm();
// }


double GEOMETRY::return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j)
{
  MATRIX rel_vec(TRAJ.N_dimension, 1, 0.);
  GEOMETRY::get_minimum_distance(TRAJ, index_t, index_i, index_j, rel_vec);
  return rel_vec.norm();
}


// double GEOMETRY::return_minimum_distance(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j)
// {
//   MATRIX rel_vec(TRAJ.N_dimension, 1, 0.);
//   GEOMETRY::get_minimum_distance(TRAJ, index_t, index_i, index_j, rel_vec);
//   return rel_vec.norm();
// }

double GEOMETRY::minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t)
{
  double time_st = dsecnd();
  for (MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
        {
          double diff = TRAJ(target_t, i, k) - 0.5*TRAJ.box_dimension[k];
          double sign = diff/fabs(diff);
          if (fabs(diff) > 0.5*TRAJ.box_dimension[k])
            {
              TRAJ(target_t, i, k) -= sign*TRAJ.box_dimension[k];
            }
        }
    }
  return dsecnd() - time_st;
}

// MKL_LONG GEOMETRY::minimum_image_convention(TRAJECTORY_HDF5& TRAJ, MKL_LONG target_t)
// {
//   for (MKL_LONG i=0; i<TRAJ.Np; i++)
//     {
//       for (MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//         {
//           // double diff = TRAJ(target_t, i, k) - 0.5*TRAJ.box_dimension[k];
//           double diff = TRAJ.data[target_t][i][k] - 0.5*TRAJ.box_dimension[k];
//           double sign = diff/fabs(diff);
//           if (fabs(diff) > 0.5*TRAJ.box_dimension[k])
//             {
//               // TRAJ(target_t, i, k) -= sign*TRAJ.box_dimension[k];
//               TRAJ.data[target_t][i][k] -= sign*TRAJ.box_dimension[k];
//             }
//         }
//     }
//   return 0;
// }
// inlined
// MKL_LONG GEOMETRY::minimum_image_convention_particle(TRAJECTORY& TRAJ, MKL_LONG target_t, MKL_LONG index_particle)
// {
//   for(MKL_LONG k=0; k<TRAJ.dimension; k++)
//     {
//       double diff = TRAJ(target_t, index_particle, k) - 0.5*TRAJ.box_dimension[k];
//       double sign = diff/fabs(diff);
//       if (fabs(diff) > 0.5*TRAJ.box_dimension[k])
// 	{
// 	  TRAJ(target_t, index_particle, k) -= sign*TRAJ.box_dimension[k];
// 	}
//     }
//   return 0;
// }


MKL_LONG GEOMETRY::get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
{
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ(index_t, given_index, k), TRAJ(index_t, target_index, k), TRAJ.box_dimension[k]);
    }
  return 0;
}

// MKL_LONG GEOMETRY::get_minimum_distance_pos_vector(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
// {
//   for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//     {
//       given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ.data[index_t][given_index][k], TRAJ.data[index_t][target_index][k], TRAJ.box_dimension[k]);
//     }
//   return 0;
// }


// MKL_LONG GEOMETRY::get_minimum_distance_pos_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check)
// {
//   for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//     {
//       // given_vec(k) = UTIL_ARR::get_minimum_image_k_from_x(TRAJ(index_t, given_index, k), TRAJ(index_t, target_index, k), TRAJ.box_dimension[k]);
//       given_vec(k) = TRAJ(index_t, target_index, k) - TRAJ(index_t, given_index, k)
//     }
//   return 0;
// }


MKL_LONG GEOMETRY::get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
{
  GEOMETRY::get_minimum_distance_pos_vector(TRAJ, index_t, given_index, target_index, given_vec);
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      // direction convention:
      // +: direction to the given bead
      // -: direction to the target bead
      given_vec(k) -= TRAJ(index_t, given_index, k);
    }
  return 0;
}
// inlined
// double GEOMETRY::get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check)
// {
//   GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJ, index_t, index_i, index_j, given_vec, beyond_box_check);
//   return given_vec.norm();
// }


// inlined
// <<<<<<< HEAD
// // MKL_LONG GEOMETRY::get_minimum_distance_rel_vector(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec)
// // {
// //   GEOMETRY::get_minimum_distance_pos_vector(TRAJ, index_t, given_index, target_index, given_vec);
// //   for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
// //     {
// //       // direction convention:
// //       // +: direction to the given bead
// //       // -: direction to the target bead
// //       // given_vec(k) -= TRAJ(index_t, given_index, k);
// //       given_vec(k) -= TRAJ.data[index_t][given_index][k];
// //     }
// //   return 0;
// // }


// MKL_LONG GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check)
// {
//   // GEOMETRY::get_minimum_distance_pos_vector_cell_list(TRAJ, index_t, given_index, target_index, given_vec, beyond_box_check);
//   // for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//   //   {
//   //     // direction convention:
//   //     // +: direction to the given bead
//   //     // -: direction to the target bead
//   //     given_vec(k) -= TRAJ(index_t, given_index, k);
//   //   }
//   for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//     {
//       double PBC_coord_target = TRAJ(index_t, target_index, k) + (double)beyond_box_check[k]*TRAJ.box_dimension[k];
//       given_vec(k) = PBC_coord_target - TRAJ(index_t, given_index, k);
//     }
//   return 0;
// }
// =======
// // inlined
// // MKL_LONG GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check)
// // {
// //   // GEOMETRY::get_minimum_distance_pos_vector_cell_list(TRAJ, index_t, given_index, target_index, given_vec, beyond_box_check);
// //   // for(MKL_LONG k=0; k<TRAJ.dimension; k++)
// //   //   {
// //   //     // direction convention:
// //   //     // +: direction to the given bead
// //   //     // -: direction to the target bead
// //   //     given_vec(k) -= TRAJ(index_t, given_index, k);
// //   //   }
// //   for(MKL_LONG k=0; k<TRAJ.dimension; k++)
// //     {
// //       double PBC_coord_target = TRAJ(index_t, target_index, k) + (double)beyond_box_check[k]*TRAJ.box_dimension[k];
// //       given_vec(k) = PBC_coord_target - TRAJ(index_t, given_index, k);
// //     }
// //   return 0;
// // }
// >>>>>>> 8692b3ce744ba0597384a10095736b1f24e127f7


// MKL_LONG GEOMETRY::get_minimum_distance_rel_vector_cell_list(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check)
// {
//   // GEOMETRY::get_minimum_distance_pos_vector_cell_list(TRAJ, index_t, given_index, target_index, given_vec, beyond_box_check);
//   // for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//   //   {
//   //     // direction convention:
//   //     // +: direction to the given bead
//   //     // -: direction to the target bead
//   //     given_vec(k) -= TRAJ(index_t, given_index, k);
//   //   }
//   for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//     {
//       // double PBC_coord_target = TRAJ(index_t, target_index, k) + (double)beyond_box_check[k]*TRAJ.box_dimension[k];
//       double PBC_coord_target = TRAJ.data[index_t][target_index][k] + (double)beyond_box_check[k]*TRAJ.box_dimension[k];
//       given_vec(k) = PBC_coord_target - TRAJ.data[index_t][given_index][k];
//     }
//   return 0;
// }


double UTIL_ARR::get_minimum_image_k_from_x(double x, double k, double dimension)
 {
   double kd[3] = {k-dimension - x, k - x, k+dimension - x};
   double re= kd[get_index_minimum_abs(kd, 3)] + x;
   return re;
 }

MKL_LONG UTIL_ARR::get_index_minimum_abs(double *k, MKL_LONG N)
 {
   MKL_LONG re = 0;
   for(MKL_LONG i=1; i<N; i++)
     {
       if (fabs(k[i]) < fabs(k[re]))
         re = i;
     }
   return re;
 }

double GEOMETRY::get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i)
{
  double distance = 0.;
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      distance += pow(TRAJ(index_t, index_i, k) - TRAJ(index_t, index_i, k), 2.0);
    }
  return sqrt(distance);

}

// double GEOMETRY::get_simple_distance(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i)
// {
//   double distance = 0.;
//   for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
//     {
//       // distance += pow(TRAJ(index_t, index_i, k) - TRAJ(index_t, index_i, k), 2.0);
//       distance += pow(TRAJ.data[index_t][index_i][k] - TRAJ.data[index_t][index_i][k], 2.0);
//     }
//   return sqrt(distance);

// }
