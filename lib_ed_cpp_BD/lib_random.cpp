
#include "lib_random.h"


long RANDOM::return_LONG_INT_rand(long SUP)
{
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  // srandom(0); // for seeding of randome generation for the seed
  gsl_rng_set(r, random());
  long re = gsl_rng_uniform_int(r, SUP);
  gsl_rng_free(r);
  return re;  
}

long RANDOM::return_LONG_INT_rand_boost(gsl_rng* r, long SUP)
{
  return gsl_rng_uniform_int(r, SUP);
}

long RANDOM::get_LONG_ARR_rand_boost(gsl_rng* r, long SUP, long* given_long_arr, long N_arr);
{
  for(long i=0; i<N_arr; i++)
    {
      given_long_arr[i] = gsl_rng_uniform_int(r, SUP);
    }
  return 0;  
}


double RANDOM::return_double_rand_SUP1_boost(gsl_rng* r)
{
  return gsl_rng_uniform(r);
}

double RANDOM::get_DOUBLE_ARR_rand_boost(gsl_rng* r, double* given_double_arr, long N_arr);
{
  for(long i=0; i<N_arr; i++)
    {
      given_double_arr[i] = gsl_rng_uniform(r);
    }
  return 0;
}


double RANDOM::return_double_rand_SUP1()
{
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  // srandom(0); // for seeding of randome generation for the seed
  gsl_rng_set(r, random());
  double re = gsl_rng_uniform(r);
  gsl_rng_free(r);
  return re;  
}

// long RANDOM::dumbbell_connectivity(TRAJECTORY& TRAJ)
// {
//   TRAJ.CONNECTIVITY.initial(TRAJ.Np, TRAJ.Np, 0.);
//   for(long i=0; i<TRAJ.Np -1; i+= 2)
//     {
//       TRAJ.CONNECTIVITY(i, i+1) = TRUE; // forward connection
//       TRAJ.CONNECTIVITY(i+1, i) = TRUE; // backward connection
//     }
//   return 0;
// }

// long RANDOM::random_position_dumbbell_generator(TRAJECTORY& TRAJ)
// {
//   RANDOM::random_position_generator(TRAJ);
//   RANDOM::dumbbell_connectivity(TRAJ);
//   MATRIX tmp(3,1, 0.);
//   for(long i=1; i<TRAJ.Np; i+= 2)
//     {
//       RANDOM::single_unit_random_vector_generator(tmp);
//       for(long k=0; k<TRAJ.dimension; k++)
//         {
//           TRAJ(0, i, k) = TRAJ(0, i-1, k) + tmp(k);
//           // *(TRAJ.R_ref[0](i,k)) = *(TRAJ.R_ref[0](i-1, k)) + tmp(k);
//         }
//     }
//   return 0;
// }


long RANDOM::single_random_vector_generator(MATRIX& given_vec)
{
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  // srandom(0); // for seeding of randome generation for the seed
  gsl_rng_set(r, random());
  for (long k=0; k<given_vec.size; k++)
    {
      given_vec(k) = 2.0*(gsl_rng_uniform(r) - 0.5); 
      // -0.5 is used for vector components from -0.5 to 0.5, then multiplied by 2 gave us -1 to 1. 
    }
  gsl_rng_free(r);
  return 0;
}

long RANDOM::single_random_vector_generator_variance(MATRIX& given_vec, double s_2)
{
  single_random_vector_generator(given_vec);
  for(long k=0; k<given_vec.size; k++)
    {
      given_vec(k) *= sqrt(3.0)*s_2; // to get variance 1
    }
  return 0;
}

long RANDOM::random_vector_generator(MATRIX& R_VEC_TRANS)
{
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  // srandom(1); // for seeding of randome generation for the seed

  for (long k=0; k<R_VEC_TRANS.cols; k++)
    {
      gsl_rng_set(r, random());
      for(long i=0; i<R_VEC_TRANS.rows; i++)
        {
          R_VEC_TRANS(i, k) = gsl_rng_uniform(r) - 0.5; // -0.5 is used for vector components from -0.5 to 0.5
        }
    }
  gsl_rng_free(r);
  return 0;
}



long RANDOM::single_unit_random_vector_generator(MATRIX& given_vec)
{
  double norm = 0.;
  do
    {
      single_random_vector_generator(given_vec);
      norm = given_vec.norm();
    } while (norm > 1.);

  matrix_mul(given_vec, 1./norm);
  // for(long k=0; k<given_vec.size; k++)
  //   {
  //     given_vec(k) /= norm;
  //   }

  return 0;
}

long RANDOM::unit_random_vector_generator(MATRIX& R_VEC_TRANS)
{
  // double norm = 2.;
  // this selection is of importance to avoid over-generating for diagonal components.
  // however, this box generation scheme is very efficiency compared with other well-known algorithm; therefore, it would be better to 
  // while(norm <= 1.0) 
  //   {
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  // srandom(0); // for seeding of randome generation for the seed
  gsl_rng_set(r, random());
  double norm = 0.;
  for(long i=0; i<R_VEC_TRANS.rows; i++)
    {
      do
        {
          norm = 0.;
          for (long k=0; k<R_VEC_TRANS.cols; k++)
            {
              R_VEC_TRANS(i, k) = gsl_rng_uniform(r) - 0.5; // -0.5 is used for vector components from -0.5 to 0.5
              norm += pow(R_VEC_TRANS(i,k),2.);
            }
          norm = sqrt(norm);
        } while (norm > 1.);
      for(long k=0; k<R_VEC_TRANS.cols; k++)
        {
          R_VEC_TRANS(i,k) /= norm;
        }

    }
  gsl_rng_free(r);

  return 0;
}

long RANDOM::unit_random_vector_generator_2D(MATRIX& R_VEC_TRANS)
{
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, random());
  double pi = M_PI;
  long i=0;
  for(i=0; i<R_VEC_TRANS.rows; i++)
    {
      double theta = gsl_rng_uniform(r)*2.0*pi;
      R_VEC_TRANS(i, 0) = cos(theta);
      R_VEC_TRANS(i, 1) = sin(theta);
    }
  gsl_rng_free(r);
  return 0;
}

