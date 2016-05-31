
#ifndef RANDOM_H
#define RANDOM_H

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_mksa.h>
}

#include <math.h>
#include <mkl.h>
#include <iostream>
#include "matrix.h"
#include "file_IO.h"

// note that the RANDOM package should not has the dependency with lib_traj.h
class RNG_BOOST
{
 public:
  MKL_LONG N_THREADS_BD;
  MKL_LONG N_THREADS_SS;
  gsl_rng **BOOST_BD;
  gsl_rng **BOOST_SS;
  bool INITIALIZATION_BD;
  bool INITIALIZATION_SS;
  RNG_BOOST()
    {
      std::cout<< "ERR: basic constructor for RNG_BOOST class is not allowed\n";
    }
  RNG_BOOST(COND& given_condition)
    {
      INITIALIZATION_BD = TRUE;
      const gsl_rng_type *T_boost;
      N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
      /* BOOST_BD = (gsl_rng**)mkl_malloc(N_THREADS_BD*sizeof(gsl_rng*), BIT); */
      BOOST_BD = new gsl_rng* [N_THREADS_BD];
      gsl_rng_env_setup();
      T_boost = gsl_rng_default;
      for(MKL_LONG i=0; i<N_THREADS_BD; i++)
        {
          BOOST_BD[i] = gsl_rng_alloc(T_boost);
          MKL_LONG seed_BD = atol(given_condition("basic_random_seed").c_str());
          gsl_rng_set(BOOST_BD[i], seed_BD + i); // it set the seed with index i
        }
      if(given_condition("Step")!="EQUILIBRATION")
        {
          INITIALIZATION_SS = TRUE;
          N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());
          /* BOOST_SS = (gsl_rng**)mkl_malloc(N_THREADS_SS*sizeof(gsl_rng*), BIT); */
          BOOST_SS = new gsl_rng* [N_THREADS_SS];
          for(MKL_LONG i=0; i<N_THREADS_SS; i++)
            {
              BOOST_SS[i] = gsl_rng_alloc(T_boost);
              MKL_LONG seed_SS = atol(given_condition("basic_random_seed_SS").c_str());
              gsl_rng_set(BOOST_SS[i], seed_SS + i); 
            }
        }
    }

  ~RNG_BOOST()
    {
      if(INITIALIZATION_BD)
        {
          for(MKL_LONG i=0; i<N_THREADS_BD; i++)
            {
              gsl_rng_free(BOOST_BD[i]);
            }
          /* mkl_free(BOOST_BD); */
          delete[] BOOST_BD;
          if(INITIALIZATION_SS)
            {
              for(MKL_LONG i=0; i<N_THREADS_SS; i++)
                {
                  gsl_rng_free(BOOST_SS[i]);
                }
              /* mkl_free(BOOST_SS); */
              delete[] BOOST_SS;
            }
        }
    }
};


namespace RANDOM
{
  MKL_LONG random_vector_generator(MATRIX& R_VEC_TRANS);
  MKL_LONG single_random_vector_generator(MATRIX& given_vec);
  MKL_LONG single_random_vector_generator_variance(MATRIX& given_vec, double s_2);
  MKL_LONG single_random_vector_generator_boost(MATRIX& given_vec, gsl_rng* r_boost);
  MKL_LONG single_random_vector_generator_variance_boost(MATRIX& given_vec, double s_2, gsl_rng* r_boost);
  MKL_LONG single_unit_random_vector_generator(MATRIX& given_vec);
  MKL_LONG unit_random_vector_generator(MATRIX& R_VEC_TRANS);
  MKL_LONG unit_random_vector_generator_2D(MATRIX& R_VEC_TRANS);

  MKL_LONG return_LONG_INT_rand(MKL_LONG SUP);
  MKL_LONG return_LONG_INT_rand_boost(gsl_rng* r, MKL_LONG SUP);
  double return_double_rand_SUP1();
  double return_double_rand_SUP1_boost(gsl_rng* r);
  MKL_LONG get_LONG_ARR_rand_boost(gsl_rng* r, MKL_LONG SUP, MKL_LONG* given_long_arr, MKL_LONG N_arr);
  MKL_LONG get_DOUBLE_ARR_rand_boost(gsl_rng* r, double* given_double_arr, MKL_LONG N_arr);
  
}


#endif 
