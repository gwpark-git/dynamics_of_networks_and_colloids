
#ifndef LIB_POTENTIAL_H
#define LIB_POTENTIAL_H
#include <iostream>
#include <math.h>
#include <mkl.h>
#include "read_file_condition.h"
#include "matrix_ed.h"
#include <string>

class POTENTIAL_SET
{
 public:
  double *force_variables;
  double (*f_connector)(double distance, double* given_variables);
  double (*e_connector)(double distance, double* given_variables);
  double (*f_repulsion)(double distance, double* given_variables);
  double (*e_repulsion)(double distance, double* given_variables);
  double (*scale_random)(MATRIX& basic_random_var_unity, double* given_variables);
  double (*w_function)(double distance, double tension, double* given_variables);

  POTENTIAL_SET()
    {
    }
  /* POTENTIAL_SET(COND& given_cond); */
  virtual ~POTENTIAL_SET()
    {
      if(force_variables)
        delete[] force_variables;
    }
};

namespace FORCE
{
  namespace DEFAULT
  {
    double EMPTY_force_contribution(double distance, double *given_variables);
    double basic_random(MATRIX& given_basic_random, double* given_variables);
    double time_scaling_random(MATRIX& given_basic_random, double scale_factor);
    MKL_LONG EMPTY_force_set(POTENTIAL_SET& given_POT, COND& given_condition);

    double Gaussian_spring_force(double distance, double N_dimension);
    double Gaussian_spring_potential(double distance, double N_dimension);
  }
  namespace NAPLE
  {
    double excluded_volume_force(double distance, double effective_distance);
    double excluded_volume_potential(double distance, double effective_distance);
    
    namespace SIMPLE_REPULSION
    {
      double MAP_excluded_volume_force(double distance, double* given_variables);
      double MAP_excluded_volume_potential(double distance, double* given_variables);
      double MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables);
      MKL_LONG MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond);
    }
    namespace MC_ASSOCIATION
    {
      double MAP_Gaussian_spring_force(double distance, double* given_variables);
      double MAP_Gaussian_spring_potential(double distance, double* given_variables);
      double MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables);
      MKL_LONG MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond);
      double Detachment_weight(double distance, double tension, double* force_variables);
    }
    
  }
}  


#endif
