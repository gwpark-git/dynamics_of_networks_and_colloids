
#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <iostream>
#include <math.h>
#include <mkl.h>
/* #include "read_file_condition.h" */
#include "file_IO.h"
#include "matrix.h"
#include <string>
#define BIT 64

struct POTENTIAL_VARIABLE
{
  /*
    This struct replace the 'double* force_variable' where it contains various variable including repulsion coefficient, effective distance for repulsion, spatial dimensionality, l_cap, and so on. 
    For instance, the current FORCE::NAPLE::MC_ASSOCIATION::MAP_potential set uses the following list (part of it):
  given_POT.force_variables[0] = atof(given_cond("repulsion_coefficient").c_str());
  given_POT.force_variables[1] = atof(given_cond("effective_distance").c_str());
  given_POT.force_variables[2] = 1./sqrt(given_POT.force_variables[0]);
  given_POT.force_variables[3] = atol(given_cond("N_dimension").c_str());
  given_POT.force_variables[4] = atof(given_cond("l_cap").c_str());
   */
  double repulsion_coefficient;  // C_rep
  double effective_distance;
  double inv_sqrt_repulsion_coefficient; // 1/C_rep
  MKL_LONG N_dimension;
  double l_cap; // l_cap
  double scale_factor_chain; // alpha
  double ratio_RM_R0; // R_max
  double cutoff_connection;
  double delta_t0;
  double modified_friction_tau0; // this is temporal method
  double energy_barrier; // this is temporal method

MKL_LONG constructor_POTENTIAL_VARIABLE()
{
  // empty constructor
  repulsion_coefficient = 0;
  effective_distance = 0;
  inv_sqrt_repulsion_coefficient=0;
  N_dimension = 0;
  l_cap = 0;
  scale_factor_chain = 0;
  ratio_RM_R0 = 0;
  cutoff_connection = 0;
  delta_t0 = 0;
  modified_friction_tau0 = 0;
  energy_barrier = 0;
    
  return 0;
}

MKL_LONG destructor_POTENTIAL_VARIABLE()
{
  // empty destructor
  return 0;
}

  
};


class POTENTIAL_SET
{
 public:
  
  // member variables
  // double *force_variables;
  POTENTIAL_VARIABLE force_variables;

  double
  (*zeta_particle)
  (MKL_LONG N_connectors, POTENTIAL_VARIABLE& given_variables);
  
  double
    (*f_connector)
    (double distance, POTENTIAL_VARIABLE& given_variables);
  double
    (*e_connector)
    (double distance, POTENTIAL_VARIABLE& given_variables);
  double
    (*PDF_connector)
    (double distance, POTENTIAL_VARIABLE& given_variables);
  double
    (*f_repulsion)
    (double distance, POTENTIAL_VARIABLE& given_variables);
  double
    (*e_repulsion)
    (double distance, POTENTIAL_VARIABLE& given_variables);
  double
    (*scale_random)
    (MATRIX& basic_random_var_unity, POTENTIAL_VARIABLE& given_variables);
  double
    (*w_function)
  (double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
  
    // (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
  double
    (*transition)
  (double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
    // (double distance, double tension, double* given_varialbes);


  POTENTIAL_SET()
    {
      force_variables.constructor_POTENTIAL_VARIABLE();
    }
  virtual
    ~POTENTIAL_SET()
    {
      // if(force_variables) // note that without initialization, it always have NULL pointer
      //   delete[] force_variables;
      // // note that the force_variables is replaced by structure variables.
      force_variables.destructor_POTENTIAL_VARIABLE();
    }
};


namespace FORCE
{
  namespace CUTOFF_ASSOCIATION
  {
    double
    cutoff_equal_probability(double distance, double effective_distance);

  }
  namespace BROWNIAN
  {

    MKL_LONG MAP_potential_variable(POTENTIAL_SET& given_POT, COND& given_cond);
    MKL_LONG
      MAP_potential_set
      (POTENTIAL_SET& given_POT, COND& given_cond);
    double
      no_time_scaling_random
      (MATRIX& given_basic_random, double scale_factor);
    double
      MAP_no_time_scaling_random
      (MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables);
  }
  namespace DUMBBELL
  {
    MKL_LONG MAP_potential_variable(POTENTIAL_SET& given_POT, COND& given_cond);
    
    MKL_LONG
      MAP_potential_set
      (POTENTIAL_SET& given_POT, COND& given_cond);
    double
      MAP_modified_Gaussian_spring_force
      (double distance, POTENTIAL_VARIABLE& given_variables);
    double
      MAP_modified_Gaussian_spring_potential
      (double distance, POTENTIAL_VARIABLE& given_variables);

  }
  namespace GAUSSIAN
  {
    double
      spring_force
      (double distance, double N_dimension);
    double
      spring_potential
      (double distance, double N_dimension);
    double
      Boltzmann_distribution
      (double distance, double N_dimension);
  }
  namespace MODIFIED_GAUSSIAN
  {
    double
      spring_force
      (double distance, double N_dimension, double scale_factor);
    double
      spring_potential
      (double distance, double N_dimension, double scale_factor);
    double
      Boltzmann_distribution
      (double distance, double N_dimension, double scale_factor);
    double
      cutoff_Boltzmann_distribution
      (double distance, double N_dimension, double scale_factor, double cutoff_radius); 
  }
  namespace FENE
  {
    double
      non_Gaussian_factor
      (double distance, double N_dimension, double ratio_RM_R0);
    double
      spring_force
    (double distance, double N_dimension, double scale_factor, double ratio_RM_R0);
// =======
//     spring_force
//       (double distance, double N_dimension, double ratio_RM_R0);
// >>>>>>> step_shear
    double
      spring_potential
    (double distance, double N_dimension, double scale_factor, double ratio_RM_R0);
    double
      Boltzmann_distribution
    (double distance, double N_dimension, double scale_factor, double ratio_RM_R0);
  }
  namespace DEFAULT
  {
    double
      EMPTY_force_contribution
      (double distance, POTENTIAL_VARIABLE& given_variables);
    double
      basic_random
      (MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables);
    double
      time_scaling_random
      (MATRIX& given_basic_random, double scale_factor);
    MKL_LONG
      EMPTY_force_set
      (POTENTIAL_SET& given_POT, COND& given_condition);

    double
    dimensionless_friction
    (MKL_LONG N_connectors, POTENTIAL_VARIABLE& given_variables);

  }
  namespace NAPLE
  {
    double
    excluded_volume_force
      (double distance, double effective_distance);
    double
      excluded_volume_potential
      (double distance, double effective_distance);
    
    namespace SIMPLE_REPULSION
    {

      MKL_LONG MAP_potential_variable(POTENTIAL_SET& given_POT, COND& given_cond);
      double
	MAP_excluded_volume_force
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
	MAP_excluded_volume_potential
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
	MAP_time_scaling_random
	(MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables);
      MKL_LONG
	MAP_potential_set
	(POTENTIAL_SET& given_POT, COND& given_cond);
    }
    namespace MC_ASSOCIATION
    {

      MKL_LONG MAP_potential_variable(POTENTIAL_SET& given_POT, COND& given_cond);
      double
      friction_LOOP_DISSOCIATION_TIME
      (MKL_LONG N_connectors, POTENTIAL_VARIABLE& given_variables);

      
      double
	MAP_Gaussian_spring_force
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
	MAP_Gaussian_spring_potential
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
      MAP_Gaussian_Boltzmann
	(double distance, POTENTIAL_VARIABLE& given_variables);

      double
      MAP_minimum_R0_Gaussian_Boltzmann
      (double distance, POTENTIAL_VARIABLE& given_variables);

      double
	MAP_modified_Gaussian_spring_force
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
	MAP_modified_Gaussian_spring_potential
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
      MAP_modified_Gaussian_Boltzmann
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
      MAP_minimum_R0_modified_Gaussian_Boltzmann
	(double distance, POTENTIAL_VARIABLE& given_variables);
      
      double
	MAP_cutoff_modified_Gaussian_Boltzmann
	(double distance, POTENTIAL_VARIABLE& given_variables);

      double
      MAP_minimum_R0_cutoff_modified_Gaussian_Boltzmann
	(double distance, POTENTIAL_VARIABLE& given_variables);
      
      double
      MAP_cutoff_equal_probability
      (double distance, POTENTIAL_VARIABLE& given_variables);

      
      double
	MAP_FENE_spring_force
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
	MAP_FENE_spring_potential
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
      MAP_FENE_Boltzmann
	(double distance, POTENTIAL_VARIABLE& given_variables);
      double
	MAP_minimum_R0_FENE_Boltzmann
	(double distance, POTENTIAL_VARIABLE& given_variables);      
      double
	MAP_time_scaling_random
	(MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables);
      MKL_LONG
	MAP_potential_set
	(POTENTIAL_SET& given_POT, COND& given_cond);
    }

  }
}  

namespace KINETICS
{
  namespace UNIFORM
  {
    double detachment_weight(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
    double transition_probability(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
  }
  namespace WEIGHTED
  {
    double detachment_weight(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
    double transition_probability(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
  }

  namespace METROPOLIS
  {
    double detachment_weight(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
    double transition_probability(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, double* given_varialbes);
  }

  double dissociation_probability(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
  double cutoff_dissociation_probability(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
  double minimum_R0_dissociation_probability(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);

  // (double distance, double tension, doouble* given_variables);
  double dissociation_probability_equal_modified_gaussian(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
  
  namespace FIRST_ORDER
  {
    double dissociation_probability(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables);

    // (double distance, double tension, POTENTIAL_VARIABLE& given_variables);
  }

}


// inline functions
// note that all the inline function 'must be' defined on above for readability.

inline double
KINETICS::WEIGHTED::
detachment_weight
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  // given_variables[3] == Nd
  // given_variables[4] == l_cap
  // in normalized scheme, the energy barrier is canceled out with normalization factor
  // therefore, it is only affected by the tension exerted on the chain
  double tension = force(distance, given_variables);
  // return exp(tension*given_variables[4]);
  return exp(tension*given_variables.l_cap);
  
}


inline double
KINETICS::METROPOLIS::
detachment_weight
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  return 1.0;
}

inline double
KINETICS::METROPOLIS::
transition_probability
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  double tension = force(distance, given_variables);
  // double tpa = exp(tension*given_variables[4] - given_variables[6]);
  double tpa = exp(tension*given_variables.l_cap - given_variables.energy_barrier);
  if (tpa > 1.0)
    return 1.0;
  return tpa;
}

inline double
KINETICS::WEIGHTED::
transition_probability
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  double tension = force(distance, given_variables);

  double tpa = exp(tension*given_variables.l_cap);
  if (tpa > 1.0)
    return 1.0;
  return tpa;
}

inline double
KINETICS::UNIFORM::
detachment_weight
// (double distance, double tension, double *given_variables)
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
{
  return 1.0;
}

inline double
KINETICS::UNIFORM::
transition_probability
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, double* given_varialbes)
{
  return 1.0;
}

inline double
KINETICS::
dissociation_probability
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  double tension = force(distance, given_variables);

  // double tpa = given_variables[6]*exp(tension*given_variables[4]);
  double tpa = given_variables.delta_t0*exp(tension*given_variables.l_cap);
  
  if (tpa > 1.0)
    return 1.0;
  return tpa;
}

inline double
KINETICS::
minimum_R0_dissociation_probability
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  double tension = 0.;
  if (distance < 1.0)
    tension = force(1.0, given_variables);
  else
    tension = force(distance, given_variables);

  // double tpa = given_variables[6]*exp(tension*given_variables[4]);
  double tpa = given_variables.delta_t0*exp(tension*given_variables.l_cap);

  if (tpa > 1.0)
    return 1.0;
  return tpa;
}


inline double
KINETICS::
cutoff_dissociation_probability
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  // double tpa = given_variables[6]*exp(tension*given_variables[4]);
  // double tpa = given_variables[6];
  double tpa = given_variables.delta_t0;
  if (tpa > 1.0 || distance >= 1.0)
    return 1.0;
  return tpa;
}

inline double
KINETICS::
dissociation_probability_equal_modified_gaussian
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  /*
    Unlike the previous approach, it set the dissociation probability as the same with Boltzman distribution without normalization.
    In consequence, the detachment frequency have the same form with the association map.
   */
  // double tpa = given_variables[6]*exp(FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_potential(distance, given_variables));
  double tpa = given_variables.delta_t0*exp(FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_potential(distance, given_variables));
  
  if (tpa > 1.0)
    return 1.0;
  return tpa;
}

inline double
KINETICS::FIRST_ORDER::
dissociation_probability
(double distance, double (*force)(double, POTENTIAL_VARIABLE&), POTENTIAL_VARIABLE& given_variables)
// (double distance, double tension, POTENTIAL_VARIABLE& given_variables)
{
  double tension = force(distance, given_variables);
  // double Dt = given_variables[6];
  // double lcap = given_variables[4];
  // double beta = exp(tension*lcap);
  // double tpa = 1.0 - exp(-beta*Dt);
  double beta = exp(tension*given_variables.l_cap);
  double tpa = 1.0 - exp(-beta*given_variables.delta_t0);

  
  // double tpa = given_variables[6]*exp(tension*given_variables[4]);
  // if (tpa > 1.0)
  //   return 1.0;
  return tpa;
}


inline double
FORCE::DEFAULT::
EMPTY_force_contribution
// (double distance, double* given_varialbes)
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return 0;
}

inline double
FORCE::DEFAULT::
basic_random
(MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables)
{
  return 0;
}


inline double
FORCE::DEFAULT::
time_scaling_random
(MATRIX& given_basic_random, double scale_factor)
{
  matrix_mul(given_basic_random, scale_factor); // this is varied due to the given time scale
  return scale_factor;
}

inline double
FORCE::DEFAULT::dimensionless_friction
(MKL_LONG N_connectors, POTENTIAL_VARIABLE& given_variables)
{
  return 1.; // dimensionless contribution in basic part will be unity.
}


inline double
FORCE::NAPLE::SIMPLE_REPULSION::
MAP_time_scaling_random
(MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables)
{
  // return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
  return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables.inv_sqrt_repulsion_coefficient);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_time_scaling_random
(MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables)
{
  // return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
  return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables.inv_sqrt_repulsion_coefficient);
  
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
friction_LOOP_DISSOCIATION_TIME
(MKL_LONG N_connectors, POTENTIAL_VARIABLE& given_variables)
{
  if (N_connectors == 0)
    return 1.;
  // return N_connectors*given_variables[9]/2.; // the denominator 2 related with correction for multiple counter
  return N_connectors*given_variables.modified_friction_tau0/2.; // the denominator 2 related with correction for multiple counter
  
}



inline double
FORCE::NAPLE::
excluded_volume_force
(double distance, double effective_distance)
{
  if (distance < effective_distance)
    {
      return -(1. - pow(distance, 2.0));
    }
  return 0.;
}

inline double
FORCE::NAPLE::
excluded_volume_potential
(double distance, double effective_distance)
{
  if (distance < effective_distance)
      return (1./3.)*pow(1. - distance, 2.)*(2. + distance);
  return 0.;
}

inline double
FORCE::NAPLE::SIMPLE_REPULSION::
MAP_excluded_volume_force
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  // return FORCE::NAPLE::excluded_volume_force(distance, given_variables[1]);
  return FORCE::NAPLE::excluded_volume_force(distance, given_variables.effective_distance);
}

inline double
FORCE::NAPLE::SIMPLE_REPULSION::
MAP_excluded_volume_potential
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::NAPLE::excluded_volume_potential(distance, given_variables.effective_distance);
}


inline double
FORCE::GAUSSIAN::
spring_force
(double distance, double N_dimension)
{
  return N_dimension*distance;
}

inline double
FORCE::GAUSSIAN::
spring_potential
(double distance, double N_dimension)
{
  return (0.5)*N_dimension*distance*distance;
}

inline double
FORCE::GAUSSIAN::
Boltzmann_distribution
(double distance, double N_dimension)
{
  return exp(-spring_potential(distance, N_dimension));
}

inline double
FORCE::MODIFIED_GAUSSIAN::
spring_force
(double distance, double N_dimension, double scale_factor)
{
  return FORCE::GAUSSIAN::spring_force(scale_factor*scale_factor*distance, N_dimension);
}

inline double
FORCE::MODIFIED_GAUSSIAN::
spring_potential
(double distance, double N_dimension, double scale_factor)
{
  return FORCE::GAUSSIAN::spring_potential(scale_factor*distance, N_dimension);
}

inline double
FORCE::MODIFIED_GAUSSIAN::
Boltzmann_distribution
(double distance, double N_dimension, double scale_factor)
{
  return FORCE::GAUSSIAN::Boltzmann_distribution(scale_factor*distance, N_dimension);
}

inline double
FORCE::MODIFIED_GAUSSIAN::
cutoff_Boltzmann_distribution
(double distance, double N_dimension, double scale_factor, double cutoff_radius)
{
  // note that this cut-off do not have any benefit at this moment
  // after implementaion for cell-list, the cut-off scheme becomes efficience
  if (distance < cutoff_radius)
    return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(distance, N_dimension, scale_factor);
  return 0.;
}

inline double
FORCE::FENE::
non_Gaussian_factor
(double distance, double N_dimension, double ratio_RM_R0)
{
  // return 1.0/(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
  if (distance < ratio_RM_R0)
    return 1.0/(1.0 - pow(distance/ratio_RM_R0, 2.0));
  else
    return 0.;
  // return ratio_RM_R0/(pow(ratio_RM_R0, 2.0) - pow(distance, 2.0));
}

inline double
FORCE::FENE::
spring_force
(double distance, double N_dimension, double scale_factor, double ratio_RM_R0)
{
  return non_Gaussian_factor(distance, N_dimension, ratio_RM_R0)*FORCE::GAUSSIAN::spring_force(scale_factor*scale_factor*distance, N_dimension);
}

inline double
FORCE::FENE::
spring_potential
(double distance, double N_dimension, double scale_factor, double ratio_RM_R0)
{
  if (distance < ratio_RM_R0)
    return (-(double)N_dimension/2.0)*pow(scale_factor*ratio_RM_R0, 2.0)*log(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
  else
    {
      printf("ERR: FENE::the length of spring (%lf) cannot be over than finite extensibility (%lf)\n", distance, ratio_RM_R0);
      return 0.;
    }
}

inline double
FORCE::FENE::
Boltzmann_distribution
(double distance, double N_dimension, double scale_factor, double ratio_RM_R0)
{
  // Note that even without cut-off range based on the ratio_RM_R0, the code is properly working.
  // However, when we applied this scheme, the sorting procedure is more stable and normalized in naturally.
  if (distance < ratio_RM_R0)
    return exp(-FORCE::FENE::spring_potential(distance, N_dimension, scale_factor, ratio_RM_R0));
  return 0.;
}

inline double
FORCE::CUTOFF_ASSOCIATION::
cutoff_equal_probability
(double distance, double effective_distance)
{
  if (distance < effective_distance)
    return 1.;
  return 0.;
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_Gaussian_spring_force
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::GAUSSIAN::spring_force(distance, given_variables.N_dimension);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_Gaussian_spring_potential
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::GAUSSIAN::spring_potential(distance, given_variables.N_dimension);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_Gaussian_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::GAUSSIAN::Boltzmann_distribution(distance, given_variables.N_dimension);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_minimum_R0_Gaussian_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  if (distance < 1.0)
    FORCE::GAUSSIAN::Boltzmann_distribution(1.0, given_variables.N_dimension);
  return FORCE::GAUSSIAN::Boltzmann_distribution(distance, given_variables.N_dimension);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_modified_Gaussian_spring_force
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::spring_force(distance, given_variables.N_dimension, given_variables.scale_factor_chain);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_modified_Gaussian_spring_potential
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::spring_potential(distance, given_variables.N_dimension, given_variables.scale_factor_chain);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_modified_Gaussian_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(distance, given_variables.N_dimension, given_variables.scale_factor_chain);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_minimum_R0_modified_Gaussian_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  if (distance < 1.0)
    return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(1.0, given_variables.N_dimension, given_variables.scale_factor_chain);
  return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(distance, given_variables.N_dimension, given_variables.scale_factor_chain);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_cutoff_equal_probability
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::CUTOFF_ASSOCIATION::cutoff_equal_probability(distance, given_variables.cutoff_connection);
}



inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_cutoff_modified_Gaussian_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::cutoff_Boltzmann_distribution(distance, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.cutoff_connection);
}


inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_minimum_R0_cutoff_modified_Gaussian_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  if (distance < 1.0)
    return FORCE::MODIFIED_GAUSSIAN::cutoff_Boltzmann_distribution(1.0, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.cutoff_connection);
  return FORCE::MODIFIED_GAUSSIAN::cutoff_Boltzmann_distribution(distance, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.cutoff_connection);
}


inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_FENE_spring_force
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::FENE::spring_force(distance, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.ratio_RM_R0);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_FENE_spring_potential
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::FENE::spring_potential(distance, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.ratio_RM_R0);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_FENE_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::FENE::Boltzmann_distribution(distance, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.ratio_RM_R0);
}

inline double
FORCE::NAPLE::MC_ASSOCIATION::
MAP_minimum_R0_FENE_Boltzmann
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  if (distance < 1.0)
    FORCE::FENE::Boltzmann_distribution(1.0, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.ratio_RM_R0);
  return FORCE::FENE::Boltzmann_distribution(distance, given_variables.N_dimension, given_variables.scale_factor_chain, given_variables.ratio_RM_R0);
}

inline double
FORCE::DUMBBELL::
MAP_modified_Gaussian_spring_force
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::spring_force(distance, given_variables.N_dimension, given_variables.scale_factor_chain);
}

inline double
FORCE::DUMBBELL::
MAP_modified_Gaussian_spring_potential
(double distance, POTENTIAL_VARIABLE& given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::spring_potential(distance, given_variables.N_dimension, given_variables.scale_factor_chain);
}


#endif
