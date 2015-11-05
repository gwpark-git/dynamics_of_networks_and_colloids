
#include "lib_potential.h"


double FORCE::DEFAULT::EMPTY_force_contribution(double distance, double* given_varialbes)
{
  return 0;
}

double FORCE::DEFAULT::basic_random(MATRIX& given_basic_random, double* given_variables)
{
  return 0;
}

long FORCE::DEFAULT::EMPTY_force_set(POTENTIAL_SET& given_POT, COND& given_condition)
{
  given_POT.f_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.f_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.scale_random = FORCE::DEFAULT::basic_random;
  return 0;
}

double FORCE::DEFAULT::time_scaling_random(MATRIX& given_basic_random, double scale_factor)
{
  matrix_mul(given_basic_random, scale_factor); // this is varied
  return scale_factor;
}

double FORCE::NAPLE::SIMPLE_REPULSION::MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables)
{
  return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
}

double FORCE::NAPLE::MC_ASSOCIATION::MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables)
{
  return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
}


double FORCE::NAPLE::excluded_volume_force(double distance, double effective_distance)
{
  if (distance < effective_distance)
    {
      return -(1. - pow(distance, 2.0));
    }
  return 0.;
}

double FORCE::NAPLE::excluded_volume_potential(double distance, double effective_distance)
{
  if (distance < effective_distance)
      return (1./3.)*pow(1. - distance, 2.)*(2. + distance);
  return 0.;
}

double FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_force(double distance, double* given_variables)
{
  return FORCE::NAPLE::excluded_volume_force(distance, given_variables[1]);
}

double FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_potential(double distance, double* given_variables)
{
  return FORCE::NAPLE::excluded_volume_potential(distance, given_variables[1]);
}

long FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond)
{
  // given_POT.force_variables = new double [3];
  given_POT.force_variables = (double*) mkl_malloc(3*sizeof(double), BIT);
  given_POT.force_variables[0] = atof(given_cond("repulsion_coefficient").c_str());
  given_POT.force_variables[1] = atof(given_cond("effective_distance").c_str());
  given_POT.force_variables[2] = 1./sqrt(given_POT.force_variables[0]);
  // cout << "======== NAPLE_repulsion_set =========\n";
  // cout << given_POT.force_variables[0] << '\t' << given_POT.force_variables[1] << '\t' << given_POT.force_variables[2] << endl;

  given_POT.f_repulsion = MAP_excluded_volume_force;
  given_POT.e_repulsion = MAP_excluded_volume_potential;

  given_POT.f_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_connector = FORCE::DEFAULT::EMPTY_force_contribution;

  given_POT.scale_random = MAP_time_scaling_random;

  return 0;
}


double FORCE::GAUSSIAN::spring_force(double distance, double N_dimension)
{
  return N_dimension*distance;
}

double FORCE::GAUSSIAN::spring_potential(double distance, double N_dimension)
{
  return (0.5)*N_dimension*distance*distance;
}

double FORCE::FENE::non_Gaussian_factor(double distance, double N_dimension, double ratio_RM_R0)
{
  return 1.0/(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
}

double FORCE::FENE::spring_force(double distance, double N_dimension, double ratio_RM_R0)
{
  return non_Gaussian_factor(distance, N_dimension, ratio_RM_R0)*FORCE::GAUSSIAN::spring_force(distance, N_dimension);
}

double FORCE::FENE::spring_potential(double distance, double N_dimension, double ratio_RM_R0)
{
  return (-(double)N_dimension/2.0)*pow(ratio_RM_R0, 2.0)*log(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
}


double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force(double distance, double* given_variables)
{
  return FORCE::GAUSSIAN::spring_force(distance, given_variables[3]);
}

double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_potential(double distance, double* given_variables)
{
  return FORCE::GAUSSIAN::spring_potential(distance, given_variables[3]);
}

double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_force(double distance, double* given_variables)
{
  return FORCE::FENE::spring_force(distance, given_variables[3], given_variables[5]);
}

double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_potential(double distance, double* given_variables)
{
  return FORCE::FENE::spring_potential(distance, given_variables[3], given_variables[5]);
}


long FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond)
{
  given_POT.force_variables = new double [6];
  given_POT.force_variables[0] = atof(given_cond("repulsion_coefficient").c_str());
  given_POT.force_variables[1] = atof(given_cond("effective_distance").c_str());
  given_POT.force_variables[2] = 1./sqrt(given_POT.force_variables[0]);
  given_POT.force_variables[3] = atol(given_cond("N_dimension").c_str());
  given_POT.force_variables[4] = atof(given_cond("l_cap").c_str());
  // cout << "======== NAPLE_repulsion_set =========\n";
  // cout << given_POT.force_variables[0] << '\t' << given_POT.force_variables[1] << '\t' << given_POT.force_variables[2] << endl;
  given_POT.f_repulsion = FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_force;
  given_POT.e_repulsion = FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_potential;

  // given_POT.f_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;
  // given_POT.e_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;

  if(given_cond("FENE")=="TRUE")
    {
      given_POT.force_variables[5] = atof(given_cond("ratio_RM_R0").c_str());
      given_POT.f_connector = MAP_FENE_spring_force;
      given_POT.e_connector = MAP_FENE_spring_force;
    }
  else
    {
      given_POT.f_connector = MAP_Gaussian_spring_force;
      given_POT.e_connector = MAP_Gaussian_spring_potential;
    }

  given_POT.w_function = Detachment_weight;
  given_POT.scale_random = MAP_time_scaling_random;
  return 0;
}

double FORCE::NAPLE::MC_ASSOCIATION::Detachment_weight(double distance, double tension, double* given_variables)
{
  // given_variables[3] == Nd
  // given_variables[4] == l_cap
  return exp(tension*given_variables[4]);
}

