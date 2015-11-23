#include "lib_traj.h"

using namespace std;

long TRAJECTORY::read_exist_traj(const char* fn_given_traj)
{
  ifstream GIVEN_FILE;
  GIVEN_FILE.open(fn_given_traj);
  long cnt = 0;
  string line;
  while(getline(GIVEN_FILE, line))
    {
      cnt ++;
    }
  GIVEN_FILE.clear();
  GIVEN_FILE.seekg(0);
  for(long i=0; i<cnt-1; i++)
    {
      getline(GIVEN_FILE, line);
    }
  // then, flag has set for the last line
  GIVEN_FILE >> (*this)(0); // time recording
  for(long i=0; i<Np; i++)
    {
      for(long i_RV=0; i_RV<=1; i_RV++)
        {
          for(long k=0; k<dimension; k++)
            {
              GIVEN_FILE >> (*this)(i_RV, 0, i, k); // modified version for all the record
            }
        }
    }
  for(long i=0; i<Np; i++)
    {
      for(long i_RV=0; i_RV<=1; i_RV++)
        {
          for(long k=0; k<dimension; k++)
            {
              cout << (*this)(i_RV, 0, i, k) << '\t';
            }
        }
      cout << endl;
    }
  GIVEN_FILE.close();
  return 0;
}

TRAJECTORY::TRAJECTORY() : MATRIX(N_t_test, 2*2*N_p_test + 1) 
{
  std::cout << "## TRAJECTORY class must used with input arguments such as total number of time and number of columns." << std::endl;
  std::cout << "## At the moment, it is temporally called by test functionality with Nt = " << N_t_test << " and Np = " << N_p_test << std::endl;
  std::cout << "As default, dimension is set to 2d" << std::endl;
  dimension = 2;
  initialization(N_t_test, N_p_test, 0.001);
}

// WARNING: Notice that this copy-constructor is not verified with the inherance of MATRIX class
TRAJECTORY::TRAJECTORY(const TRAJECTORY& TRAJ)// this is copy constructor
{
  initialization(TRAJ.Nt, TRAJ.Np, TRAJ.dt);
  for(long i=0; i<size; i++)
    {
      data[i] = TRAJ.data[i];
    }
}

TRAJECTORY::TRAJECTORY(COND& given_condition, long N_basic) : MATRIX(N_basic, atol(given_condition("N_dimension").c_str())*2*atol(given_condition("Np").c_str()) + 1, 0.)
{
  initialization_COND(given_condition);
}


long TRAJECTORY::initialization(long N_time, long N_particle, double given_dt)
{
  Nt = N_time;
  c_t = 0;
  Np = N_particle;
  box_dimension = (double*)mkl_malloc(dimension*sizeof(double), BIT);
  for (long k=0; k<dimension; k++)
    {
      box_dimension[k] = 10.0;
    }
  dt = given_dt;
  srandom(0);

  return 0;
}

long TRAJECTORY::initialization_COND(COND& given_condition)
{
  Method = given_condition("Method");
  if(given_condition("Integrator") == "Euler")
    Nt = 2;
  c_t = 0;
  Np = atol(given_condition("Np").c_str());
  dimension = atol(given_condition("N_dimension").c_str());
  box_dimension = (double*)mkl_malloc(dimension*sizeof(double), BIT);
  for (long k=0; k<dimension; k++)
    {
      box_dimension[k] = atof(given_condition("box_dimension").c_str());
    }
  dt = atof(given_condition("dt").c_str());

  N_energy_frequency = atol(given_condition("N_energy_frequency").c_str());
  srandom(0);

  if (given_condition("CONTINUATION_TRAJ")=="TRUE")
    {
      read_exist_traj(given_condition("CONTINUATION_TRAJ_FN").c_str());
    }
  else
    {
      GENERATOR::random_position_generator(*this);
    }
  return 0;
}


double& TRAJECTORY::operator()(long time_t)
{
  // it will return the reference of time
  return data[index(time_t, 0)];
}

double& TRAJECTORY::operator()(long time_t, long bead_i, long dimension_k)
{
  long index_position = 2*dimension*bead_i + 1 + dimension_k;
  // printf("INDEX(%ld, %ld, %ld) = %ld", time_t, bead_i, dimension_k, index_position);
  return data[index(time_t, index_position)];
}

double& TRAJECTORY::operator()(long i_RV, long time_t, long bead_i, long dimension_k)
{
  long index_position = 2*dimension*bead_i + 1 + dimension_k + 2*i_RV;
  return data[index(time_t, index_position)];
}



long GENERATOR::random_position_generator_REF(TRAJECTORY& TRAJ, MATRIX& R_VEC_TRANS)
{
  if (R_VEC_TRANS.rows != TRAJ.Np || R_VEC_TRANS.cols != TRAJ.dimension)
    {
      R_VEC_TRANS.initial(TRAJ.Np, TRAJ.dimension, 0.);
    }
  RANDOM::random_vector_generator(R_VEC_TRANS);
  for (long k=0; k<TRAJ.dimension; k++)
    {
      for (long i=0; i<TRAJ.Np; i++)
        {
          TRAJ(0, i, k) = (R_VEC_TRANS(i, k) + 0.5)*TRAJ.box_dimension[k];
        }
    }
  return 0;
}


long GENERATOR::random_position_generator(TRAJECTORY& TRAJ)
{
  MATRIX R_VEC_TRANS(TRAJ.Np, TRAJ.dimension, 0.);
  GENERATOR::random_position_generator_REF(TRAJ, R_VEC_TRANS);
  return 0;
}

