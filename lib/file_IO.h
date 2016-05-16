
#ifndef READ_FILE_CONDITION_H
#define READ_FILE_CONDITION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <mkl.h>
#include "H5Cpp.h"

/* #include <stdlib.h> */
using namespace std;


class COND
{
public:
  
  string** arg;
  string ERR;
  MKL_LONG N_arg;
  ifstream GIVEN_FILE;

  int cond_print();
  
  COND()
  {
    std::cout << "COND CLASS must be contructed with valid fileformat\n";
  }
  COND(char* fn);
  virtual ~COND()
  {
    if(arg)
      {
        for(MKL_LONG i=0; i<N_arg; i++)
          {
            delete[] arg[i];
          }
        delete[] arg;
      }
  }

 string& operator()(string option_type);
};


class RECORD_DATA
{
 public:
  string filename_trajectory, filename_energy, filename_HASH, filename_weight, filename_chain, filename_MC_LOG;
  ofstream traj, ener, hash, weight, chain, MC_LOG;

  
  RECORD_DATA()
    {
      std::cout << "ERR: Basic Constructor is not supported for RECORD_DATA class\n";
    }

  RECORD_DATA(COND& given_condition)
    {
      filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_base") + ".traj").c_str();
      traj.open(filename_trajectory, std::ios_base::app);
      
      filename_energy = (given_condition("output_path") + '/' + given_condition("filename_base") + ".ener").c_str();
      ener.open(filename_energy, std::ios_base::app);
      
      filename_HASH = (given_condition("output_path") + '/' + given_condition("filename_base") + ".hash").c_str();
      hash.open(filename_HASH, std::ios_base::app);
      
      filename_weight = (given_condition("output_path") + '/' + given_condition("filename_base") + ".weight").c_str();
      weight.open(filename_weight, std::ios_base::app);

      if(given_condition("tracking_individual_chain")=="TRUE")
        {
          filename_chain = (given_condition("output_path") + '/' + given_condition("filename_base") + ".chain").c_str();
          chain.open(filename_chain, std::ios_base::app);
        }

      if(given_condition("MC_LOG")=="TRUE")
        {
          filename_MC_LOG = (given_condition("output_path") + '/' + given_condition("filename_base") + ".MC_LOG").c_str();
          MC_LOG.open(filename_MC_LOG, std::ios_base::app);
        }
    }
  
  ~RECORD_DATA()
    {
      if(traj)
        traj.close();
      if(ener)
        ener.close();
      if(hash)
        hash.close();
      if(weight)
        weight.close();
      if(chain)
        chain.close();
      if(MC_LOG)
        MC_LOG.close();
    }
};



/* class HDF5_FILE */
/* { */
/*  public: */
/*   H5File file; */
/*   Group data; */
/*   DataSet traj; */
/*   DataSet adj_list, adj_weight, adj_token; */
  
/*   HDF5_FILE() */
/*     { */
/*       std::cout << "HDF5_FILE is not support basic constructor\n"; */
/*     } */

/*   HDF5_FILE(COND& given_condition) */
/*     { */
/*       file(given_condition("output_path") + '/' + given_condition("filename_base"), H5F_ACC_TRUNC); */
/*       // H5F_ACC_TRUNC is default file creation properties and default file access properties. */
/*       // [CK] it is necessary to check the exact meaning for this flag */
/*       group(file.createGroup("/Data")); */

/*       // this is for trajectory  */
/*       MKL_LONG data_traj_RANK = 3; // time, particle index, dimension */
/*       MKL_LONG data_dim[3]; */
/*       data_dim[0] = atoi(given_condition("Nt").c_str())/(atoi(given_condition("N_skip").c_str())); // number of time for output */
/*       data_dim[1] = atoi(given_condition("Np").c_str()); // number of particles */
/*       data_dim[2] = atoi(given_condition("N_dimension").c_str()); // spatial dimensionality */
/*       DataSpace space_traj(data_traj_RANK, data_dim); */
/*       // the dataset made without any chunk as default */
/*       // it means the dataset cannot be compressed at this moment */
/*       traj = DataSet(file.createDataSet("/Data/traj", PredType::NATIVE_DOUBLE, space_traj)); */

/*       // this is for connectivity */
/*       MKL_LONG data_adj_RANK = 3; // time, particle index, adjustable dimensionlity for connections */
/*       data_dim[2] = 5; // default number for connections */
/*       MKL_LONG data_max_dim[3]; // it making fore resizeable dimensionality */
/*       data_max_dim[0] = data_dim[0]; // inherit the number of times */
/*       data_max_dim[1] = data_dim[1]; // inherit the number of particles */
/*       data_max_dim[2] = 2*atoi(given_condition("N_chains_per_particle").c_str()) + atoi(given_condition("tolerance_allowing_connections")); // this is maximally possible number for the distingushable connections */
/*       DataSpace space_adj(data_adj_RANK, data_dim, data_max_dim); */
/*       adj_list = DataSet(file.createDataSet("/Data/adj_list", PredType::NATIVE_INT, space_adj)); */
/*       adj_weight = DataSet(file.createDataSet("/Data/adj_weight", PredType::NATIVE_INT, space_adj)); */

/*       // this is for token */
/*       MKL_LONG data_adj_token_RANK = 2; */
/*       MKL_LONG data_adj_token_dim[2]; */
/*       data_adj_token_dim[0] = data_dim[0]; // inherit the number of times */
/*       data_adj_token_dim[1] = data_dim[1]; // inherit the number of particles */
/*       DataSpace space_adj_token(data_adj_token_RANK, data_adj_token_dim); */
/*       adj_token = DataSet(file.createDataSet("/Data/adj_token", PredType::NATIVE_INT, space_adj_token)); */

/*       // Below DataSpace will not be used during computation. */
/*       delete space_traj; */
/*       delete space_adj; */
/*       delete space_adj_token; */
/*     } */

/*   ~HDF5_FILE() */
/*     { */
/*       delete traj; */
/*       delete adj_list; */
/*       delete adj_weight; */
/*       delete adj_token; */
/*       delete traj; */
/*       delete group; */
/*       file.close(); */
/*     } */
/* }; */


#endif
