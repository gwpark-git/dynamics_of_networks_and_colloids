
#ifndef READ_FILE_CONDITION_H
#define READ_FILE_CONDITION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <mkl.h>
#include "H5Cpp.h"

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
  COND
  (char* fn);
  virtual
  ~COND()
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

  RECORD_DATA
  (COND& given_condition)
  {
    if(given_condition("output_path")=="FALSE")
      {
        filename_trajectory = (given_condition("filename_base") + ".traj").c_str();
        traj.open(filename_trajectory.c_str(), std::ios_base::app);
      
        filename_energy = (given_condition("filename_base") + ".ener").c_str();
        ener.open(filename_energy.c_str(), std::ios_base::app);
      
        filename_HASH = (given_condition("filename_base") + ".hash").c_str();
        hash.open(filename_HASH.c_str(), std::ios_base::app);
      
        filename_weight = (given_condition("filename_base") + ".weight").c_str();
        weight.open(filename_weight.c_str(), std::ios_base::app);

        if(given_condition("tracking_individual_chain")=="TRUE")
          {
            filename_chain = (given_condition("filename_base") + ".chain").c_str();
            chain.open(filename_chain.c_str(), std::ios_base::app);
          }

        if(given_condition("MC_LOG")=="TRUE")
          {
            filename_MC_LOG = (given_condition("filename_base") + ".MC_LOG").c_str();
            MC_LOG.open(filename_MC_LOG.c_str(), std::ios_base::app);
          }
	  
      }
    else
      {
        filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_base") + ".traj").c_str();
        traj.open(filename_trajectory.c_str(), std::ios_base::app);
      
        filename_energy = (given_condition("output_path") + '/' + given_condition("filename_base") + ".ener").c_str();
        ener.open(filename_energy.c_str(), std::ios_base::app);
      
        filename_HASH = (given_condition("output_path") + '/' + given_condition("filename_base") + ".hash").c_str();
        hash.open(filename_HASH.c_str(), std::ios_base::app);
      
        filename_weight = (given_condition("output_path") + '/' + given_condition("filename_base") + ".weight").c_str();
        weight.open(filename_weight.c_str(), std::ios_base::app);

        if(given_condition("tracking_individual_chain")=="TRUE")
          {
            filename_chain = (given_condition("output_path") + '/' + given_condition("filename_base") + ".chain").c_str();
            chain.open(filename_chain.c_str(), std::ios_base::app);
          }

        if(given_condition("MC_LOG")=="TRUE")
          {
            filename_MC_LOG = (given_condition("output_path") + '/' + given_condition("filename_base") + ".MC_LOG").c_str();
            MC_LOG.open(filename_MC_LOG.c_str(), std::ios_base::app);
          }
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




#endif
