
#ifndef READ_FILE_CONDITION_H
#define READ_FILE_CONDITION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <mkl.h>
/* #include <stdlib.h> */
using namespace std;


class COND
{
public:
  /* MKL_LONG N_dimension, Np, Nt; */
  /* double dt; */
  
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

#endif
