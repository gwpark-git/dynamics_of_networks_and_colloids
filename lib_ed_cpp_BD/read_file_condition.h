
#ifndef READ_FILE_CONDITION_H
#define READ_FILE_CONDITION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
/* #include <stdlib.h> */
using namespace std;


class COND
{
public:
  /* long N_dimension, Np, Nt; */
  /* double dt; */
  
  string** arg;
  string ERR;
  long N_arg;
  ifstream GIVEN_FILE;

  int cond_print();
  
  COND()
  {
    std::cout << "COND CLASS must be contructed with valid fileformat\n";
  }
  COND(char* fn);
  ~COND()
  {
    if(arg)
      {
        for(long i=0; i<N_arg; i++)
          {
            delete[] arg[i];
          }
        delete[] arg;
      }
  }

 string& operator()(string option_type);
};

#endif
