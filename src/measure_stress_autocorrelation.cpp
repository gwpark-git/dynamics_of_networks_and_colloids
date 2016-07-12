
#include <iostream>
#include "../lib/matrix.h"
#include "../lib/analysis.h" // related with analysis and post-processing
using namespace std;

int
help()
{
  cout << "USAGE of measuring correlation function\n";
  cout << "argv[1] == given ener file\n";
  cout << "argv[2] == output file\n";
  return 0;
}

int
main(int argc, char* argv[])
{
  if(argc==1)
    {
      help();
      return 0;
    }
  MATRIX ener_data;
  ener_data.read_exist_data(argv[1]);
  printf("READ_FILE_DONE\n");
  MATRIX result(ener_data.rows/2, 7, 0.);
  for(MKL_LONG t=0; t<result.rows; t++)
    result(t, 0) = ener_data(t, 0);
  MKL_LONG N_THREADS = 6;

  using namespace POST_PROCESSING;
  compute_autocorrelation_brute_force_OMP(ener_data, 21, N_THREADS, result, 1);
  printf("PROCESS_21\n");
  compute_autocorrelation_brute_force_OMP(ener_data, 22, N_THREADS, result, 2);
  printf("PROCESS_22\n");      
  compute_autocorrelation_brute_force_OMP(ener_data, 23, N_THREADS, result, 3);
  printf("PROCESS_23\n");      
  compute_autocorrelation_brute_force_OMP(ener_data, 27, N_THREADS, result, 4);
  printf("PROCESS_27\n");
  compute_autocorrelation_brute_force_OMP(ener_data, 28, N_THREADS, result, 5);
  printf("PROCESS_28\n");
  compute_autocorrelation_brute_force_OMP(ener_data, 29, N_THREADS, result, 6);
  printf("PROCESS_29\n");      

  ofstream FILE;
  FILE.open(argv[2]);
  result.fprint(FILE);
  FILE.close();
  return 0;
  return 0;
}

