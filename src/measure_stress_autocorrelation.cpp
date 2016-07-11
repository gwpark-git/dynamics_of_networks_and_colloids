
#include <iostream>
#include "../lib/matrix.h"
#include "../lib/analysis.h"

int
main(int argc, char* argv[])
{
  MATRIX data;
  data.read_exist_data(argv[1]);
  for(MKL_LONG i=0; i<data.rows; i++)
    {
      for(MKL_LONG j=0; j<data.cols; j++)
        {
          printf("%3.2lf\t", data(i,j));
        }
      printf("\n");
    }

  return 0;
}

