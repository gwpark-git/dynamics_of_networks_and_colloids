

#inlcude <iostream>
#include "../lib/read_file_condition.h"
#include "../lib/trajectory.h"

int
help()
{
  cout << "DESIGN post processing code\n";
  return 0;
}

MKL_LONG
check_line_number
(const char* fn_given_file)
{
  ifstream GIVEN_FILE;
  GIVEN_FILE.open(fn_given_file);
  MKL_LONG cnt = 0;
  string line;
  while(getline(GIVEN_FILE, line))
    cnt++;
  GIVEN_FILE.close();
  return cnt;
}

int
main
(int argc, char* argv[])
{
  if(argc==1)
    {
      help();
      return 0;
    }
  COND given_condition(argv[1]);

  string filename_trajectory = (given_condition("output_path") + '/' + given_condition("filename_base") + ".traj").c_str();
  string filename_chain = (given_condition("output_path") + '/' + given_condition("filename_base") + ".chain").c_str();
  MKL_LONG Nt = check_line_number(fn_given_file);
  
  return 0;
}
