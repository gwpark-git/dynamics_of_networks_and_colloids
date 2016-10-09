// #include "read_file_condition.h"
#include "file_IO.h"
#include <string>
#include <iomanip> // for setw
// using namespace std;

bool
COND::
identify_effective_parameters(string given_argument)
{
  if (given_argument.length() == 0)
    return FALSE;
  // if (given_argument[0] == '/' && given_argument[1] == '/')
  // the comments in inp file will be changed from // to #
  if(given_argument[0] == '#')
    return FALSE;
  return TRUE;
}

COND::
COND
(char* fn)
  {
    GIVEN_FILE.open(fn);
    MKL_LONG cnt = 0;
    string line;
    string cond, val;
    while(getline(GIVEN_FILE, line))
      {
        cnt ++;
      }
    GIVEN_FILE.clear(); // since the previous get-line reach the EOF, this is bad-status. So, it need clear to use file object.
    GIVEN_FILE.seekg(0);
    N_arg = cnt; 

    arg = (string**) new string* [N_arg];
    for(MKL_LONG i=0; i<N_arg; i++)
      {
        arg[i] = (string*) new string [2];
      }

    cnt = 0;
    while(getline(GIVEN_FILE, line))
      {
        stringstream iss(line);
        getline(iss, cond, '='); arg[cnt][0] = cond;
        getline(iss, val, '\n'); arg[cnt][1] = val;
        cnt ++;
      }
    GIVEN_FILE.close();
    ERR = "ERR";
  }

string&
COND::
operator()
  (string option_type)
{
  for(MKL_LONG i=0; i<N_arg; i++)
    {
      // conditional phrase
      // if(arg[i][0] == option_type && (arg[i][0][0] != '/' && arg[i][0][1] != '/'))
      if(arg[i][0] == option_type && COND::identify_effective_parameters(arg[i][0]))
        {
          return arg[i][1];
        }
    }
  cout << "Bad condition " << option_type << endl;
  return ERR;
}

int
COND::
cond_print()
{
  if (arg)
    {
      cout << "### PRINTING ARGUMENTS ###\n";
      for(MKL_LONG i=0; i<N_arg; i++)
        {
          // if(arg[i][0][0] != '/' && arg[i][0][1] != '/')
          if(COND::identify_effective_parameters(arg[i][0]))
            cout << setw(50) << left << arg[i][0] << ": " << setw(50) << left  << arg[i][1] << endl;
        }
      cout << "### END PRINT ###\n";
    }
  return 0;
}


