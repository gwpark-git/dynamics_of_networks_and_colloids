#include "read_file_condition.h"
#include <string>
// using namespace std;

COND::COND(char* fn)
  {
    GIVEN_FILE.open(fn);
    long cnt = 0;
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
    for(long i=0; i<N_arg; i++)
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

string& COND::operator()(string option_type)
{
  for(long i=0; i<N_arg; i++)
    {
      if(arg[i][0] == option_type)
        {
          return arg[i][1];
        }
    }
  cout << "Bad condtion " << option_type << endl;
  return ERR;
  // string out = "ERR";
  // return out;
}

int COND::cond_print()
{
  if (arg)
    {
      cout << "### PRINTING ARGUMENTS ###\n";
      for(long i=0; i<N_arg; i++)
        {
          cout << arg[i][0] << '\t' << arg[i][1] << endl;
        }
      cout << "### END PRINT ###\n";
    }
  return 0;
}


// long main(long argc, char* argv[])
// {
//   // string fn(argv[1]);
//   COND given_condition(argv[1]);
//   cout << "###READ_TEST###\n";
//   // cout << type(given_condition.arg[0][0]) << endl;
//   // printf("%s\n", given_condition.arg[0][0]);
//   cout << given_condition.arg[0][0] << '\t' << given_condition.arg[0][1] << endl;
//   // cout << "dt = " << given_condition("dt") << endl;
//   // cout << "Np = " << given_condition("Np") << endl;
//   // string data_condition[6][2] = {};
//   // string cond, line, val;
//   // ifstream FILE(argv[1]);
//   // // ifstream FILE.open(argv[1], 'rt');
//   // long cnt = 0;
//   // while(getline(FILE,line))
//   //   {
//   //     stringstream iss(line);
//   //     getline(iss, cond, '=');
//   //     getline(iss, val, '\n');
//   //     data_condition[cnt][0] = cond;
//   //     data_condition[cnt][1] = val;
//   //     cnt++;
//   //     // cout << cond << '\t' << val << endl;
//   //   }
//   // FILE.close();


//   // for(long i=0; i<6; i++)
//   //   {
//   //     // cout << data_condition[i][0] << '\t' << data_condition[i][1] << endl;
//   //     cout << data_condition[i][1] << '\t' << atof(data_condition[i][1].c_str()) << endl;

//   //   }
//   return 0;
// }
