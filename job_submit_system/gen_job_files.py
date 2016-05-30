import sys
from numpy import *
from collections import OrderedDict 

def COND(fn):
    given_condition = OrderedDict() # it keep the given order of input file (note that the typical dictionary in Python will be automatically sorted based on internal hash table)
    with open(fn, 'r') as f:
        while(1):
            try:
                cond_val = f.readline().rstrip('\n').split('=') # rstrip for preventing \n character into data
                given_condition[cond_val[0]] = cond_val[1]
            except:
                break
    return given_condition

def write_COND(given_condition, output_file):
    with open(output_file, 'w') as f_COND:
        for key, value in given_condition.items():
            f_COND.write('%s=%s\n'%(key,value))
    return 0
        


def gen_jdl(given_condition, out_path, script_file, jdl_file_name, inp_file_name, arg_str):
    base_fn = given_condition['filename_base']
    with open(jdl_file_name, 'w') as f:
        f.write('Type = "Job";\n')
        f.write('Executable = "%s";\n'%(base_fn + '.sh'))
        f.write('Arguments = "%s %s %s";\n'%(inp_file_name, inp_file_name.rstrip('.inp') + '.tgz', int(given_condition['N_THREADS_BD'])))
        f.write('StdOutput = "%s";\n'%(base_fn + '.log'))
        f.write('StdError = "%s";\n'%(base_fn + '.log'))
        f.write('Requirement = RegExp("emi2-ce0[1-2].scope.unina.it:8443/cream-pbs-unina_hpc",other.GlueCEUniqueID);\n')
        f.write('SMPGranularity = %d;\n'%(int(given_condition['N_THREADS_BD'])))
        f.write('CpuNumber = %d;\n'%(int(given_condition['N_THREADS_BD'])))
        # f.write('InputSandbox = {"%s", "%s", "%s", "%s"};\n'%(script_file, 'stochastic_simulation', base_fn + '.inp', given_condition['CONTINUATION_TRAJ_FN']))
        f.write('InputSandbox = {"%s", "%s", "%s"'%(script_file, 'stochastic_simulation', base_fn + '.inp'))
        for i in range(size(arg_str)):
            f.write(', "%s"'%(arg_str[i]))
        f.write('};')
        f.write('OutputSandbox = {"%s"};\n'%(base_fn + '.log'))
        f.write('PerusalFileEnable = true;\n')
        f.write('PerusalTimeInterval = 1800;')

if __name__=="__main__":
    if(size(sys.argv)< 4):
        print 'USAGE:'
        print 'argv[1] == given input file'
        print 'argv[2] == number of runs'
        print 'argv[3] == given script file'
    else:
        file_inp = sys.argv[1]
        out_path = file_inp + '.job'
        N_runs = int(sys.argv[2])
        script_file = sys.argv[3]

        arg_str = [] # it will be passed for jdl file generation
        for i in range(4, size(sys.argv)):
            arg_str.append(sys.argv[i])

        given_condition = COND(file_inp)
        print '\n* JOB SUBMITTION:'
        basic_filename = given_condition['filename_base']
        for i in range(N_runs):
            if(N_runs <> 1):
                given_condition['basic_random_seed']=i*100
                given_condition['basic_random_seed_SS']=i*10
            given_condition['filename_base'] = basic_filename + '_%02d'%(i)
            jdl_file_name = out_path + '/' + basic_filename + '_%02d.jdl'%(i)
            inp_file_name = jdl_file_name.rstrip('.jdl') + '.inp'
            print '** TODO JOB_ID = %d: seeds=(%d, %d), NT=%d, %s'%(i, int(given_condition['basic_random_seed']), int(given_condition['basic_random_seed_SS']), int(given_condition['N_THREADS_BD']), jdl_file_name)
            write_COND(given_condition, inp_file_name)
            gen_jdl(given_condition, out_path, script_file, jdl_file_name, inp_file_name, arg_str)
            
