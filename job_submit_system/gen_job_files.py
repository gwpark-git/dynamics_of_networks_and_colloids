
def COND(fn):
    given_condition = {}
    with open(fn, 'r') as f:
        while(1):
            try:
                cond_val = f.readline().rstrip('\n').split('=') # rstrip for preventing \n character into data
                given_condition[cond_val[0]] = cond_val[1]
            except:
                break
    return given_condition

def gen_jdl(given_condition):
    base_fn = given_condition['filename_base']
    with open(base_fn + '.jdl', 'w') as f:
        f.write('Type = "Job";\n')
        f.write('Executable = "%s";\n'%(base_fn + '.sh'))
        f.write('Arguments = "%s %s %s";\n'%(base_fn + '.inp', base_fn + '.tgz', int(given_condition['N_THREADS_BD'])))
        f.write('StdOutput = "%s";\n'%(base_fn + '.log'))
        f.write('StdError = "%s";\n'%(base_fn + '.log'))
        f.write('Requirement = RegExp("emi2-ce0[1-2].scope.unina.it:8443/cream-pbs-unina_hpc",other.GlueCEUniqueID);\n')
        f.write('SMPGranularity = %d;\n'%(int(given_condition['N_THREADS_BD'])))
        f.write('CpuNumber = %d;\n'%(int(given_condition['N_THREADS_BD'])))
        f.write('InputSandbox = {"%s", "%s", "%s", "%s"};\n'%('stochastic_HEUR.sh', 'stochastic_siulation', base_fn + '.inp', given_cond['CONTINUATION_TRAJ_FN']))
        f.write('OutputSandbox = {"%s"};\n'%(base_fn + '.log'))
