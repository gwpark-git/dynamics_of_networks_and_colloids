
from numpy import *
import sys

def sign(x):
    if x < 0.:
        return -1.
    return 1.

if size(sys.argv) < 6:
    print 'USAGE: Converting Trajectory without PBC'
    print 'Note: this is purpose to measure dynamic qunatities'
    print 'argv[1] == input trajectory'
    print 'argv[2] == output data'
    print 'argv[3] == spatial dimension'
    print 'argv[4] == number of particles'
    print 'argv[5] == box dimension'
    print 'argv[6] == (OPTIONAL) Wi_tau_C: tau_C is characteristic time to update Langevin Equation'
    print 'argv[7] == (OPTIONAL: FELLOW argv[6]) shear_direction (if Wi > 0)'
    print 'argv[8] == (OPTIONAL: FELLOW argv[6]) shear_grad_direction (if Wi > 0)'
    print 'argv[9] == (OPTIONAL) remove history (TRUE/FALSE), default is FALSE'
    print 'argv[10] == (OPTIONAL) tracking single particle trajectory (TRUE/FALSE), default is TRUE'
else:
    dat = loadtxt(sys.argv[1])
    ND = int(sys.argv[3])
    NP = long(sys.argv[4])
    LB = float(sys.argv[5])
    Wi_tau_C = 0
    shear_direction = 0
    shear_grad_direction = 1
    try:
        dt = dat[1, 0] - dat[0, 0]
        Wi_tau_C = float(sys.argv[6])
        shear_direction = int(sys.argv[7])
        shear_grad_direction = int(sys.argv[8])
    except:
        Wi_tau_C = 0

    remove_history = 0
    try:
        if (sys.argv[9] == 'TRUE'):
            remove_history == 1
    except:
        remove_history = 0

    single_particle_trajectory = 1
    try:
        if (sys.argv[10] == 'FALSE'):
            single_particle_trajectory = 0
    except:
        single_particle_trajectory = 1
        
    Nt = shape(dat)[0]
    
    # def inv_PBC(x_now, x_next, LB):
    #     dX = x_next - x_now
    #     if abs(dX) > 0.5*LB:
    #         return inv_PBC(x_now, x_next - sign(dX)*LB, LB)
    #     return x_next

    def inv_PBC(x_now, x_next, LB):
        dX = x_next - x_now
        if abs(dX) > 0.5*LB: 
            # it check the new trajectory jump the PBC boundary.
            # note that without jumping PBC boundary, if the abs(dX) is larger than half of box dimension,
            # it is used wrong time step
            return inv_PBC(x_now, x_next - sign(dX)*LB, LB)
        return x_next

    def inv_PBC_loop(x_now, x_next, LB):
        # it is verified to have the same functionality with inv_PBC.
        # note that the recursion form for inv_PBC is already verified in previous time for equilibrium simulation
        # It is of importance to aware that this inverse mapping cannot identify boundary jump
        # since it is affected by history for jumping
        # therefore, boundary jump should be identified in elsewhere
        dX = x_next - x_now
        # count = 0
        while (abs(dX) > 0.5*LB):
            x_next -= sign(dX)*LB
            dX = x_next - x_now
            # count += sign(dX) # note that count will sum over all the sign, it will prevent multiple counting
        return x_next#, count


    # temporal backup to be on the safe side
    # dat_record = zeros(shape(dat))
    # for i in range(shape(dat)[0]):
    #     for j in range(shape(dat)[1]):
    #         dat_record[i,j] = dat[i,j]

    
    if Wi_tau_C > 0 and remove_history == 1:
        # it will remove the history of shear flow into coordination change
        for i in range(NP):
            index_shear_direction = 2*ND*i + 1 + shear_direction
            index_shear_grad_direction = 2*ND*i + 1 + shear_grad_direction
            coord_mov_shear = 0
            for t in range(Nt):
                # the sequence have been changed
                # since the history remove is related with previous time
                # not about the current time
                dat[t, index_shear_direction] -= coord_mov_shear
                coord_mov_shear += Wi_tau_C*dat[t, index_shear_grad_direction]*dt
                
            
    for i in range(NP):
        print 'processing particle %d'%(i)
        for k in range(ND):
            # the sequence of processing is of importance
            # it affect to the identification of inverse PBC boundary condition
            # here, the first re-mapping is set with shear_grad_direction
            # therefore, we can re-mapping the coordination for shear_direction without any distracted affect
            # then, the shear_direction coordination is applied with inversion of PBC condtion.
            # in consequence, we have succesful re-mapping
            # note again, the inv_PBC functions are for equilibrium
            ind_k = (k + shear_grad_direction) %ND
            # ind_k = (k + shear_direction) % ND
            index_Rik = 2*ND*i + 1 + ind_k
            record_pre_coord = 0
            record_now_coord = dat[0, index_Rik]
            coord_shift_factor = 0
            inv_PBC_shift_factor = 0
            for t in range(1, Nt):
                # ident_t = Nt/10
                # if t%ident_t==0:
                #     print 'processing %d out of %d'%(int(t/ident_t), 10)
                record_pre_coord = record_now_coord
                record_now_coord = dat[t, index_Rik]
                
                rev_map_coord = inv_PBC_loop(dat[t-1, index_Rik], dat[t, index_Rik], LB)
                dat[t, index_Rik] = rev_map_coord
                if Wi_tau_C > 0 and ind_k == shear_grad_direction:
                    # in this part, the shifting in shear_direction (default: x-direction) due to boundary jump on shear_grad_direction (default: y-direction) will be removed
                    # in consequence of this reverse mapping procedure, the coordinate on shear_direction will NOT track the real particle, but track PARTICLE in PBC BOX (like simulation).
                    # It is of importance to aware that mechanical perturbation reveals the degenerated cases left-right PBC box condition in shear_grad_direction since it should have different velocity contribution
                    # Even if our simulation will take configurational Langevin equation, the velocity gradient part make distingushable between particles in its image
                    index_shear_direction = 2*ND*i + 1 + shear_direction
                    diff = record_now_coord - record_pre_coord
                    if abs(diff) > 0.5*LB:
                        inv_PBC_shift_factor = Wi_tau_C*LB*dat[t, 0] 
                        pre_coord = dat[t, index_shear_direction]
                        for tmp_t in range(t, Nt):
                            # it will remove all the shear history for the future
                            dat[tmp_t, index_shear_direction] -= sign(diff)*inv_PBC_shift_factor

                    if(single_particle_trajectory == 1):
                        # this second part is related with correction for coordination of shear_grad_direction in order to track the real particle case
                        # the correction factor stored in coord_shift_factor sum over all the history
                        
                        diff = rev_map_coord - record_now_coord #differece between transition of y-coordinate
                        if abs(diff) > 0.5*LB:
                            coord_shift_factor += sign(diff)*Wi_tau_C*LB*dt # record the history of the coordinate change due to tracking single particle
                        dat[t, index_shear_direction] += coord_shift_factor
                        
                        
    savetxt(sys.argv[2], dat)
