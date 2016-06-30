
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
    print 'argv[6] == Wi_tau_C: tau_C is characteristic time to update Langevin Equation'
    print 'argv[7] == shear_direction (if Wi > 0)'
    print 'argv[8] == shear_grad_direction (if Wi > 0)'
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
    # inv_PBC_shift_factor_wo_dt = Wi_tau_C*LB;
    # tmp = box_dimension*int(2*inv_PBC_shift_factor_wo_dt)
    # ND = 2
    # NP = 100
    Nt = shape(dat)[0]
    # LB = 10.0

    # dat_conv = zeros([Nt, 2*ND*NP + 1])

    
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

    if Wi_tau_C > 0:
        # it will remove the history of shear flow into coordination change
        for i in range(NP):
            index_shear_direction = 2*ND*i + 1 + shear_direction
            index_shear_grad_direction = 2*ND*i + 1 + shear_grad_direction
            coord_mov_shear = 0
            for t in range(Nt):
                coord_mov_shear += Wi_tau_C*dat[t, index_shear_grad_direction]*dt
                dat[t, index_shear_direction] -= coord_mov_shear
                
            
    for i in range(NP):
        for k in range(ND):
            # the sequence of processing is of importance
            # it affect to the identification of inverse PBC boundary condition
            # here, the first re-mapping is set with shear_grad_direction
            # therefore, we can re-mapping the coordination for shear_direction without any distracted affect
            # then, the shear_direction coordination is applied with inversion of PBC condtion.
            # in consequence, we have succesful re-mapping
            # note again, the inv_PBC functions are for equilibrium
            ind_k = (k + shear_grad_direction) %ND
            index_Rik = 2*ND*i + 1 + ind_k
            record_pre_coord = 0
            record_now_coord = dat[0, index_Rik]
            coord_shift_factor = 0
            for t in range(1, Nt):
                record_pre_coord = record_now_coord
                record_now_coord = dat[t, index_Rik]
                
                rev_map_coord = inv_PBC_loop(dat[t-1, index_Rik], dat[t, index_Rik], LB)
                dat[t, index_Rik] = rev_map_coord

                if Wi_tau_C > 0 and ind_k == shear_grad_direction:
                    # when simple shear is implemented into the system
                    index_shear_direction = 2*ND*i + 1 + shear_direction
                    dat[t, index_shear_direction] -= coord_shift_factor
                    # temporal backup to be on the safe side
                    # diff = dat_record[t, index_Rik] - dat_record[t-1, index_Rik]
                    diff = record_now_coord - record_pre_coord
                    if abs(diff) > 0.5*LB:
                        inv_PBC_shift_factor = Wi_tau_C*LB*dat[t, 0] % LB
                        # coord_shift_factor = Wi_tau_C*record_now_coord*dt
                        pre_coord = dat[t, index_shear_direction]
                        for tmp_t in range(t, Nt):
                            # it will remove all the shear history for the future
                            dat[tmp_t, index_shear_direction] -= sign(diff)*inv_PBC_shift_factor + coord_shift_factor

                            
                    
                # if Wi_tau_C > 0 and ind_k == shear_direction:
                #     # it will remove the history of shear flow
                #     dat[t, index_Rik] -= coord_shift_factor 
            # if Wi_tau_C > 0:
            #     # remove history of simple shear
            #     for t in range(1, Nt):
            #         for tmp_t in range(t, Nt):
                        
    savetxt(sys.argv[2], dat)
