import matplotlib.pyplot as plt
import numpy as np

import tdma_solver as tdma
import cds

if __name__ == '__main__':

    #CONSTANTS
    U = 4
    L = 5
    RHO = 900
    GAMMA = 650

    PHILEFT = 3
    PHIRIGHT = 2

    N = 100
    EX_RATIO = 0.9

    #discretize the environment
    domain = cds.discretize_domain(len_domain=L, expansion_ratio=EX_RATIO, number_of_elements=N, debug=False)

    #set up the variables ap, ae, and aw
    apc = cds.a_point_coefficients(GAMMA=GAMMA, domain=domain, debug=False)
    aec, aec_convective, aec_diffusive = cds.a_east_coefficients(RHO=RHO, U=U, GAMMA=GAMMA, domain=domain, debug=False)
    awc, awc_convective, awc_diffusive = cds.a_west_coefficients(RHO=RHO, U=U, GAMMA=GAMMA, domain=domain, debug=False)

    #SANITY CHECK
    #APCONV + AECONV + AWCONV = 0
    san_check_conv = aec_convective + awc_convective
    #APDIF + AEDIF + AWDIF = 0
    san_check_dif = apc + aec_diffusive + awc_diffusive

    #RHS of equation
    b = np.zeros(len(domain) - 2)

    #TDMA solver
    solution_tdma = tdma.tdma_solver(ap=apc, ae=aec, aw=awc, phi_left=PHILEFT, phi_right=PHIRIGHT, b=b)

    #Present solution
    print(solution_tdma)



    plt.scatter(domain, solution_tdma)
    plt.show()

    print("END")
