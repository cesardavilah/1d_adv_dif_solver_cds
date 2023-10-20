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
    PHIRIGHT = 1

    N = 10
    EX_RATIO = 1.1



    #discretize the environment
    domain = cds.discretize_domain(len_domain=L, expansion_ratio=EX_RATIO, number_of_elements=N, debug=True)

    #set up the variables aw ap and ae
    apc = cds.a_point_coefficients(GAMMA=GAMMA, domain=domain, debug=False)
    aec = cds.a_east_coefficients(RHO=RHO, U=U, GAMMA=GAMMA, domain=domain, debug=False)
    awc = cds.a_west_coefficients(RHO=RHO, U=U, GAMMA=GAMMA, domain=domain, debug=False)


    #set up the matrix
    # num_rows = N-1
    # num_cols = N-1
    # p_matrix = np.empty((num_rows, num_cols))
    # e_matrix = np.empty((num_rows, num_cols))
    # w_matrix = np.empty((num_rows, num_cols))
    p_matrix = np.diag(apc, 0)
    e_matrix = np.diag(aec[:-1], 1)
    w_matrix = np.diag(awc[1:], -1)

    b = np.zeros(N-1)
    b[0] = -1*awc[0]
    b[-1] = -1*aec[-1]

    matrix = p_matrix + e_matrix + w_matrix
    #TODO CHECK THE MATRIX MANUALLY



    print("END")



    #use tdma
    solution = tdma.tdma_solver(matrix=matrix)


    #present solution



