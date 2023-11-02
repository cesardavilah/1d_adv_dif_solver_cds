import numpy as np



def tdma_solver_neumann_dirichlet_second_order(ap: np.array, ae: np.array, aw: np.array, phi_left_slope, phi_right, delta_x, debug=False):
    N = len(ap)+2

    q = np.zeros(N)
    q_star = np.zeros(N)
    solution = np.zeros(N)



    # solution[0] =
    solution[N-1] = phi_right

    #do some padding

    apn = np.zeros(N)
    apn[1:N-1] = ap.copy()


    aen = np.zeros(N)
    aen[1:N-1] = ae.copy()


    awn = np.zeros(N)
    awn[1:N-1] = aw.copy()



    #SETUP
    q[0] = 2 * aen[1] * delta_x * phi_left_slope + q[1]
    q_star[0] = 2 * aen[1] * delta_x * phi_left_slope + q[1]

    awn[0] = 0
    aen[0] = (4 * aen[1]) + (apn[1])
    apn[0] = (awn[1]) - (3 * aen[1])

    #FORWARD ELIMINATION

    for i in range(1, N-1):
        apn[i] = apn[i] - ((awn[i] * aen[i - 1])/(apn[i - 1]))

    for i in range(1, N-1):
        q_star[i] = q[i] - ((awn[i] * q_star[i - 1])/(apn[i - 1]))

    #BACKWARD SUBSTITUTION
    for i in range(N-2, -1, -1):
        solution[i] = (q_star[i] - (aen[i] * solution[i+1]))/(apn[i])





    if debug:
        print("DEBUGGINS INSIDE @tdma_solver")
        print(ap)
        print(ae)
        print(aw)

    return solution

def tdma_solver_neumann_dirichlet_first_order(ap: np.array, ae: np.array, aw: np.array, phi_left_slope, phi_right, delta_x, debug=False):
    N = len(ap)+2

    q = np.zeros(N)
    q_star = np.zeros(N)
    solution = np.zeros(N)

    q[0] = phi_left_slope * delta_x
    q_star[0] = phi_left_slope * delta_x

    # solution[0] =
    solution[N-1] = phi_right

    #do some padding

    apn = np.zeros(N)
    apn[1:N-1] = ap.copy()
    apn[0] = -1

    aen = np.zeros(N)
    aen[1:N-1] = ae.copy()
    aen[0] = 1

    awn = np.zeros(N)
    awn[1:N-1] = aw.copy()
    awn[0] = 0

    #FORWARD ELIMINATION

    for i in range(1, N-1):
        apn[i] = apn[i] - ((awn[i] * aen[i - 1])/(apn[i - 1]))

    for i in range(1, N-1):
        q_star[i] = q[i] - ((awn[i] * q_star[i - 1])/(apn[i - 1]))

    #BACKWARD SUBSTITUTION
    for i in range(N-2, -1, -1):
        solution[i] = (q_star[i] - (aen[i] * solution[i+1]))/(apn[i])





    if debug:
        print("DEBUGGINS INSIDE @tdma_solver")
        print(ap)
        print(ae)
        print(aw)

    return solution

def tdma_solver(ap: np.array, ae: np.array, aw: np.array, phi_left, phi_right, debug=False):

    N = len(ap)+2

    q = np.zeros(N)
    q_star = np.zeros(N)
    solution = np.zeros(N)

    q[0] =  phi_left
    q_star[0] = phi_left

    solution[0] = phi_left
    solution[N-1] = phi_right

    #do some padding

    apn = np.zeros(N)
    apn[1:N-1] = ap.copy()
    apn[0] = 1

    aen = np.zeros(N)
    aen[1:N-1] = ae.copy()

    awn = np.zeros(N)
    awn[1:N-1] = aw.copy()

    #FORWARD ELIMINATION

    for i in range(1, N-1):
        apn[i] = apn[i] - ((awn[i] * aen[i - 1])/(apn[i - 1]))

    for i in range(1, N-1):
        q_star[i] = q[i] - ((awn[i] * q_star[i - 1])/(apn[i - 1]))

    #BACKWARD SUBSTITUTION
    for i in range(N-2, 0, -1):
        solution[i] = (q_star[i] - (aen[i] * solution[i+1]))/(apn[i])





    if debug:
        print("DEBUGGINS INSIDE @tdma_solver")
        print(ap)
        print(ae)
        print(aw)

    return solution