import matplotlib.pyplot as plt
import numpy as np

import tdma_solver as tdma
import cds

if __name__ == '__main__':


    save_to = "G:\My Drive\PHD\Fall2023\Computational_environmental_fluid_mechanics\hw3\PLOTS"
    #CONSTANTS
    U = 4
    L = 5
    RHO = 900
    GAMMA = 650

    PHILEFT = 3
    PHIRIGHT = 1

    # NN = np.arange(2,80, 1)
    NN = [10, 20, 40, 80, 160]
    NUM_ER = 1.7
    # NUM_ER = 0.7
    EERR = [NUM_ER, NUM_ER**(1/2), NUM_ER**(1/4), NUM_ER**(1/8), NUM_ER**(1/16)]
    UNIFORM_EERR = [1, 1, 1, 1, 1]
    # EX_RATIO =
    FLAG_N = 0




    for i in range(0, len(NN)):
        N = NN[i]
        EX_RATIO = EERR[i]

        if EX_RATIO == 1:
            spacing = "uniform"
        else:
            spacing = "non-uniform"

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
        solution_tdma = tdma.tdma_solver(ap=apc, ae=aec, aw=awc, phi_left=PHILEFT, phi_right=PHIRIGHT)
        exact_solution = cds.exact_solution(phi_left=PHILEFT, phi_right=PHIRIGHT, u=U, L=L, rho=RHO, gamma=GAMMA, domain=domain)

        #Present solution
        print(solution_tdma)
        print(exact_solution)

        # CONVERGENCE STUDIES
        plt.scatter(domain, solution_tdma, color="red", label="Numerical solution", marker='x')
        plt.scatter(domain, exact_solution, color="green", label="Exact solution", marker="^")
        plt.legend()
        plt.title(f"c) Convergence studies in {spacing} spacing for N={N}")
        # plt.savefig(rf"{save_to}\Convergence studies in {spacing} spacing for N={N}.png")
        plt.show()



        #


        # ABSOLUTE ERROR
    #     max_error = np.abs(solution_tdma - exact_solution)
    #     max_error_overall = max(max_error)
    #     print(max_error_overall)
    #
    #     plt.plot(domain, max_error, label=f"N={N}")
    #
    #     plt.title(f"Max error in {spacing} spacing for different N")
    #     plt.xlabel("X")
    #     plt.ylabel("Error")
    #
    #
    #     if max_error_overall <= 10e-3:
    #         print("Error less than 10e-3")
    #         FLAG_N = N
    #         # break
    #     else:
    #         print("Error higher than 10e-3")
    #
    # plt.plot(domain, np.ones(len(domain)) * 10e-3, color='red', label="Error 10e-3", linestyle='dashed')
    # plt.legend()
    # # plt.savefig(rf"{save_to}\Max error in {spacing} spacing.png")
    # plt.show()
    # print(f"flag n {FLAG_N}")
    #
    #
    # #RATE OF CONVERGENESS
    # #
    # #
    # #
