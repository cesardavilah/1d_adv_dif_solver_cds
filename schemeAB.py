import matplotlib.pyplot as plt
import numpy as np

import tdma_solver as tdma
import cds

if __name__ == '__main__':


    save_to = "G:\My Drive\PHD\Fall2023\Computational_environmental_fluid_mechanics\hw4\PLOTS"
    #CONSTANTS
    U = -4
    L = 5
    RHO = 900
    GAMMA = 650

    PE = (U * RHO * L) / GAMMA

    PHILEFTSLOPE = 2 * (PE/L) * (2/(1 - np.exp(-360/13)))
    PHIRIGHT = 1

    # NN = np.arange(2,80, 1)
    NN = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]
    # NN = np.zeros(100)
    # NN[0] = 10


    # def generate_doubling_list(N):
    #     if N <= 0:
    #         return []
    #
    #     result = []
    #     value = 10
    #
    #     for _ in range(N):
    #         result.append(value)
    #         value *= 2
    #
    #     return result
    #
    # NN = generate_doubling_list(15)

    # NUM_ER = 1.7
    # NUM_ER = 0.7
    # EERR = [NUM_ER, NUM_ER**(1/2), NUM_ER**(1/4), NUM_ER**(1/8), NUM_ER**(1/16)]
    UNIFORM_EERR = np.ones(len(NN))
    # EX_RATIO =
    FLAG_N = 0



    max_err_first = []
    max_err_second = []

    rel_err_first = []
    rel_err_second = []
    for i in range(0, len(NN)):
        N = NN[i]
        EX_RATIO = UNIFORM_EERR[i]

        if EX_RATIO == 1:
            spacing = "uniform"
        else:
            spacing = "non-uniform"

        #discretize the environment
        domain, delta_x = cds.discretize_domain(len_domain=L, expansion_ratio=EX_RATIO, number_of_elements=N, debug=False)

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
        # solution_tdma = tdma.tdma_solver(ap=apc, ae=aec, aw=awc, phi_left=PHILEFT, phi_right=PHIRIGHT)
        solution_tdma_first = tdma.tdma_solver_neumann_dirichlet_first_order(ap=apc, ae=aec, aw=awc, phi_left_slope=PHILEFTSLOPE, phi_right=PHIRIGHT, delta_x=delta_x)
        solution_tdma_second = tdma.tdma_solver_neumann_dirichlet_second_order(ap=apc, ae=aec, aw=awc, phi_left_slope=PHILEFTSLOPE, phi_right=PHIRIGHT, delta_x=delta_x)
        exact_solution = cds.exact_solution_general(A=4, B=1, domain=domain, PE=PE, L=L)



        # exact_solution = cds.exact_solution(phi_left=PHILEFT, phi_right=PHIRIGHT, u=U, L=L, rho=RHO, gamma=GAMMA, domain=domain)
        #
        # #Present solution
        # print(solution_tdma_first)
        print(exact_solution)


        # print(exact_solution)
        error_first = np.abs(solution_tdma_first - exact_solution)
        error_second = np.abs(solution_tdma_second - exact_solution)

        max_err_first.append(np.max(error_first))
        max_err_second.append(np.max(error_second))

        relative_error_first = np.abs(solution_tdma_first[0] - exact_solution[0]) / np.abs(exact_solution[0]) * 100

        relative_error_second = np.abs(solution_tdma_second[0] - exact_solution[0]) / np.abs(exact_solution[0]) * 100

        rel_err_first.append(relative_error_first)
        rel_err_second.append(relative_error_second)












        # # CONVERGENCE STUDIES
        plt.scatter(domain, solution_tdma_first, color="red", label="Numerical solution First Order", marker='x')
        plt.scatter(domain, solution_tdma_second, color="blue", label="Numerical solution Second Order", marker='x')
        plt.scatter(domain, exact_solution, color="green", label="Exact solution", marker="^")
        plt.legend()
        plt.title(f"c) Convergence studies in {spacing} spacing for N={N}")
        plt.savefig(rf"{save_to}\Convergence studies in {spacing} spacing for N={N}.png")
        plt.show()
        # # #

    # plt.clf()
    # plt.loglog(NN, max_err_first, label="Max error from Numerical solution First order")


    n1 = cds.power(NN, -1)
    n2 = cds.power(NN, -2)

    one_over_N = cds.oneover(NN)


    # plt.loglog(NN, n1, label="log(N) vs -1log(N)")
    plt.loglog(one_over_N, rel_err_first, label="Relative Error Numerical solution First order")
    plt.loglog(one_over_N, rel_err_second, label="Relative Error Numerical solution Second order")

    plt.loglog(one_over_N, n2, label="log(N) vs -2log(N)")

    plt.title(f"Relative % error vs 1/N values for uniform spacing")
    plt.legend()
    plt.savefig(rf"{save_to}\Relative_percentage_error")
    plt.show()






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
