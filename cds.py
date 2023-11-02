import numpy as np
import matplotlib.pyplot as plt
import math


def power(my_list, pow):
    return [x**pow for x in my_list]

def oneover(my_list):
    return [1/x for x in my_list]
def exact_solution_general(A, B, domain, PE, L):
    return A*np.exp(domain * (PE/L)) + B


def exact_solution(phi_left, phi_right, u, L, rho, gamma, domain):

    pe = (rho * u * L) / (gamma)
    phi = np.zeros(len(domain))
    for i in range(0, len(domain)):
        phi[i] = ( phi_left + ((phi_right - phi_left) * ((math.exp(domain[i]*(pe/L)) - 1)/((math.exp(pe)) - 1))))

    return phi

def a_point_coefficients(GAMMA: float, domain: np.array, debug=False): #TODO
    '''
    :param GAMMA: This is the cosntant GAMMA, given.
    :param domain: This is the discretized domain as a list of xi points. This can be obtained from the function @discretize_domain
    :return: Returns a list of point coefficients from the second point python index [1] to the last point in the domain.
    '''
    if debug:
        print("DEBUGGING INSIDE @a_point_coefficients")
    #HELPER
    def a_point_at_idx(GAMMA: float, delta_xi: float, delta_xi_plus_one: float):
        p_coef = GAMMA * (2/(delta_xi * delta_xi_plus_one))
        return p_coef
    #HELPER

    apc = np.zeros(len(domain)-2)

    counter = 0
    for i in range(1, len(domain)-1):

        delta_xi = domain[i] - domain[i-1]
        delta_xi_plus_one = domain[i+1] - domain[i]
        apc[counter] = a_point_at_idx(GAMMA=GAMMA, delta_xi=delta_xi, delta_xi_plus_one=delta_xi_plus_one)
        counter += 1

        if debug:
            print(f"This is the index: {i}")
            print(f"This is delta_xi: {delta_xi} from: {domain[i]} minus {domain[i-1]}")
            print(f"This is delta_xi_plus_one: {delta_xi_plus_one} from: {domain[i+1]} minus {domain[i]}")


    if debug:
        print(domain)
        print(apc)






    return apc

def a_east_coefficients(RHO: float, U: float, GAMMA: float, domain: np.array, debug=False): #TODO
    '''

    :param RHO: This is the constant RHO, given.
    :param U: This is the constant U, given.
    :param domain: This is the discretized domain as a list of xi points. This can be obtained from the function @discretize_domain
    :return: Returns a list of east coefficients from the second point python index [1] to the last point in the domain/
    '''
    if debug:
        print("DEBUGGING INSIDE @a_east_coefficients")
    #HELPER
    def a_east_at_idx(GAMMA: float, RHO: float, U: float,
                      xi_minus_one: float,
                      xi_plus_one: float,
                      delta_xi: float,
                      delta_xi_plus_one: float,
                      ):

        convective = (RHO * U) / (xi_plus_one - xi_minus_one)

        diffusive = -1*((GAMMA * 2) / (delta_xi_plus_one * (delta_xi + delta_xi_plus_one)))

        e_coef = convective + diffusive
        return e_coef, convective, diffusive
    #HELPER


    aec = np.zeros(len(domain) - 2)
    aec_convective = np.zeros(len(domain) - 2)
    aec_diffusive = np.zeros(len(domain) - 2)

    counter = 0
    for i in range(1, len(domain) - 1):
        xi_minus_one = domain[i-1]
        xi_plus_one = domain[i+1]
        delta_xi = domain[i] - domain[i - 1]
        delta_xi_plus_one = domain[i + 1] - domain[i]

        aec[counter], aec_convective[counter], aec_diffusive[counter]  = a_east_at_idx(
            GAMMA=GAMMA,
            RHO=RHO,
            U=U,
            xi_minus_one=xi_minus_one,
            xi_plus_one=xi_plus_one,
            delta_xi=delta_xi,
            delta_xi_plus_one=delta_xi_plus_one,
        )
        counter += 1

        if debug:
            print(f"This is the index: {i}")
            print(f"This is xi_minus_one: {domain[i - 1]}")
            print(f"This is xi_plus_one: {domain[i + 1]}")
            print(f"This is delta_xi: {delta_xi} from: {domain[i]} minus {domain[i - 1]}")
            print(f"This is delta_xi_plus_one: {delta_xi_plus_one} from: {domain[i + 1]} minus {domain[i]}")
    if debug:
        print(domain)
        print(aec)



    return aec, aec_convective, aec_diffusive

def a_west_coefficients(RHO: float, U: float, GAMMA: float, domain: np.array, debug=False): #TODO
    '''

    :param RHO: This is the constant RHO, given.
    :param U: This is the constant U, given.
    :param domain: This is the discretized domain as a list of xi points. This can be obtained from the function @discretize_domain
    :return: Returns a list of west coefficients from the second point python index [1] to the last point in the domain/
    '''
    if debug:
        print(print("DEBUGGING INSIDE @a_west_coefficients"))

    #HELPER
    def a_west_at_idx(GAMMA: float, RHO: float, U: float,
                      xi_minus_one: float,
                      xi_plus_one: float,
                      delta_xi: float,
                      delta_xi_plus_one: float,
                      ):
        convective = -1*((RHO * U) / (xi_plus_one - xi_minus_one))

        diffusive = -1*((GAMMA * 2) / (delta_xi * (delta_xi + delta_xi_plus_one)))

        w_coef = convective + diffusive
        return w_coef, convective, diffusive
    #HELPER

    awc = np.zeros(len(domain) - 2)
    awc_convective = np.zeros(len(domain) - 2)
    awc_diffusive = np.zeros(len(domain) - 2)

    counter = 0
    for i in range(1, len(domain) - 1):
        xi_minus_one = domain[i - 1]
        xi_plus_one = domain[i + 1]
        delta_xi = domain[i] - domain[i - 1]
        delta_xi_plus_one = domain[i + 1] - domain[i]

        awc[counter], awc_convective[counter], awc_diffusive[counter] = a_west_at_idx(
            GAMMA=GAMMA,
            RHO=RHO,
            U=U,
            xi_minus_one=xi_minus_one,
            xi_plus_one=xi_plus_one,
            delta_xi=delta_xi,
            delta_xi_plus_one=delta_xi_plus_one,
        )
        counter += 1

        if debug:
            print(f"This is the index: {i}")
            print(f"This is xi_minus_one: {domain[i - 1]}")
            print(f"This is xi_plus_one: {domain[i + 1]}")
            print(f"This is delta_xi: {delta_xi} from: {domain[i]} minus {domain[i - 1]}")
            print(f"This is delta_xi_plus_one: {delta_xi_plus_one} from: {domain[i + 1]} minus {domain[i]}")
    if debug:
        print(domain)
        print(awc)

    return awc, awc_convective, awc_diffusive


def initial_deltax(len_domain: float, expansion_ratio: float, number_of_elements: int):
    '''
    :param len_domain: The lenght of the domain in meters
    :param expansion_ratio: The expansion ratio re as defined in lecture 9. For re>1 more elements on left. For re<1 more elements on right.
    :param number_of_elements: The number of elements desired, aka N
    :return: Initial delta_x as a value. This can be used to calculate the rest of the points in the domain.
    '''

    if expansion_ratio == 1:
        deltax = len_domain/number_of_elements
    else:
        deltax = len_domain * ((expansion_ratio - 1)/((expansion_ratio**number_of_elements) - 1))


    return deltax

def discretize_domain(len_domain: float, expansion_ratio: float, number_of_elements: int, debug=False):
    '''
    :param len_domain:
    :param expansion_ration:
    :param number_of_elements:
    :param debug:
    :return: Returns the discretized domain as a 1d numpy array with N number of points.
    '''

    #determine initial deltax
    delta_x = initial_deltax(len_domain=len_domain, expansion_ratio=expansion_ratio, number_of_elements=number_of_elements)

    domain = np.empty(number_of_elements+1)

    for i in range(0, len(domain)):
        if i == 0:
            domain[i] = 0
        elif i == 1:
            domain[i] = 0 + delta_x
        else:
            domain[i] = domain[i-1] + ((expansion_ratio**(i-1)) * delta_x)




    if debug:
        print(f"This is deltax {delta_x}")
        for d in domain:
            print(d)
        plt.scatter(domain, np.ones(len(domain)))
        plt.show()


    return domain, delta_x