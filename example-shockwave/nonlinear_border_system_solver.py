import numpy as np
from scipy.optimize import fsolve


def solver(left_params):
    
    v_0, rho_0, p_0, E_0 = left_params

    def func(x):
        return [
            x[0] * x[1] - v_0 * rho_0,
            x[1] * x[0]^2 + x[2] - rho_0 * v_0^2 - p_0,
            x[1] * x[0] * (x[3] + x[0]^2/2 + x[2]/x[1]) - rho_0 * v_0 * (E_0 + v_0^2/2 + p_0/rho_0)
            ]

    ans = fsolve(func, [1, 1, 1, 1])
    return ans