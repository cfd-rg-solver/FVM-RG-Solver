import numpy as np
from scipy.optimize import fsolve


velocity_left = 0
density_left = 0.800773

T_left = 1000
UniversalGasConstant = 8.3144598
kB=1.38064852e-23
argon_molarMass = 0.039948
pressure_left = UniversalGasConstant * T_left * density_left / argon_molarMass

argon_mass = 6.633521356992e-26
energy_left = 3. * kB * T_left / (2. * argon_mass)


def solver(velocity_left, density_left, T_left):
    
    v_0, rho_0, T_0 = velocity_left, density_left, T_left
    
    R = 8.3144598
    kB = 1.38064852e-23
    argon_molarMass = 0.039948
    # p_0 = R*T_0*rho_0/argon_molarMass

    argon_mass = 6.633521356992e-26
    # E_0 = 3*kB*T_0/(2*argon_mass)

    def func(x):
        return [
            x[0] * x[1] - v_0 * rho_0,
            x[1] * np.power(x[0],2) + R*x[2]*x[1]/argon_molarMass - rho_0 * np.power(v_0,2) - R*T_0*rho_0/argon_molarMass,
            x[1] * x[0] * (3*kB*x[2]/(2*argon_mass) + np.power(x[0],2)/2 + R*x[2]*x[1]/(argon_molarMass*x[1])) - rho_0 * v_0 * (3*kB*T_0/(2*argon_mass) + np.power(v_0,2)/2 + R*T_0*rho_0/(argon_molarMass*rho_0))
            ]

    ans = fsolve(func, [1, 1, 1])
    return ans

ans = solver(velocity_left, density_left, T_left)
print("v_n = ", ans[0])
print("rho_n = ", ans[1])
print("T_n = ", ans[2])