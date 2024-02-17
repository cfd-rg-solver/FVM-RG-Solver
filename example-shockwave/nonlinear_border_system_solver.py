import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# consts

gamma = 1.667 # argon
argon_molarMass = 0.039948 # argon
argon_mass = 6.633521356992e-26 # argon

R = 8.3144598
kB = 1.38064852e-23

# set parameters
speed_of_sound_ar_room = 323 # m/s
Ma = 5
T_left = 300 # K
pressure = 100000 # Pa

# calculated parameters
velocity_left = Ma * speed_of_sound_ar_room
density_left = pressure*argon_molarMass / (R * T_left)


def solver(velocity_left, density_left, T_left):
    
    v_0, rho_0, T_0 = velocity_left, density_left, T_left
    
    # R = 8.3144598
    # kB = 1.38064852e-23
    # argon_molarMass = 0.039948
    # p_0 = R*T_0*rho_0/argon_molarMass

    # argon_mass = 6.633521356992e-26
    # E_0 = 3*kB*T_0/(2*argon_mass)

    def func(x):
        return [
            x[0] * x[1] - v_0 * rho_0, # x[0] - velocity, x[1] - density, x[2] - temperature
            x[1] * np.power(x[0],2) + R*x[2]*x[1]/argon_molarMass - rho_0 * np.power(v_0,2) - R*T_0*rho_0/argon_molarMass,
            x[1] * x[0] * (3*kB*x[2]/(2*argon_mass) + np.power(x[0],2)/2 + R*x[2]*x[1]/(argon_molarMass*x[1])) \
                 - rho_0 * v_0 * (3*kB*T_0/(2*argon_mass) + np.power(v_0,2)/2 + R*T_0*rho_0/(argon_molarMass*rho_0))
            ]
    
    vs, rhos, Ts = [], [], []
    for x in range(1, 800):
        cur_ans = fsolve(func, [x, x/100, x])
        if any(np.isclose(func(cur_ans), [0.0,0.0,0.0])):
            vs.append(cur_ans[0])
            rhos.append(cur_ans[1])
            Ts.append(cur_ans[2])
    ans = [np.median(vs), np.median(rhos), np.median(Ts)]
    
    print(func(ans)) # ����������� ���������� ������� � �������, �� ����� ������ ������� ������������ �����
    
    # plt.plot(vs)
    # plt.show()
    # plt.plot(rhos)
    # plt.show()
    # plt.plot(Ts)
    # plt.show()
    
    return ans

print("Initial conditions:")
print("v_0 = ", velocity_left)
print("rho_0 = ", density_left)
print("T_0 = ", T_left)

print("Getting answer via numeric scipy.fsolver:")
ans = solver(velocity_left, density_left, T_left)
print("v_n = ", ans[0])
print("rho_n = ", ans[1])
print("T_n = ", ans[2])

################################################################################################

def solver_approx(velocity_left, density_left, T_left):

    density_right = ((gamma + 1) * pow(Ma,2))/(2 + (gamma-1)*pow(Ma,2))*density_left
    velocity_right = velocity_left*density_left/density_right # �� �������

    pressure_left = R * T_left * density_left / argon_molarMass # �� ��������� ���������
    pressure_right = (pow(Ma,2)*2*gamma - (gamma-1))/(gamma+1)*pressure_left
    T_right = pressure_right/(density_right*R/argon_molarMass)

    return velocity_right, density_right, T_right

print("Getting answer via approximate solver:")
ans2 = solver_approx(velocity_left, density_left, T_left)
print("v_n = ", ans2[0])
print("rho_n = ", ans2[1])
print("T_n = ", ans2[2])
