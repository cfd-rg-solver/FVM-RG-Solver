import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

# ARGON DATA
# gamma = 1.667 # argon
# molarMass = 0.039948 # argon
# mass = 6.633521356992e-26 # argon
# speed_of_sound_room = 322.6 # argon, m/s t T = 300 K and p = 6.66 Pa conditions (https://webbook.nist.gov/cgi/fluid.cgi?ID=C7440371&Action=Page)

# METHANE DATA
gamma = 1.3084 # calculated, approx for methane from NIST: 1.302, other data: https://www.mem50212.com/MDME/iTester/get-info/thermodynamics.html 
molarMass = 0.016043 # methane
mass = 2.663732314e-26 # methane
speed_of_sound_room = 450.06 # methane, m/s at T = 300 K and p = 100 Pa conditions (https://webbook.nist.gov/cgi/fluid.cgi?ID=C74828&Action=Page)

# consts
R = 8.3144598
kB = 1.38064852e-23
Nav = 6.02214129e23
hc = 6.62559e-34 * 2.99792458e8

Ma = 5
T_left = 300 # K
pressure = 100 # Pa

# calculated parameters
velocity_left = Ma * speed_of_sound_room
density_left = pressure*molarMass / (R * T_left)
om_e = np.array([302550, 158270, 315680, 136740])
ds = np.array([1, 2, 3, 3])
es = hc*om_e # e_0
e_0000 = sum(hc * (om_e*ds)/2)
# print("e_0000 = ", e_0000)
D_diss = 3668582.3189 # m^-1, converted from 438.86 kJ/mol https://www.weizmann.ac.il/oc/martin/tools/hartree.html
max_vibr_lvls =  [9, 17, 9, 20] # [4, 4, 4, 4]

possible_inds = []
for i1 in range(max_vibr_lvls[0]):
    for i2 in range(max_vibr_lvls[1]):
        for i3 in range(max_vibr_lvls[2]):
            for i4 in range(max_vibr_lvls[3]):

                e = (
                    om_e[0] * (i1 + ds[0] / 2.) + 
                    om_e[1] * (i2 + ds[1] / 2.) + 
                    om_e[2] * (i3 + ds[2] / 2.) + 
                    om_e[3] * (i4 + ds[3] / 2.)
                )    

                if e <= D_diss:
                    possible_inds.append([i1, i2, i3, i4])

# possible_inds = [[0,0,0,0]] # ground state case
################################################################################################

print("Initial conditions:")
print("v_0 = ", velocity_left)
print("rho_0 = ", density_left)
print("T_0 = ", T_left)
print("______________________________________")

################################################################################################

def solver_approx(velocity_left, density_left, T_left):
    # Rankine-Hugoniot boundary conditions
    density_right = ((gamma + 1) * pow(Ma,2))/(2 + (gamma-1)*pow(Ma,2))*density_left
    velocity_right = velocity_left*density_left/density_right 

    pressure_left = R * T_left * density_left / molarMass 
    pressure_right = (pow(Ma,2)*2*gamma - (gamma-1))/(gamma+1)*pressure_left
    T_right = pressure_right/(density_right*R/molarMass)

    return velocity_right, density_right, T_right

print("Getting answer via approximate solver:")
ans2 = solver_approx(velocity_left, density_left, T_left)
print("v_n = ", ans2[0])
print("rho_n = ", ans2[1])
print("T_n = ", ans2[2])
print("______________________________________")


################################################################################################


def solver(velocity_left, density_left, T_left):
    
    v_0, rho_0, T_0 = velocity_left, density_left, T_left
    
    def func(x):
        # Energy calculated for polyatomic gas:
        # U_tr = 3 / 2 * kB * temp * n / density, where n = Nav * density / molarMass
        # U_rot = kB * temp / mass
        # U_vibr = sum((e_0 + e_0000) / (Z) * exp(-e_0 / (kB * temp)))
        
        Zvibr_0, Zvibr_1 = 0, 0
        for inds in possible_inds:
            s = ((inds[1]+1)*(inds[2]+1)*(inds[2]+2)*(inds[3]+1)*(inds[3]+2)/4)
            e_0 = sum(es*inds)
            Zvibr_0 += s*np.exp(-e_0/(kB*T_0))
            Zvibr_1 += s*np.exp(-e_0/(kB*x[2]))
        
        Uvibr_0, Uvibr_1 = 0, 0
        for inds in possible_inds:
            s = ((inds[1]+1)*(inds[2]+1)*(inds[2]+2)*(inds[3]+1)*(inds[3]+2)/4)
            e_0 = sum(es*inds)
            Uvibr_0 += s * e_0 * np.exp(-e_0/(kB*T_0))
            Uvibr_1 += s * e_0 * np.exp(-e_0/(kB*x[2]))
        
        Uvibr_0 = Uvibr_0/Zvibr_0 + e_0000
        # Uvibr_0 = 0 # only ground state case
        Uvibr_1 = Uvibr_1/Zvibr_1 + e_0000
        # Uvibr_1 = 0 # only ground state case
            
        return [
            x[0] * x[1] - v_0 * rho_0, # x[0] - velocity, x[1] - density, x[2] - temperature
            x[1] * np.power(x[0],2) + R*x[2]*x[1]/molarMass - rho_0 * np.power(v_0,2) - R*T_0*rho_0/molarMass,
            x[1] * x[0] * ( 3/2*kB*x[2]*Nav/molarMass +  3/2*kB*x[2]/mass + Uvibr_1/mass + np.power(x[0],2)/2 + R*x[2]*x[1]/(molarMass*x[1]))
            - rho_0 * v_0 * ( 3/2*kB*T_0*Nav/molarMass + 3/2*kB*T_0/mass + Uvibr_0/mass + np.power(v_0,2)/2 + R*T_0*rho_0/(molarMass*rho_0))
            ]
    
    vs, rhos, Ts = [], [], []
    
    for x in tqdm(range(1, 200)):
        cur_ans = fsolve(func, [x, x/100, x])
        if any(np.isclose(func(cur_ans), [0.0,0.0,0.0])):
            vs.append(cur_ans[0])
            rhos.append(cur_ans[1])
            Ts.append(cur_ans[2])
    ans = [np.median(vs), np.median(rhos), np.median(Ts)]
    
    print("Must be near zero values =", func(ans))

    return ans


print("Getting answer via numeric scipy.fsolver:")
ans = solver(velocity_left, density_left, T_left)
print("v_n = ", ans[0])
print("rho_n = ", ans[1])
print("T_n = ", ans[2])
print("______________________________________")


# Mach 3.8

# Initial conditions:
# v_0 =  1710.2279999999998
# rho_0 =  0.0006431766819856015
# T_0 =  300

# Post-shock conditions:
# v_n =  263.5241
# rho_n =  0.004174111
# T_n =  781.843

# Mach 3

# Initial conditions:
# v_0 =  1350.18
# rho_0 =  0.0006431766819856015
# T_0 =  300

# Post-shock conditions:
# v_n =  271.0525325102617
# rho_n =  0.0032038228325672277
# T_n =  624.6138539532955

# Mach 5

# Initial conditions:
# v_0 =  2250.3
# rho_0 =  0.0006431766819856015
# T_0 =  300

# Post-shock conditions:
# v_n =  262.1277941014496
# rho_n =  0.0055215071429467605
# T_n =  1040.5303197231733