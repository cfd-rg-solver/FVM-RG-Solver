import math

# Constants for Ar
gamma = 1.67
kb = 1.38064852e-23
T = 296.65
L = 1
R_ar = 8.314 / 0.039948
rho = 0.00000187 #0.000115 # 0.00000478
v = 1377.77

# Pressure calculation
P = rho * R_ar * T

# Reynolds number calculation
mu = 1.5755e-05 # 2.0988e-05
nu = mu / rho
Re = v * L / nu

# Mach number calculation
a = math.sqrt(P * gamma / rho)
M = v / a

# Knudsen number calculations
# v1
# sigma = 71e-12  # WRONG !!!
# Kn_v1 = kb * T / (math.sqrt(2) * math.pi * sigma**2 * P * L)

# # v2
# Kn_v2 = M / Re * math.sqrt(gamma * math.pi / 2)

# v3 (MFP / L)
Kn_v3 = mu / P * math.sqrt(math.pi * R_ar * T / 2) / L

# Print results
print(f"Pressure (P): {P}")
print(f"Reynolds number (Re): {Re}")
print(f"Mach number (M): {M}")
# print(f"Knudsen number v1 (Kn_v1): {Kn_v1}")
# print(f"Knudsen number v2 (Kn_v2): {Kn_v2}")
print(f"Knudsen number v3 (Kn_v3): {Kn_v3}")