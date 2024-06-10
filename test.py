import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the parameters
omega_1 = 1000
omega_2 = 1000
k_12 = 0.001
k_21 = 0.001
gamma_1 = 1.001
gamma_2 = 0.999

# Define the system of differential equations
def system(t, y):
    a1, a2 = y
    da1_dt = -1j * omega_1 * a1 + k_12 * a2+ gamma_1 * a1
    da2_dt = -1j * omega_2 * a2 + k_21 * a1 + gamma_2 * a2
    return [da1_dt, da2_dt]

# Initial conditions
a1_0 = 1.0 + 0j
a2_0 = 0.0 + 0j
y0 = [a1_0, a2_0]

# Time span
t_span = (0, 100)
t_eval = np.linspace(t_span[0], t_span[1], 4000)

# Solve the system of ODEs
sol = solve_ivp(system, t_span, y0, t_eval=t_eval)

# Extract the solutions
a1_sol = sol.y[0]
a2_sol = sol.y[1]

# Compute the amplitude (magnitude)
a1_amp = np.abs(a1_sol)
a2_amp = np.abs(a2_sol)

# Plot the results
plt.figure(figsize=(12, 6))

plt.plot(sol.t, a1_amp, label='|a1(t)|')
plt.plot(sol.t, a2_amp, label='|a2(t)|')
plt.title('Amplitude of a1(t) and a2(t)')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

plt.tight_layout()
plt.show()