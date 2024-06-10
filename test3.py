import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

gamma_factor=1.0

kappa=100
omega_1=20*np.pi*kappa
omega_2=20*np.pi*kappa
gamma_1=gamma_factor*kappa
gamma_2=-gamma_factor*kappa

# Define the system of differential equations
def system(t, y):
    a1, a2 = y
    da1_dt = ((omega_1-1j*gamma_1)*a1+kappa*a2)*(-1j)
    da2_dt = (kappa*a1+(omega_2-1j*gamma_2)*a2)*(-1j)
    return [da1_dt, da2_dt]

# Initial conditions
a1_0 = 1.0 + 0j
a2_0 = 0 + 0j
y0 = [a1_0, a2_0]

t_span = (0, 0.1)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# Solve the system of ODEs

sol = solve_ivp(system, t_span, y0, t_eval=t_eval)

# Extract the solutions
a1_sol = sol.y[0]
a2_sol = sol.y[1]

# Compute the amplitude (magnitude)
a1_amp = np.real(a1_sol)
a2_amp = np.real(a2_sol)

# Plot the results
plt.figure(figsize=(12, 6))

plt.rcParams['text.usetex'] = True


plt.plot(sol.t, a1_amp, label='a1(t)')
plt.plot(sol.t, a2_amp, label='a2(t)')
plt.title(r'Amplitude of $a_1(t)$ and $a_2(t)$ for $\gamma/\kappa=0.5$')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

plt.tight_layout()
plt.show()
