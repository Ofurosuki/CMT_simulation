import numpy as np
import matplotlib.pyplot as plt

gamma_list=np.linspace(0,2,100)
eigenvalues_list=np.zeros((len(gamma_list),2),dtype=np.complex128)

for gamma_factor in gamma_list:
    kappa=100
    kappa_12=kappa
    kappa_21=kappa
    omega_1=20.005*np.pi*kappa
    omega_2=20*np.pi*kappa
    gamma_1=gamma_factor*kappa
    gamma_2=-gamma_factor*kappa

    eigen_values=np.linalg.eigvals([
        [omega_1-1j*gamma_1,kappa_12],
        [kappa_21,omega_2-1j*gamma_2]
    ])
    eigenvalues_list[gamma_list==gamma_factor]=eigen_values

plt.figure(figsize=(10, 5))

plt.subplot(121)
plt.plot(gamma_list, np.real(eigenvalues_list[:, 0]), label='Eigenvalue 1')
plt.plot(gamma_list, np.real(eigenvalues_list[:, 1]), label='Eigenvalue 2')
plt.xlabel('Gamma')
plt.ylabel('Real Part')
plt.legend()

plt.subplot(122)
plt.plot(gamma_list, np.imag(eigenvalues_list[:, 0]), label='Eigenvalue 1')
plt.plot(gamma_list, np.imag(eigenvalues_list[:, 1]), label='Eigenvalue 2')
plt.xlabel('Gamma')
plt.ylabel('Imaginary Part')

plt.tight_layout()
plt.show()


