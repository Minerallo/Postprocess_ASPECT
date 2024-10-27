import numpy as np
import matplotlib.pyplot as plt



def compute_mu(mu_d, mu_s, dynamic_characteristic_strain_rate, X, strain_rate_dev_inv2):
    # Compute the value of mu using the formula provided
    term = (1 + (strain_rate_dev_inv2 / dynamic_characteristic_strain_rate)) ** X
    mu = mu_d + (mu_s - mu_d) / term
    return mu

# Parameters
mu_d = 0.045
mu_s = 0.30
dynamic_characteristic_strain_rate = 2e-16
X = 2

angle_in_degrees = np.degrees(np.arctan(0.3))
print('Internal angle mu_d =', angle_in_degrees)

# Create a vector of strain rate deviations spanning several orders of magnitude
strain_rate_dev_inv2_values = np.logspace(-18, -14, 100)  # 100 points from 10^-18 to 10^-14

# Compute mu for each value in the strain rate deviation vector
mu_values = [compute_mu(mu_d, mu_s, dynamic_characteristic_strain_rate, X, strain_rate) for strain_rate in strain_rate_dev_inv2_values]

# Plotting
plt.figure(figsize=(10, 6))
plt.semilogx(strain_rate_dev_inv2_values, mu_values, marker='', linestyle='-', label='Mu vs. Strain Rate')
plt.xlabel('Strain Rate  (s-1)')
plt.ylabel('Mu Value')
plt.title('Mu Value vs. Strain Rate')
plt.grid(True)
plt.show()
