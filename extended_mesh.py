import numpy as np
import pandas as pd
from porousZones.plot import plot_forces
from porousZones.analytical import compute_forces_xy, compute_force_z, alpha_v_const, inflow_uvw

# Define training data from HiFi CFD or Wind Tunnel
u_mag = 10
rho = 1
porous_dimensions = [0.2, 0.0842527, 0.25]

f_ij = np.array([[134.2093102, -125.64895069, 0], [-101.56379806, 95.26388569, 0], [0, 0, 1.04757228]])

psi = np.array([-80, -60, -40, -20, 1E-6, 20, 40, 60, 80])
theta = np.array([10, 30, 50, 70, 89.999999, 110, 130, 150, 170])
psi_c2 = np.rad2deg(alpha_v_const(f_ij[1, 1], f_ij[1, 0]))
print(psi_c2)

# HiFi data
fx_train = np.array([0.656813, 2.074982, 3.209272, 3.574992, 2.968339, 1.698684, 0.482158, -0.31434, -0.19178])
fy_train = np.array([-0.67436, -2.10182, -3.16394, -3.41387, -2.72193, -1.42857, -0.27297, 0.408222, 0.250917])
fz_train = np.array([0.097, 0.146, 0.137, 0.076, 0, -0.076, -0.137, -0.146, -0.097])


# Case 1
u0, v0 = inflow_uvw(u_mag, psi, 90)[0:2]
fx_a, fy_a = compute_forces_xy(porous_dimensions, rho, u0, v0, f_ij)

# Case 2
u0, v0, w0 = inflow_uvw(u_mag, psi_c2, theta)
fz_a = compute_force_z(porous_dimensions, rho, u0, v0, w0, f_ij)

# Read porous CFD data
porous_cfd = pd.read_csv("psi_sweep_results_extended_mesh.csv", comment='#')
porous_cfd_2 = pd.read_csv("theta_sweep_results_extended_mesh.csv", comment='#')


# Do some Plotting
fx_to_plot = {"HiFidelity CFD": [psi, fx_train], "2D Analytical": [psi, fx_a],
              "3D Porous CFD": [porous_cfd['Psi'], porous_cfd['Fx']]}
fy_to_plot = {"HiFidelity CFD": [psi, fy_train], "2D Analytical": [psi, fy_a],
              "3D Porous CFD": [porous_cfd['Psi'], porous_cfd['Fy']]}
fz_to_plot = {"HiFidelity CFD": [theta, fz_train], "2D Analytical": [theta, fz_a],
              "3D Porous CFD": [porous_cfd_2['Theta'], porous_cfd_2['Fz']]}

plot_forces("Fx", fx_to_plot)
plot_forces("Fy", fy_to_plot)
plot_forces("Fz", fz_to_plot)
