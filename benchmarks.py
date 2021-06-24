import numpy as np
import pandas as pd
from porousZones.plot import plot_forces
from porousZones.analytical import compute_forces_xy, compute_force_z, inflow_uvw, alpha_v_const

u_mag = 10
rho = 1
porous_dimensions = [0.5, 1, 1]

f_ij = np.array([[20, -18, 0], [-18, 16, 0], [0, 0, 2]])

psi = np.linspace(-89, 89, 179)
theta = np.linspace(1, 179, 179)
psi_c2 = np.rad2deg(alpha_v_const(f_ij[1, 1], f_ij[1, 0]))
print(psi_c2)


#Case 1
u0, v0 = inflow_uvw(u_mag, psi, 90)[0:2]
fx_a, fy_a = compute_forces_xy(porous_dimensions, rho, u0, v0,f_ij)

#Case 2
u0, v0, w0 = inflow_uvw(u_mag, psi_c2, theta)
fz_a = compute_force_z(porous_dimensions,rho,u0,v0,w0,f_ij)

# Read sweep_results
porous_cfd_1 = pd.read_csv("psi_sweep_results_benchmark.csv")

# Read sweep_results
porous_cfd_2 = pd.read_csv("theta_sweep_results_benchmark.csv")

fx_to_plot = {"2D Analytical": [psi, fx_a], "3D Porous CFD": [porous_cfd_1['Psi'], porous_cfd_1['Fx']]}
fy_to_plot = {"2D Analytical": [psi, fy_a], "3D Porous CFD": [porous_cfd_1['Psi'], porous_cfd_1['Fy']]}
fz_to_plot = {"2D Analytical": [theta, fz_a], "3D Porous CFD Fy": [porous_cfd_2['Theta'], porous_cfd_2['Fy']],
              "3D Porous CFD Fz": [porous_cfd_2['Theta'], porous_cfd_2['Fz']]}

f_all = {"Fx Analytical - Case 1": [psi, fx_a],
         "Fx Porous CFD - Case 1": [porous_cfd_1['Psi'], porous_cfd_1['Fx']],
         "Fy Analytical - Case 1": [psi, fy_a],
         "Fy Porous CFD - Case 1": [porous_cfd_1['Psi'], porous_cfd_1['Fy']],
         "Fz Analytical - Case 2": [theta - 90, fz_a],
         "Fz Porous CFD - Case 2": [porous_cfd_2['Theta']-90, porous_cfd_2['Fz']]}

plot_forces("Fx", fx_to_plot)
plot_forces("Fy", fy_to_plot)
plot_forces("Fz", fz_to_plot)
plot_forces("Analytical vs. CFD Forces", f_all)
