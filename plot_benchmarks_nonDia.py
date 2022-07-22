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


# Case 1
u0, v0 = inflow_uvw(u_mag, psi, 90)[0:2]
fx_a, fy_a = compute_forces_xy(porous_dimensions, rho, u0, v0, f_ij)

# Case 2
u0, v0, w0 = inflow_uvw(u_mag, psi_c2, theta)
fz_a = compute_force_z(porous_dimensions, rho, u0, v0, w0, f_ij)

# Read sweep_results
porous_cfd_1 = pd.read_csv("data/psi_sweep_results_benchmark.csv", comment='#')
porous_cfd_2 = pd.read_csv("data/theta_sweep_results_benchmark.csv", comment='#')
porous_cfd_3 = pd.read_csv("data/psi_sweep_results_benchmark_fluent.csv", comment='#')
porous_cfd_4 = pd.read_csv("data/theta_sweep_results_benchmark_fluent.csv", comment='#')

porous_cfd_1.index = porous_cfd_1['Psi']
porous_cfd_2.index = porous_cfd_2['Theta']

porous_cfd_1 = porous_cfd_1.loc[[-89, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70]]
porous_cfd_2 = porous_cfd_2[::5]

fx_to_plot = {"2D Analytical": [psi, fx_a],
              "3D Porous CFD OpenFOAM": [porous_cfd_1['Psi'], porous_cfd_1['Fx']],
              "3D Porous CFD FLUENT": [porous_cfd_3['Psi'], porous_cfd_3['Fx']]}
fy_to_plot = {"2D Analytical": [psi, fy_a],
              "3D Porous CFD OpenFOAM": [porous_cfd_1['Psi'], porous_cfd_1['Fy']],
              "3D Porous CFD FLUENT": [porous_cfd_3['Psi'], porous_cfd_3['Fy']]}
fz_to_plot = {"2D Analytical": [theta, fz_a],
              "3D Porous CFD Fy": [porous_cfd_2['Theta'], porous_cfd_2['Fy']],
              "3D Porous CFD Fz OpenFOAM": [porous_cfd_2['Theta'], porous_cfd_2['Fz']],
              "3D Porous CFD Fz FLUENT": [porous_cfd_4['Theta'], porous_cfd_4['Fz']]}

f_xy = {"Fx Analytical - Case 1": [psi, fx_a],
         "Fx Porous CFD - Case 1": [porous_cfd_1['Psi'], porous_cfd_1['Fx']],
         "Fy Analytical - Case 1": [psi, fy_a],
         "Fy Porous CFD - Case 1": [porous_cfd_1['Psi'], porous_cfd_1['Fy']],
        }

analytical_forces = pd.DataFrame({"Psi":psi, "Fx":fx_a, "Fy":fy_a})
analytical_forces.to_csv("./data/psi_sweep_results_benchmark_non_dia_analytical.csv")
analytical_forces = pd.DataFrame({"Theta":theta, "Fz":fz_a})
analytical_forces.to_csv("./data/theta_sweep_results_benchmark_non_dia_analytical.csv")

plot_forces("Analytical vs. CFD Forces", f_xy, save_eps=False)
plot_forces("Analytical vs. CFD Forces", fz_to_plot, save_eps=False)
