import numpy as np
import pandas as pd
from porousZones.analytical import compute_forces_xy, compute_force_z, inflow_uvw, alpha_v_const
from porousZones.plot import plot_forces

u_mag = 10
rho = 1
porous_dimensions=[0.5, 1, 1]

f_ij=np.array([[20, -18, 0],[-18, 16, 0],[0,0,2]])
alpha=np.linspace(-89, 89, 179)
beta=np.linspace(1,179,179)
alpha_beta=np.rad2deg(alpha_v_const(f_ij[1,1],f_ij[1,0]))
print(alpha_beta)


#Case 1
u0, v0 = inflow_uvw(u_mag, alpha, 90)[0:2]
fx_a, fy_a = compute_forces_xy(porous_dimensions, rho, u0, v0,f_ij)

#Case 2
u0, v0, w0 = inflow_uvw(u_mag, alpha_beta, beta)
fz_a = compute_force_z(porous_dimensions,rho,u0,v0,w0,f_ij)

# Read sweep_results
porous_cfd_1 = pd.read_csv("psi_sweep_results_benchmark", sep="[\s|\t|(|)]", header=None, engine="python")
fy_porous = -(porous_cfd_1[6] / (porous_dimensions[1] * porous_dimensions[2]) -
              u_mag * np.sin(np.deg2rad(porous_cfd_1[0]))) \
            * 10 * np.cos(np.deg2rad(porous_cfd_1[0])) * (porous_dimensions[1] * porous_dimensions[2])


# Read sweep_results
porous_cfd_2 = pd.read_csv("theta_sweep_results_benchmark", sep="[\s|\t|(|)]", header=None, engine="python")

beta_porous = porous_cfd_2[0]
fy_porous_2 = -(porous_cfd_2[6] / (porous_dimensions[1] * porous_dimensions[2]) -
                u_mag * np.sin(np.deg2rad(alpha_beta)) * np.sin(np.deg2rad(beta_porous))) \
              * 10 * np.cos(np.deg2rad(alpha_beta))*np.sin(np.deg2rad(beta_porous))\
              * rho * (porous_dimensions[1] * porous_dimensions[2])
fz_porous = -(porous_cfd_2[7] / (porous_dimensions[1] * porous_dimensions[2]) -
              u_mag * np.cos(np.deg2rad(beta_porous))) \
              * 10 * np.cos(np.deg2rad(alpha_beta))*np.sin(np.deg2rad(beta_porous))\
              * rho * (porous_dimensions[1] * porous_dimensions[2])


fx_to_plot = {"2D Analytical": [alpha, fx_a], "3D Porous CFD": [porous_cfd_1[0], porous_cfd_1[2]]}
fy_to_plot = {"2D Analytical": [alpha, fy_a], "3D Porous CFD": [porous_cfd_1[0], fy_porous]}
fz_to_plot = {"2D Analytical": [beta, fz_a], "3D Porous CFD Fy": [beta_porous, fy_porous_2],
              "3D Porous CFD": [beta_porous, fz_porous]}
f_all={"Fx Analytical - Case 1": [alpha, fx_a], "Fx Porous CFD - Case 1": [porous_cfd_1[0], porous_cfd_1[2]],
       "Fy Analytical - Case 1": [alpha, fy_a], "Fy Porous CFD - Case 1": [porous_cfd_1[0], fy_porous],
       "Fz Analytical - Case 2": [beta-90, fz_a],  "Fz Porous CFD - Case 2": [beta_porous-90, fz_porous]}

plot_forces("Fx", fx_to_plot)
plot_forces("Fy", fy_to_plot)
plot_forces("Fz", fz_to_plot)
plot_forces("Analytical vs. CFD Forces", f_all)