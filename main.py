# This script finds Fij for the full porous matrix, if the forces on the porous screen are given.
# Copyright: Yuriy Marykovskiy, ETH Zurich, 2020

import numpy as np
import pandas as pd
from porousZones.plot import plot_forces
from porousZones.analytical import compute_forces_xy, compute_force_z, alpha_v_const, inflow_uvw
from porousZones.residual_functions import residuals_function, residuals_function_fz, apply_weights
from scipy.optimize import least_squares

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Define training data from HiFi CFD or Wind Tunnel
    u_mag_train = 10
    rho_train = 1
    porous_dimensions = [0.2, 0.0842527, 0.25]

    #Case 1
    alpha_train = np.array([-80, -60, -40, -20, 1E-6, 20, 40, 60, 80])
    fx_train = np.array([0.656813, 2.074982, 3.209272, 3.574992, 2.968339, 1.698684, 0.482158, -0.31434, -0.19178])
    fx_weights = np.ones(fx_train.size)
    fy_train = np.array([-0.67436, -2.10182, -3.16394, -3.41387, -2.72193, -1.42857, -0.27297, 0.408222, 0.250917])
    fy_weights = np.array([1, 1, 1, 1, 1, 1, 10, 10, 1])


    # Append and apply weights
    train_weights=np.concatenate((fx_weights, fy_weights), axis=None)
    f_train=np.concatenate((fx_train, fy_train), axis=None)
    fi_train = apply_weights(f_train, train_weights)

    f_coefficients = np.ones(4)
    print("Finding best fit:")
    res_lsq = least_squares(residuals_function, f_coefficients,
                            args=(alpha_train, fi_train, u_mag_train, rho_train,
                                  porous_dimensions, train_weights),
                            verbose=1, ftol=1e-12)  # ,
                            # loss='soft_l1')
                            # bounds=([10,-1000,-1000,10],[1000,-10,-10,1000]))

    # Compute forces based on the best fit coef.
    alpha_a = np.linspace(-80, 80, 161)
    u0, v0 = inflow_uvw(u_mag_train, alpha_train, 90)[0:2]
    fx_a, fy_a = compute_forces_xy(porous_dimensions, rho_train, u0, v0,
                                    np.array([[res_lsq.x[0], res_lsq.x[1], 0],
                                              [res_lsq.x[2], res_lsq.x[3], 0],
                                              [0, 0, 0]]))

    print(res_lsq.x)



    #Case 2 condition
    beta_train = np.array([10, 30, 50, 70, 89.999999, 110, 130, 150, 170])
    alpha_beta_train = np.rad2deg(alpha_v_const(res_lsq.x[3], res_lsq.x[2]))*np.ones(len(beta_train))
    print(alpha_beta_train[0])

    fz_train = np.array([0.097, 0.146, 0.137, 0.076, 0, -0.076, -0.137, -0.146, -0.097])
    fz_weights = np.ones(fz_train.size)
    f_zz=1

    u0, v0, w0 = inflow_uvw(u_mag_train, alpha_beta_train, beta_train)
    res_lsq_z = least_squares(residuals_function_fz, f_zz,
                              args=(res_lsq.x, alpha_beta_train,  beta_train, fz_train, u_mag_train, rho_train,
                                    porous_dimensions, fz_weights),
                              verbose=1, ftol=1e-12)
    print(res_lsq_z.x)
    fz_a = compute_force_z(porous_dimensions, rho_train, u0, v0, w0,
                                    np.array([[res_lsq.x[0], res_lsq.x[1], 0],
                                              [res_lsq.x[2], res_lsq.x[3], 0],
                                              [0, 0, res_lsq_z.x[0]]]))




    # Read porous CFD data
    porous_cfd=pd.read_csv("psi_sweep_results_extended_mesh", sep="[\s|\t|(|)]", header=None, engine="python")
    porous_cfd_2 = pd.read_csv("theta_sweep_results_extended_mesh", sep="[\s|\t|(|)]", header=None, engine="python")

    fy_porous = -(porous_cfd[6]/(porous_dimensions[1]*porous_dimensions[2]) -
                  u_mag_train*np.sin(np.deg2rad(porous_cfd[0])))\
                  *10*np.cos(np.deg2rad(porous_cfd[0]))*rho_train*(porous_dimensions[1]*porous_dimensions[2])

    fz_porous = -(porous_cfd_2[7] / (porous_dimensions[1] * porous_dimensions[2]) -
                  u_mag_train * np.cos(np.deg2rad(porous_cfd_2[0]))) \
                  *10*np.cos(np.deg2rad(alpha_beta_train[0])) * np.sin(np.deg2rad(porous_cfd_2[0])) \
                  * rho_train * (porous_dimensions[1] * porous_dimensions[2])

    # Do some Plotting
    fx_to_plot = {"HiFidelity CFD": [alpha_train, fx_train], "2D Analytical": [alpha_train, fx_a],
                  "3D Porous CFD": [porous_cfd[0], porous_cfd[2]]}
    fy_to_plot = {"HiFidelity CFD": [alpha_train, fy_train], "2D Analytical": [alpha_train, fy_a],
                  "3D Porous CFD": [porous_cfd[0], fy_porous]}
    fz_to_plot = {"HiFidelity CFD": [beta_train, fz_train], "2D Analytical": [beta_train, fz_a],
                  "3D Porous CFD": [porous_cfd_2[0], fz_porous]}

    plot_forces("Fx", fx_to_plot)
    plot_forces("Fy", fy_to_plot)
    plot_forces("Fz", fz_to_plot)




