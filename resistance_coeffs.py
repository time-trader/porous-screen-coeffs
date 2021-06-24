# This script finds Fij for the full porous matrix, if the forces on the porous screen are given.
# Copyright: Yuriy Marykovskiy, ETH Zurich, 2020

import numpy as np
from porousZones.analytical import alpha_v_const, inflow_uvw
from porousZones.residual_functions import residuals_function, residuals_function_fz, apply_weights
from scipy.optimize import least_squares

# Define training data from HiFi CFD or Wind Tunnel
u_mag_train = 10
rho_train = 1
porous_dimensions = [0.2, 0.0842527, 0.25]

#Case 1
psi_train_c1 = np.array([-80, -60, -40, -20, 1E-6, 20, 40, 60, 80])
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
                        args=(fi_train, psi_train_c1, u_mag_train, rho_train,
                              porous_dimensions, train_weights),
                        verbose=1, ftol=1e-12)  # ,
                        # loss='soft_l1')
                        # bounds=([10,-1000,-1000,10],[1000,-10,-10,1000]))

print(res_lsq.x)

# Case 2 condition
theta_train_c2 = np.array([10, 30, 50, 70, 89.999999, 110, 130, 150, 170])
psi_train_c2 = np.rad2deg(alpha_v_const(res_lsq.x[3], res_lsq.x[2])) * np.ones(len(theta_train_c2))
print(psi_train_c2[0])

fz_train = np.array([0.097, 0.146, 0.137, 0.076, 0, -0.076, -0.137, -0.146, -0.097])
fz_weights = np.ones(fz_train.size)
f_zz = 1

u0, v0, w0 = inflow_uvw(u_mag_train, psi_train_c2, theta_train_c2)
res_lsq_z = least_squares(residuals_function_fz, f_zz,
                          args=(res_lsq.x, psi_train_c2, theta_train_c2, fz_train, u_mag_train, rho_train,
                                porous_dimensions, fz_weights),
                          verbose=1, ftol=1e-12)
print(res_lsq_z.x)





