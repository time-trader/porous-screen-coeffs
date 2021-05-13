import numpy as np
from porousZones.analytical import inflow_uvw, compute_forces_xy, compute_force_z


def residuals_function(f_opt, alpha, f, u_mag, rho, porous_dims, weights):

    u0, v0 = inflow_uvw(u_mag, alpha, 90)[0:2]
    fx, fy = compute_forces_xy(porous_dims, rho, u0, v0, np.array([[f_opt[0], f_opt[1], 0],
                                                                         [f_opt[2], f_opt[3], 0],
                                                                         [0, 0, 0]]))
    fi = apply_weights(np.concatenate((fx, fy), axis=None), weights)
    return fi-f


def residuals_function_fz(f_zz, f_opt, alpha, beta, f, u_mag, rho, porous_dims, weights):
    u0, v0, w0 = inflow_uvw(u_mag, alpha, beta)
    fz = compute_force_z(porous_dims, rho, u0, v0, w0, np.array([[f_opt[0], f_opt[1], 0],
                                                                         [f_opt[2], f_opt[3], 0],
                                                                         [0, 0, f_zz[0]]]))
    fi = apply_weights(fz, weights)
    return fi-f


def apply_weights(array, weights, scale=1):
    array *= weights
    array *= scale
    return array