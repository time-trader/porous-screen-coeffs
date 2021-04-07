# This script finds Fij for the full porous matrix, if the forces on the porous screen are given.
# Copyright: Yuriy Marykovskiy, ETH Zurich, 2020


import numpy as np
import pandas as pd
from plot import plot_forces
from scipy.optimize import least_squares


def constants(u_mag, alpha, f_yy, f_yx):

    u0 = u_mag * np.cos(alpha)
    v0 = u_mag * np.sin(alpha)
    c11 = 0.5*(f_yy*u0 - f_yx*v0)
    c12 = 0.5*(f_yx*u0 + f_yy*v0)
    c1 = (-c11 + np.sqrt(c11 ** 2 + c12 ** 2)) / c12
    lamb = 0.5 * np.sqrt(f_yy ** 2 + f_yx ** 2)

    return lamb, c1


def velocity_u(u_mag, alpha):
    return u_mag * np.cos(alpha)


def velocity_v(x, u_mag, alpha, f_yy, f_yx):

    lamb, c1 = constants(u_mag, alpha, f_yy, f_yx)

    v_dxp = u_mag * np.cos(alpha)*\
            (2*f_yy*c1*np.exp(lamb*x)-f_yx*np.exp(2*lamb*x)+f_yx*c1**2)/\
            (2*f_yx*c1*np.exp(lamb*x)+ f_yy*np.exp(2*lamb*x)-f_yy*c1**2)
    return v_dxp


def pressure(x, rho, u_mag, alpha, f_ij):

    lamb, c1 = constants(u_mag, alpha, f_ij[1, 1], f_ij[1, 0])

    integral1 = lamb*x - \
              np.log(2*f_ij[1, 0]*c1*np.exp(lamb*x) + f_ij[1, 1]*np.exp(2*lamb*x) - f_ij[1, 1]*c1**2) + \
              np.log(2)
    integral2=0.5*f_ij[1,0]*integral1-\
              (4*c1*np.exp(lamb*x)*lamb**2)/\
              (2*f_ij[1, 0]*c1*np.exp(lamb*x) + f_ij[1, 1]*np.exp(2*lamb*x) - f_ij[1, 1]*c1**2)
    p_x = rho * velocity_u(u_mag, alpha)**2 *\
            (-f_ij[0, 0]/f_ij[1, 1]*integral1+2*f_ij[0, 1]/f_ij[1, 1]**2 * integral2)

    return p_x


def compute_field(u_mag, rho, dxp, alpha_range, f_ij):
    p_x0 = []
    u_dxp = []
    v_dxp = []
    alpha_range = np.deg2rad(alpha_range)

    for alpha in alpha_range:
        u_dxp.append(velocity_u(u_mag, alpha))
        v_dxp.append(velocity_v(dxp, u_mag, alpha, f_ij))
        p_x0.append(pressure(dxp, rho, u_mag, alpha, f_ij)-
                    pressure(0, rho, u_mag, alpha, f_ij))

    return p_x0, u_dxp, v_dxp


def compute_forces(porous_volume_dimensions, rho, u_mag, alpha, beta, a_beta, f_ij):
    alpha=np.deg2rad(alpha)
    beta=np.deg2rad(beta)
    a_beta=np.deg2rad(a_beta)
    fx=[]
    fy=[]
    fz=[]

    # Assume independent sweeps in different planes
    u0=u_mag * np.cos(alpha)
    v0=u_mag * np.sin(alpha)

    u00=u_mag * np.cos(a_beta)*np.cos(beta)
    w0=u_mag * np.sin(beta)

    v_dxp=velocity_v(porous_volume_dimensions[0], u_mag, alpha, f_ij[1,1], f_ij[1,0])
    w_dxp = velocity_v(porous_volume_dimensions[0], u_mag, beta, f_ij[2,2], 0)
    dv=v_dxp-v0
    dw=w_dxp-w0
    dp=(pressure(porous_volume_dimensions[0], rho, u_mag, alpha, f_ij) -
        pressure(0, rho, u_mag, alpha, f_ij))

    fx = dp * porous_volume_dimensions[1] * porous_volume_dimensions[2]
    fy = -rho * u0 * dv * porous_volume_dimensions[1] * porous_volume_dimensions[2]
    fz = -rho * u00 * dw * porous_volume_dimensions[1] * porous_volume_dimensions[2]

    return fx, fy, fz


def alpha_v_const(beta, f_yy, f_yx):
    alpha=2*np.arctan((f_yy-np.sqrt((f_yx*np.cos(beta))**2+f_yy**2))/(f_yx*np.cos(beta)))
    return alpha

def residuals_function(f_ij, alpha, beta, a_beta, f, u_mag, rho, porous_dims, weights):
    fx, fy, fz = compute_forces(porous_dims, rho, u_mag, alpha, beta, a_beta, np.array([[f_ij[0], f_ij[1], 0],
                                                                                [f_ij[2], f_ij[3], 0],
                                                                                [0, 0, f_ij[4]]]))
    fi = apply_weights(np.concatenate((fx, fy, fz), axis=None), weights)
    return fi-f


def apply_weights(array, weights, scale=1):
    array*=weights
    array*=scale
    return array


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Define training data from HiFi CFD or Wind Tunnel
    u_mag_train = 10
    rho_train = 1
    porous_dimensions = [0.2, 0.0842527, 0.25]
    alpha_train = np.array([-80, -60, -40, -20, 1E-6, 20, 40, 60, 80])
    beta_train = np.array([-80, -60, -40, -20, 1E-6, 20, 40, 60, 80])

    alpha_beta = 48

    fx_train = np.array([0.656813, 2.074982, 3.209272, 3.574992, 2.968339, 1.698684, 0.482158, -0.31434, -0.19178])
    fx_weights = np.ones(fx_train.size)
    fy_train = np.array([-0.67436, -2.10182, -3.16394, -3.41387, -2.72193, -1.42857, -0.27297, 0.408222, 0.250917])
    fy_weights = np.array([1, 1, 1, 1, 1, 1, 10, 10, 1])
    fz_train = np.array([-0.097, -0.146, -0.137, -0.076, 0, 0.076, 0.137, 0.146, 0.097])
    fz_weights = np.ones(fx_train.size)

    # Append and apply weights
    weights=np.concatenate((fx_weights, fy_weights, fz_weights), axis=None)
    f_train=np.concatenate((fx_train, fy_train, fz_train), axis=None)
    fi_train = apply_weights(f_train, weights)

    f_ij=np.ones(5)
    print("Finding best fit:")
    res_lsq = least_squares(residuals_function, f_ij,
                            args=(alpha_train, beta_train, alpha_beta, fi_train, u_mag_train, rho_train,
                                  porous_dimensions, weights),
                            verbose=1, ftol=1e-12)  # ,
                            # loss='soft_l1')
                            # bounds=([10,-1000,-1000,10],[1000,-10,-10,1000]))

    # Compute forces based on the best fit coef.
    alpha_a = np.linspace(-80, 80, 161)
    fx_a, fy_a, fz_a = compute_forces(porous_dimensions, rho_train, u_mag_train, alpha_train, beta_train, alpha_beta,
                                np.array([[res_lsq.x[0], res_lsq.x[1], 0],
                                          [res_lsq.x[2], res_lsq.x[3], 0],
                                          [0, 0, res_lsq.x[4]]]))

    #Case 2 condition
    print(np.rad2deg(alpha_v_const(beta_train, res_lsq.x[3], res_lsq.x[2])))


    # Read porous CFD data
    porous_cfd=pd.read_csv("sweep_results", sep="[\s+|\t|(|)]", header=None, engine="python")
    fy_porous=-(porous_cfd[6]/(porous_dimensions[1]*porous_dimensions[2]) -
                u_mag_train*np.sin(np.deg2rad(porous_cfd[0])))\
              *10*np.cos(np.deg2rad(porous_cfd[0]))*rho_train*(porous_dimensions[1]*porous_dimensions[2])

    # Do some Plotting
    fx_to_plot = {"HiFidelity CFD": [alpha_train, fx_train], "2D Analytical": [alpha_train, fx_a],
                  "3D Porous CFD": [porous_cfd[0], porous_cfd[2]]}
    fy_to_plot = {"HiFidelity CFD": [alpha_train, fy_train], "2D Analytical": [alpha_train, fy_a],
                  "3D Porous CFD": [porous_cfd[0], fy_porous]}
    fz_to_plot = {"HiFidelity CFD": [beta_train, fz_train], "2D Analytical": [beta_train, fz_a]}

    plot_forces("Fx", fx_to_plot)
    plot_forces("Fy", fy_to_plot)
    plot_forces("Fz", fz_to_plot)

    print(res_lsq.x)




