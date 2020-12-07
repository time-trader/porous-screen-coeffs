# This script finds Fij for the full porous matrix, if the forces on the porous screen are given.
# Copyright: Yuriy Marykovskiy, ETH Zurich, 2020


import numpy as np
from plot import plot_forces
from scipy.optimize import least_squares

def constants(u_mag, alpha, f_ij):

    u0 = u_mag * np.cos(alpha)
    v0 = u_mag * np.sin(alpha)
    c11 = 0.5*(f_ij[1, 1]*u0 - f_ij[1, 0]*v0)
    c12 = 0.5*(f_ij[1, 0]*u0 + f_ij[1, 1]*v0)
    c1 = (-c11 + np.sqrt(c11 ** 2 + c12 ** 2)) / c12
    lamb = 0.5 * np.sqrt(f_ij[1, 1] ** 2 + f_ij[1, 0] ** 2)

    return lamb, c1


def velocity_u(u_mag, alpha):
    return u_mag * np.cos(alpha)


def velocity_v(x, u_mag, alpha, f_ij):

    lamb, c1 = constants(u_mag, alpha, f_ij)

    v_dxp = u_mag * np.cos(alpha)*\
            (2*f_ij[1, 1]*c1*np.exp(lamb*x)-f_ij[1, 0]*np.exp(2*lamb*x)+f_ij[1, 0]*c1**2)/\
            (2*f_ij[1, 0]*c1*np.exp(lamb*x)+f_ij[1, 1]*np.exp(2*lamb*x)-f_ij[1, 1]*c1**2)
    return v_dxp


def pressure(x, rho, u_mag, alpha, f_ij ):

    lamb, c1 = constants(u_mag, alpha, f_ij)

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

def compute_forces(porous_volume_dimensions, rho, u_mag, alpha, f_ij):
    alpha=np.deg2rad(alpha)
    fx=[]
    fy=[]

    u0=u_mag * np.cos(alpha)
    v0=u_mag * np.sin(alpha)
    v_dxp=velocity_v(porous_volume_dimensions[0], u_mag, alpha, f_ij)
    dv=v_dxp-v0
    dp=(pressure(porous_volume_dimensions[0], rho, u_mag, alpha, f_ij) -
        pressure(0, rho, u_mag, alpha, f_ij))

    fx = dp * porous_volume_dimensions[1] * porous_volume_dimensions[2]
    fy = -rho * u0 * dv * porous_volume_dimensions[1] * porous_volume_dimensions[2]
    return fx, fy


def residuals_function(f_ij, alpha, f):
    fx, fy = compute_forces([0.00778, 0.06, 0.00027], 1, 10, alpha, np.array([[f_ij[0], f_ij[1]], [f_ij[2], f_ij[3]]]))
    fi = apply_weights(np.append(fx, fy))
    return fi-f


def apply_weights(array):
    scale=1E4
    weights = np.array([1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1])
    array*=weights
    array*=scale
    return array


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    alpha_train = np.array([-45, -30, -15, 1E-6, 15, 30, 45])
    fx_train = np.array([2.201, 2.437, 2.474, 2.252, 1.849, 1.086, 2.620E-01]) * 1E-4
    fy_train = np.array([-2.265, -2.447, -2.325, -1.941, -1.507, -7.660E-01, 6.900E-03]) * 1E-4

    fx_porous_cfd = np.array([2.201, 2.437, 2.474, 2.252, 1.849, 1.086, 2.620E-01]) * 1E-4
    fy_porous_cfd = np.array([-2.265, -2.447, -2.325, -1.941, -1.507, -7.660E-01, 6.900E-03]) * 1E-4

    fi_train=apply_weights(np.append(fx_train, fy_train))
    f_ij=np.ones(4)
    print("Finding best fit:")
    res_lsq = least_squares(residuals_function, f_ij, args=(alpha_train, fi_train),
                            verbose=1)  # ,
                            # loss='soft_l1',
                            # bounds=([10,-1000,-1000,10],[1000,-10,-10,1000]))
    #all_fx=[fx_train, res_lsq.y[0:np.size(alpha_train)], fx_porous_cfd]
    #plot_forces(alpha_train, all_fx)

    print(res_lsq.x)




