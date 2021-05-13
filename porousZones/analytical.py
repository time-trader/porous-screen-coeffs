# Copyright: Yuriy Marykovskiy, ETH Zurich, 2020

import numpy as np


def inflow_uvw(u_mag, psi, theta):
    """
    Transforms U from Spherical to Cartesian coordinates
    :param u_mag:
    :param psi:
    :param theta:
    :return: u0, v0, w0 - Inflow vector in cartesian coordinates (u_x(0), u_y(0), u_z(0))
    """

    psi = np.deg2rad(psi)
    theta = np.deg2rad(theta)

    u0 = u_mag * np.cos(psi) * np.sin(theta)
    v0 = u_mag * np.sin(psi) * np.sin(theta)
    w0 = u_mag * np.cos(theta)

    return u0, v0, w0


def constants_1(u0, v0, f_ij):

    c11 = 0.5*(f_ij[1,1]*u0 - f_ij[1,0]*v0)
    c12 = 0.5*(f_ij[1,0]*u0 + f_ij[1,1]*v0)
    c1 = (-c11 + np.sqrt(c11 ** 2 + c12 ** 2)) / c12
    lamb = 0.5 * np.sqrt(f_ij[1,1] ** 2 + f_ij[1,0] ** 2)

    return lamb, c1


def constants_2(u0, v0, w0, f_ij):

    c11 =np.sqrt(u0**2+v0**2)
    c12 =np.sqrt(u0**2+v0**2+w0**2)
    c1 = - (c11 - c12) / w0
    lamb = 0.5 *(1/u0) * np.sqrt(u0 ** 2 + v0 ** 2)*f_ij[2,2]

    return lamb, c1


def velocity_v(x, u0, v0, f_ij):

    lamb, c1 = constants_1(u0, v0,  f_ij)

    v_dxp = u0*\
            (2*f_ij[1, 1]*c1*np.exp(lamb*x)-f_ij[1, 0]*np.exp(2*lamb*x)+f_ij[1, 0]*c1**2)/\
            (2*f_ij[1, 0]*c1*np.exp(lamb*x)+f_ij[1, 1]*np.exp(2*lamb*x)-f_ij[1, 1]*c1**2)
    return v_dxp


def velocity_w(x, u0, v0, w0, f_ij):

    lamb, c1 = constants_2(u0, v0, w0, f_ij)

    w_dxp = np.sqrt(u0**2+v0**2)*(2*c1*np.exp(lamb*x))/(np.exp(2*lamb*x)-c1**2)
    return w_dxp


def pressure(x, rho, u0, v0, f_ij):

    lamb, c1 = constants_1(u0, v0, f_ij)

    integral1 = lamb*x - \
              np.log(2*f_ij[1, 0]*c1*np.exp(lamb*x) + f_ij[1, 1]*np.exp(2*lamb*x) - f_ij[1, 1]*c1**2) + \
              np.log(2)
    integral2=0.5*f_ij[1,0]*integral1-\
              (4*c1*np.exp(lamb*x)*lamb**2)/\
              (2*f_ij[1, 0]*c1*np.exp(lamb*x) + f_ij[1, 1]*np.exp(2*lamb*x) - f_ij[1, 1]*c1**2)
    p_x = rho * u0**2 *\
            (-f_ij[0, 0]/f_ij[1, 1]*integral1+2*f_ij[0, 1]/f_ij[1, 1]**2 * integral2)

    return p_x


def compute_forces_xy(porous_volume_dimensions, rho, u0, v0, f_ij):
    """
    Compute forces in Case 1 flow: w=0!

    :param porous_volume_dimensions:
    :param rho:
    :param u0:
    :param v0:
    :param f_ij:
    :return:
    """
    v_dxp = velocity_v(porous_volume_dimensions[0], u0, v0, f_ij)

    dv = v_dxp-v0

    dp = (pressure(porous_volume_dimensions[0], rho, u0, v0, f_ij) -
          pressure(0, rho, u0, v0, f_ij))

    fx = dp * porous_volume_dimensions[1] * porous_volume_dimensions[2]
    fy = -rho * u0 * dv * porous_volume_dimensions[1] * porous_volume_dimensions[2]
    return fx, fy


def compute_force_z(porous_volume_dimensions, rho, u0, v0, w0, f_ij):

    w_dxp = velocity_w(porous_volume_dimensions[0], u0, v0, w0, f_ij)
    dw = w_dxp - w0
    fz = -rho * u0 * dw * porous_volume_dimensions[1] * porous_volume_dimensions[2]
    return fz


def alpha_v_const(f_yy, f_yx):
    """
    Finds angle psi, to verify Case 2 condition
    :param f_yy:
    :param f_yx:
    :return:
    """
    psi = 2*np.arctan((f_yy-np.sqrt(f_yx**2+f_yy**2))/f_yx)
    return psi

