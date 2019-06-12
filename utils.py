"""Airfoil calculation functions."""

# pylint: disable=R0902
# - allow more than 7 instance variables
# pylint: disable=R0913
# - allow more than 5 arguments to a function
# pylint: disable=R0914
# - allow more than 15 local variables
# pylint: disable=R0915
# - allow more than 50 statemnts in class or function
# pylint: disable=C0103
# - allow variable names less than 3 characters
# pylint: disable=C0412
# - allow imports to be sorted by isort and not grouped

import math as m
from math import radians as r

import numpy as np


def cot(rad):
    """Return cotangent of an angle in radians."""
    return m.cos(rad)/m.sin(rad)


def sec(rad):
    """Return Secant of an angle in radians."""
    return 1/m.cos(rad)


def get_vars(stage, blade: str, station: str):
    """Obtain flow variables from stage class for blade locations."""
    if station == 'root':
        loc = stage.root

    elif station == 'mean':
        loc = stage.mean

    elif station == 'tip':
        loc = stage.tip

    radius = loc.radius
    if blade == 'IGV':
        avle, avte, rvle, rvte, v1, v2 = (0,
                                          loc.alpha1,
                                          loc.beta2,
                                          loc.beta1,
                                          loc.cx,
                                          loc.c1)

    elif blade == 'Rotor':
        avle, avte, rvle, rvte, v1, v2 = (loc.beta1,
                                          loc.beta2,
                                          loc.alpha1,
                                          loc.alpha2,
                                          loc.w1,
                                          loc.w2)

    elif blade == 'OSV' or 'Stator':
        avle, avte, rvle, rvte, v1, v2 = (loc.alpha2,
                                          loc.alpha1,
                                          loc.beta2,
                                          loc.beta1,
                                          loc.c2,
                                          loc.c1)

    return avle, avte, rvle, rvte, v1, v2, radius


def calcchord(radius, hub, z, aoa):
    """Calculate max chord.

    Chord is to be definied such that either:
        Chord touches edge of blade row, or
        Chord touches edge of adjacent blade
    """
    theta = m.atan(hub/(2*radius*m.tan(m.pi/z)))
    if aoa <= theta:
        chord = (2*radius*m.tan(m.pi/z))/(m.cos(aoa))
    elif aoa > theta:
        chord = hub/m.sin(aoa)

    return chord


def calcincidence(thickness, solidity, relative_inlet_angle):
    """Calculate incidence of flow over airfoil."""
    k_sh = 1.0
    k_ti = (-0.0214
            + 19.17*(thickness)
            - 122.3*(thickness)**2
            + 312.5*(thickness)**3)
    i_010 = ((0.0325 - 0.0674*solidity)
             + (-0.002364 + 0.0913*solidity)*relative_inlet_angle
             + (.0000164 - 0.000238*solidity)*relative_inlet_angle**2)
    n_slope = ((-0.063 - 0.02274*solidity)
               + (-0.0035 + 0.0029*solidity)*relative_inlet_angle
               - (0.0000379 + 0.0000111*solidity)**relative_inlet_angle**2)
    incidence = k_sh*k_ti*i_010  # + n_slope*camber
    return incidence, n_slope


def calcdeviation(thickness, solidity, relative_inlet_angle):
    """Calculate deviation of flow over airfoil."""
    k_sh = 1.0
    k_td = (0.0142
            + 6.172*(thickness)
            + 36.61*(thickness)**2)
    d_010 = ((-0.0443 + 0.1057*solidity)
             + (0.0209 - 0.0186*solidity)*relative_inlet_angle
             + (-0.0004 + 0.0007*solidity)*relative_inlet_angle**2)
    m_prime = (0.17 - (3.33*10**(-4))*(1.0 - 0.1
                                       * relative_inlet_angle)
               * relative_inlet_angle)
    b = (0.9655
         + (2.538*10**-3)*relative_inlet_angle
         + (4.221*10**-5)*relative_inlet_angle**2
         - (1.3*10**-6)*relative_inlet_angle**3)
    m_slope = m_prime/(solidity**b)
    deviation = k_sh*k_td*d_010  # + m_slope*camber
    return deviation, m_slope
    # Carter's Rule - only valid at optimal incidence
    # m_dev = ((0.23*(2*self.xc))**2 + rvte/500)
    # if blade != 'IGV':
    #     self.deviation = m_dev*self.camber*(1/self.sigma)**0.5
    # else:
    #     self.deviation = m_dev*self.camber*(1/self.sigma)


def cambercurvature(xc, yc):
    """Calculate camber curvature based on modified NACA 4."""
    kc = -yc*(1/(xc**2) + 1/((1 - xc)**2))

    return kc


def thicknesscurvature(xt, yt):
    """Calculate thickness curvature based on modified NACA 4."""
    d1 = (2.24 - 5.42*xt + 12.3*xt**2)/(10*(1 - 0.878*xt))
    d2 = (0.294 - 2*(1 - xt)*d1)/((1 - xt)**2)
    d3 = (-0.196 + (1 - xt)*d1)/((1 - xt)**3)

    p1 = (1/5)*((1 - xt)**2)/(0.588 - 2*(1 - xt)*d1)

    a0 = 0.296904
    # a1 = 0.3/xt - (15/8)*(a0/m.sqrt(xt)) - xt/(10*p1)
    a2 = -0.3/(xt**2) + (5/4)*(a0/(xt**(3/2))) + 1/(5*p1)
    a3 = 0.1/(xt**3) - (3/8)*(a0/(xt**(5/2))) - 1/(10*p1*xt)

    kt = ((5*(2*yt)*(-a0/(4*xt**(3/2)) + 2*a2 + 6*a3*xt))
          + (5*(2*yt)*(2*d2 + 6*d3*(1 - xt))))/2

    return kt


def curvatures(xc, yc, xt, yt):
    """Calculate camber and thickness curvatures based on modified NACA 4."""
    d1 = (2.24 - 5.42*xt + 12.3*xt**2)/(10*(1 - 0.878*xt))
    d2 = (0.294 - 2*(1 - xt)*d1)/((1 - xt)**2)
    d3 = (-0.196 + (1 - xt)*d1)/((1 - xt)**3)

    p1 = (1/5)*((1 - xt)**2)/(0.588 - 2*(1 - xt)*d1)

    a0 = 0.296904
    # a1 = 0.3/xt - (15/8)*(a0/m.sqrt(xt)) - xt/(10*p1)
    a2 = -0.3/(xt**2) + (5/4)*(a0/(xt**(3/2))) + 1/(5*p1)
    a3 = 0.1/(xt**3) - (3/8)*(a0/(xt**(5/2))) - 1/(10*p1*xt)

    kc = -yc*(1/(xc**2) + 1/((1 - xc)**2))
    kt = ((5*(2*yt)*(-a0/(4*xt**(3/2)) + 2*a2 + 6*a3*xt))
          + (5*(2*yt)*(2*d2 + 6*d3*(1 - xt))))/2

    return kc, kt


def camberbezier(xc, yc, kc, cle, cte):
    """Calculate camber bezier variable in BP3333 airfoil."""
    bc = 0
    co = cot(r(cle)) + cot(r(cte))
    bcu = (16 + 3*kc*co + 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)
    bcl = (16 + 3*kc*co - 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)

    cbounds = (bcl, bcu)
    bci_list = []

    for bci in cbounds:
        bci_list.append(round(bci, 4))
        if 0 < bci < yc:
            bc = float(bci)
            xlc2 = xc - m.sqrt((2/3)*(bci - yc)/(kc))

            if xlc2 <= 0:
                print(f'xc = {xc}')
                print(f'yc = {yc}')
                print(f'kc = {kc:.4f}')
                print(f'bc = {bci:.4f}')
                print(f'x - sqrt((2/3)(b - y)/k) = {xlc2}')
                raise ValueError(' Invalid Airfoil Shape'
                                 f'x - sqrt((2/3)(b - y)/k) = {xlc2}')

    if bc == 0:
        print(f'xc = {xc}'.format(xc=xc))
        print(f'yc = {yc}'.format(yc=yc))
        print(f'kc = {kc}'.format(kc=round(kc, 4)))
        print(f'bc = {bci}'.format(bci=bci_list))
        raise ValueError(f'No bc found within bounds 0 < bc < {yc: .4f}.')


def thicknessbezier(xt, yt, kt, rle, blade: str, station: str):
    """Calculate thickness bezier variable in BP3333 airfoil."""
    bt = 0
    b_poly = ((27 / 4) * kt**2,
              -27 * kt**2 * xt,
              9 * kt * yt + (81 / 2) * kt**2 * xt**2,
              2 * -rle - 18 * kt * xt * yt - 27 * kt**2 * xt**3,
              3 * yt**2 + 9 * kt * xt**2 * yt + (27 / 4) * kt**2 * xt**4)

    real_troots = [root.real for root in np.roots(b_poly) if root.imag == 0]
    low_b = max(0, xt - m.sqrt((-2/3)*(yt/kt)))
    bti_list = []

    for bti in real_troots:
        bti_list.append(round(bti, 4))
        if not real_troots:
            print(f'xt = {xt}')
            print(f'yt = {yt}')
            print(f'kt = {kt:.4f}')
            raise ValueError('Invalid Airfoil Shape.'
                             f'{blade} {station} bt is complex.')

        if low_b < bti < xt:
            bt = float(bti)
            ylt1 = yt + (3/2)*kt*(xt - bti)**2

            if ylt1 <= 0:
                print(f'xt = {xt}')
                print(f'yt = {yt}')
                print(f'kt = {kt:.4f}')
                print(f'bt = {bti:.4f}')
                print(f'(y + (3/2)k(x - b)^2) = {ylt1}')
                raise ValueError(' Invalid Airfoil Shape'
                                 f'(y + (3/2)k(x - b)^2) = {ylt1}.')

    if bt == 0:
        print(f'xt = {xt}')
        print(f'yt = {yt}')
        print(f'kt = {kt:.4f}')
        print(f'bt = {bti_list}')
        raise ValueError(f'No bt found within bounds {low_b:.4f} < bt < {xt}.')

    return bt


def beziers(xc, yc, kc, xt, yt, kt, cle, cte, rle, blade: str, station: str):
    """Calculate camber and thickness Bezier variables in BP3333 airfoil."""
    bc = 0
    bt = 0

    co = cot(r(cle)) + cot(r(cte))
    bcu = (16 + 3*kc*co + 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)
    bcl = (16 + 3*kc*co - 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)

    cbounds = (bcl, bcu)
    bci_list = []

    for bci in cbounds:
        bci_list.append(round(bci, 4))
        if 0 < bci < yc:
            bc = float(bci)
            xlc2 = xc - m.sqrt((2/3)*(bci - yc)/(kc))

            if xlc2 <= 0:
                print(f'xc = {xc}')
                print(f'yc = {yc}')
                print(f'kc = {kc:.4f}')
                print(f'bc = {bci:.4f}')
                print(f'x - sqrt((2/3)(b - y)/k) = {xlc2}')
                raise ValueError('Invalid Airfoil Shape'
                                 f'x - sqrt((2/3)(b - y)/k) = {xlc2}')

    if bc == 0:
        print(f'xc = {xc}')
        print(f'yc = {yc}')
        print(f'kc = {kc:.4f}')
        print(f'bc = {bci_list}')
        raise ValueError(f'No bc found within bounds 0 < bc < {yc: .4f}.')

    b_poly = ((27 / 4) * kt**2,
              -27 * kt**2 * xt,
              9 * kt * yt + (81 / 2) * kt**2 * xt**2,
              2 * -rle - 18 * kt * xt * yt - 27 * kt**2 * xt**3,
              3 * yt**2 + 9 * kt * xt**2 * yt + (27 / 4) * kt**2 * xt**4)

    real_troots = [root.real for root in np.roots(b_poly) if root.imag == 0]
    low_b = max(0, xt - m.sqrt((-2/3)*(yt/kt)))
    bti_list = []

    for bti in real_troots:
        bti_list.append(round(bti, 4))
        if not real_troots:
            print(f'xt = {xt}')
            print(f'yt = {yt}')
            print(f'kt = {kt:.4f}')
            raise ValueError('Invalid Airfoil Shape.'
                             f'{blade} {station} bt is complex.')

        if low_b < bti < xt:
            bt = float(bti)
            ylt1 = yt + (3/2)*kt*(xt - bti)**2

            if ylt1 <= 0:
                print(f'xt = {xt}')
                print(f'yt = {yt}')
                print(f'kt = {kt:.4f}')
                print(f'bt = {bti:.4f}')
                print(f'(y + (3/2)k(x - b)^2) = {ylt1}')
                raise ValueError('Invalid Airfoil Shape'
                                 f'(y + (3/2)k(x - b)^2) = {ylt1:.4f}')

    if bt == 0:
        print(f'xt = {xt}')
        print(f'yt = {yt}')
        print(f'kt = {kt:.4f}')
        print(f'bt = {bti_list}')
        raise ValueError(f'No bt found within bounds {low_b:.4f} < bt < {xt}.')

    return bc, bt
