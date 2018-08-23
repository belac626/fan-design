"""Airfoil calculation functions."""
import math as m
from math import radians as r

import numpy as np


def Cot(rad):
    """Return Cotangent of an angle in radians."""
    return m.cos(rad)/m.sin(rad)


def Sec(rad):
    """Return Secant of an angle in radians."""
    return 1/m.cos(rad)


def GetRootFlowVars(stage, blade: str, station: str):
    """Obtain flow variables from stage class for root blade locations."""
    rr = stage.root
    radius = rr.radius
    if blade == 'IGV':
        avle, avte, rvle, rvte, v1, v2 = (0,
                                          rr.alpha1,
                                          rr.beta2,
                                          rr.beta1,
                                          rr.cx,
                                          rr.c1)
    elif blade == 'Rotor':
        avle, avte, rvle, rvte, v1, v2 = (rr.beta1,
                                          rr.beta2,
                                          rr.alpha1,
                                          rr.alpha2,
                                          rr.w1,
                                          rr.w2)
    elif blade == 'OSV' or 'Stator':
        avle, avte, rvle, rvte, v1, v2 = (rr.alpha2,
                                          rr.alpha1,
                                          rr.beta2,
                                          rr.beta1,
                                          rr.c2,
                                          rr.c1)
    return avle, avte, rvle, rvte, v1, v2, radius


def GetMeanFlowVars(stage, blade: str, station: str):
    """Obtain flow variables from stage class for mean blade locations."""
    mm = stage.mean
    radius = mm.radius
    if blade == 'IGV':
        avle, avte, rvle, rvte, v1, v2 = (0,
                                          mm.alpha1,
                                          mm.beta2,
                                          mm.beta1,
                                          mm.cx,
                                          mm.c1)
    elif blade == 'Rotor':
        avle, avte, rvle, rvte, v1, v2 = (mm.beta1,
                                          mm.beta2,
                                          mm.alpha1,
                                          mm.alpha2,
                                          mm.w1,
                                          mm.w2)
    elif blade == 'OSV' or 'Stator':
        avle, avte, rvle, rvte, v1, v2 = (mm.alpha2,
                                          mm.alpha1,
                                          mm.beta2,
                                          mm.beta1,
                                          mm.c2,
                                          mm.c1)
    return avle, avte, rvle, rvte, v1, v2, radius


def GetTipFlowVars(stage, blade: str, station: str):
    """Obtain flow variables from stage class for tip blade locations."""
    tt = stage.tip
    radius = tt.radius
    if blade == 'IGV':
        avle, avte, rvle, rvte, v1, v2 = (0,
                                          tt.alpha1,
                                          tt.beta2,
                                          tt.beta1,
                                          tt.cx,
                                          tt.c1)
    elif blade == 'Rotor':
        avle, avte, rvle, rvte, v1, v2 = (tt.beta1,
                                          tt.beta2,
                                          tt.alpha1,
                                          tt.alpha2,
                                          tt.w1,
                                          tt.w2)
    elif blade == 'OSV' or 'Stator':
        avle, avte, rvle, rvte, v1, v2 = (tt.alpha2,
                                          tt.alpha1,
                                          tt.beta2,
                                          tt.beta1,
                                          tt.c2,
                                          tt.c1)
    return avle, avte, rvle, rvte, v1, v2, radius


def Chord(radius, hub, z, aoa):
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


def Incidence(thickness, solidity, relative_inlet_angle):
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


def Deviation(thickness, solidity, relative_inlet_angle):
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


def CamberCurvature(xc, yc):
    """Calculate camber curvature based on modified NACA 4."""
    kc = -yc*(1/(xc**2) + 1/((1 - xc)**2))

    return kc


def ThicknessCurvature(xt, yt):
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


def Curvatures(xc, yc, xt, yt):
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


def CamberBezier(xc, yc, kc, cle, cte, rle):
    """Calculate camber bezier variable in BP3333 airfoil."""
    bc = 0
    co = Cot(r(cle)) + Cot(r(cte))
    bcu = (16 + 3*kc*co + 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)
    bcl = (16 + 3*kc*co - 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)

    cbounds = (bcl, bcu)
    bci_list = []

    for bci in cbounds:
        bci_list.append(round(bci, 4))
        if bci > 0 and bci < yc:
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


def ThicknessBezier(xt, yt, kt, rle, blade: str, station: str):
    """Calculate thickness bezier variable in BP3333 airfoil."""
    bt = 0
    poly = ((27 / 4) * kt**2,
            -27 * kt**2 * xt,
            9 * kt * yt + (81 / 2) * kt**2 * xt**2,
            2 * -rle - 18 * kt * xt * yt - 27 * kt**2 * xt**3,
            3 * yt**2 + 9 * kt * xt**2 * yt + (27 / 4) * kt**2 * xt**4)

    real_troots = [root.real for root in np.roots(poly) if root.imag == 0]
    min = max(0, xt - m.sqrt((-2/3)*(yt/kt)))
    bti_list = []

    for bti in real_troots:
        bti_list.append(round(bti, 4))
        if not real_troots:
            print(f'xt = {xt}')
            print(f'yt = {yt}')
            print(f'kt = {kt:.4f}')
            raise ValueError('Invalid Airfoil Shape.'
                             f'{blade} {station} bt is complex.')

        if bti > max(0, xt - m.sqrt((-2/3)*(yt/kt))) and bti < xt:
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
        raise ValueError(f'No bt found within bounds {min:.4f} < bt < {xt}.')


def Beziers(xc, yc, kc, xt, yt, kt, cle, cte, rle,  # noqa R701
            blade: str, station: str):
    """Calculate camber and thickness Bezier variables in BP3333 airfoil."""
    bc = 0
    bt = 0

    co = Cot(r(cle)) + Cot(r(cte))
    bcu = (16 + 3*kc*co + 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)
    bcl = (16 + 3*kc*co - 4*m.sqrt(16 + 6*kc*co*(1 - yc*co)))/(3*kc*co**2)

    cbounds = (bcl, bcu)
    bci_list = []

    for bci in cbounds:
        bci_list.append(round(bci, 4))
        if bci > 0 and bci < yc:
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

    poly = ((27 / 4) * kt**2,
            -27 * kt**2 * xt,
            9 * kt * yt + (81 / 2) * kt**2 * xt**2,
            2 * -rle - 18 * kt * xt * yt - 27 * kt**2 * xt**3,
            3 * yt**2 + 9 * kt * xt**2 * yt + (27 / 4) * kt**2 * xt**4)

    real_troots = [root.real for root in np.roots(poly) if root.imag == 0]
    min = max(0, xt - m.sqrt((-2/3)*(yt/kt)))
    bti_list = []

    for bti in real_troots:
        bti_list.append(round(bti, 4))
        if not real_troots:
            print(f'xt = {xt}')
            print(f'yt = {yt}')
            print(f'kt = {kt:.4f}')
            raise ValueError('Invalid Airfoil Shape.'
                             f'{blade} {station} bt is complex.')

        if bti > max(0, xt - m.sqrt((-2/3)*(yt/kt))) and bti < xt:
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
        raise ValueError(f'No bt found within bounds {min:.4f} < bt < {xt}.')

    return bc, bt
