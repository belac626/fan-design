"""Airfoil calculation functions."""
import math as m
from math import degrees as d
from math import radians as r

import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from matplotlib import collections as mc


def Cot(rad):
    """Return Cotangent of an angle in radians."""
    return m.cos(rad)/m.sin(rad)


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
    # location = {
    #     ('IGV', 'root'): (0,
    #                       stage.root.alpha1,
    #                       stage.root.cx,
    #                       stage.root.c1,
    #                       stage.root.radius),
    #     ('IGV', 'mean'): (0,
    #                       stage.mean.alpha1,
    #                       stage.mean.cx,
    #                       stage.mean.c1,
    #                       stage.mean.radius),
    #     ('IGV', 'tip'): (0,
    #                       stage.tip.alpha1,
    #                       stage.tip.cx,
    #                       stage.tip.c1,
    #                       stage.tip.radius),
    #     ('Rotor', 'root'): (stage.root.beta1,
    #                         stage.root.beta2,
    #                         stage.root.w1,
    #                         stage.root.w2,
    #                         stage.root.radius),
    #     ('Rotor', 'mean'): (stage.mean.beta1,
    #                         stage.mean.beta2,
    #                         stage.mean.w1,
    #                         stage.mean.w2,
    #                         stage.mean.radius),
    #     ('Rotor', 'tip'): (stage.tip.beta1,
    #                        stage.tip.beta2,
    #                        stage.tip.w1,
    #                        stage.tip.w2,
    #                        stage.tip.radius),
    #     ('Stator', 'root'): (stage.root.alpha2,
    #                          stage.root.alpha1,
    #                          stage.root.c2,
    #                          stage.root.c1,
    #                          stage.root.radius),
    #     ('Stator', 'mean'): (stage.mean.alpha2,
    #                          stage.mean.alpha1,
    #                          stage.mean.c2,
    #                          stage.mean.c1,
    #                          stage.mean.radius),
    #     ('Stator', 'tip'): (stage.tip.alpha2,
    #                         stage.tip.alpha1,
    #                         stage.tip.c2,
    #                         stage.tip.c1,
    #                         stage.tip.radius),
    #     ('OSV', 'root'): (stage.root.alpha2,
    #                       stage.root.alpha1,
    #                       stage.root.c2,
    #                       stage.root.c1,
    #                       stage.root.radius),
    #     ('OSV', 'mean'): (stage.mean.alpha2,
    #                       stage.mean.alpha1,
    #                       stage.mean.c2,
    #                       stage.mean.c1,
    #                       stage.mean.radius),
    #     ('OSV', 'tip'): (stage.tip.alpha2,
    #                      stage.tip.alpha1,
    #                      stage.tip.c2,
    #                      stage.tip.c1,
    #                      stage.tip.radius),
    #     }
    # avle, avte, rvle, rvte, v1, v2, radius = location[(blade, station)]


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
                print('xc = {xc}'.format(xc=xc))
                print('yc = {yc}'.format(yc=yc))
                print('kc = {kc}'.format(kc=round(kc, 4)))
                print('bc = {bci}'.format(bci=round(bci, 4)))
                print('x - sqrt((2/3)(b - y)/k) = {xlc2}'.format(xlc2=xlc2))
                raise ValueError('x - sqrt((2/3)(b - y)/k) = {xlc2}'
                                 .format(xlc2=xlc2) + ' Invalid Airfoil Shape')

    if bc == 0:
        print('xc = {xc}'.format(xc=xc))
        print('yc = {yc}'.format(yc=yc))
        print('kc = {kc}'.format(kc=round(kc, 4)))
        print('bc = {bci}'.format(bci=bci_list))
        raise ValueError('No bc found within bounds 0 < bc < {yc}. \n\
              '.format(yc=round(yc, 4)))


def ThicknessBezier(xt, yt, kt, rle, blade: str, station: str):
    """Calculate thickness bezier variable in BP3333 airfoil."""
    bt = 0
    poly = ((27 / 4) * kt**2,
            -27 * kt**2 * xt,
            9 * kt * yt + (81 / 2) * kt**2 * xt**2,
            2 * -rle - 18 * kt * xt * yt - 27 * kt**2 * xt**3,
            3 * yt**2 + 9 * kt * xt**2 * yt + (27 / 4) * kt**2 * xt**4)

    real_troots = [root.real for root in np.roots(poly) if root.imag == 0]
    min = xt - m.sqrt((-2/3)*(yt/kt))
    bti_list = []

    for bti in real_troots:
        bti_list.append(round(bti, 4))
        if not real_troots:
            print('xt = {xt}'.format(xt=xt))
            print('yt = {yt}'.format(yt=yt))
            print('kt = {kt}'.format(kt=round(kt, 4)))
            raise ValueError('{r} {s} bt is complex. \n\
                             Invalid Airfoil Shape.'.format(r=blade,
                                                            s=station))

        if bti > max(0, xt - m.sqrt((-2/3)*(yt/kt))) and bti < xt:
            bt = float(bti)
            ylt1 = yt + (3/2)*kt*(xt - bti)**2

            if ylt1 <= 0:
                print('xt = {xt}'.format(xt=xt))
                print('yt = {yt}'.format(yt=yt))
                print('kt = {kt}'.format(kt=round(kt, 4)))
                print('bt = {bti}'.format(bti=round(bti, 4)))
                print('(y + (3/2)k(x - b)^2) = {ylt1}'.format(ylt1=ylt1))
                raise ValueError('(y + (3/2)k(x - b)^2) = {ylt1}.'
                                 .format(ylt1=ylt1) + ' Invalid Airfoil Shape')

    if bt == 0:
        print('xt = {xt}'.format(xt=xt))
        print('yt = {yt}'.format(yt=yt))
        print('kt = {kt}'.format(kt=round(kt, 4)))
        print('bt = {bti}'.format(bti=bti_list))
        raise ValueError('No bt found within bounds {min} < bt < {xt}. \n\
              '.format(min=round(max(0, min), 4),
                       xt=xt))


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
                print('xc = {xc}'.format(xc=xc))
                print('yc = {yc}'.format(yc=yc))
                print('kc = {kc}'.format(kc=round(kc, 4)))
                print('bc = {bci}'.format(bci=round(bci, 4)))
                print('x - sqrt((2/3)(b - y)/k) = {xlc2}'.format(xlc2=xlc2))
                raise ValueError('x - sqrt((2/3)(b - y)/k) = {xlc2}'
                                 .format(xlc2=xlc2) + ' Invalid Airfoil Shape')

    if bc == 0:
        print('xc = {xc}'.format(xc=xc))
        print('yc = {yc}'.format(yc=yc))
        print('kc = {kc}'.format(kc=round(kc, 4)))
        print('bc = {bci}'.format(bci=bci_list))
        raise ValueError('No bc found within bounds 0 < bc < {yc}. \n\
              '.format(yc=round(yc, 4)))

    poly = ((27 / 4) * kt**2,
            -27 * kt**2 * xt,
            9 * kt * yt + (81 / 2) * kt**2 * xt**2,
            2 * -rle - 18 * kt * xt * yt - 27 * kt**2 * xt**3,
            3 * yt**2 + 9 * kt * xt**2 * yt + (27 / 4) * kt**2 * xt**4)

    real_troots = [root.real for root in np.roots(poly) if root.imag == 0]
    min = xt - m.sqrt((-2/3)*(yt/kt))
    bti_list = []

    for bti in real_troots:
        bti_list.append(round(bti, 4))
        if not real_troots:
            print('xt = {xt}'.format(xt=xt))
            print('yt = {yt}'.format(yt=yt))
            print('kt = {kt}'.format(kt=round(kt, 4)))
            raise ValueError('{r} {s} bt is complex. \n\
                             Invalid Airfoil Shape.'.format(r=blade,
                                                            s=station))

        if bti > max(0, xt - m.sqrt((-2/3)*(yt/kt))) and bti < xt:
            bt = float(bti)
            ylt1 = yt + (3/2)*kt*(xt - bti)**2

            if ylt1 <= 0:
                print('xt = {xt}'.format(xt=xt))
                print('yt = {yt}'.format(yt=yt))
                print('kt = {kt}'.format(kt=round(kt, 4)))
                print('bt = {bti}'.format(bti=round(bti, 4)))
                print('(y + (3/2)k(x - b)^2) = {ylt1}'.format(ylt1=ylt1))
                raise ValueError('(y + (3/2)k(x - b)^2) = {ylt1}.'
                                 .format(ylt1=ylt1) + ' Invalid Airfoil Shape')

    if bt == 0:
        print('xt = {xt}'.format(xt=xt))
        print('yt = {yt}'.format(yt=yt))
        print('kt = {kt}'.format(kt=round(kt, 4)))
        print('bt = {bti}'.format(bti=bti_list))
        raise ValueError('No bt found within bounds {min} < bt < {xt}. \n\
              '.format(min=round(max(0, min), 4),
                       xt=xt))

    return bc, bt


def CreateAirfoilFile(filename: str, plot: bool,  # noqa R701
                      xc, yc, kc, bc,
                      xt, yt, kt, bt,
                      cle, cte, rle, wte):
    """Generate airfoil coordinate file.

    Depends on aeropy.xfoil_module.
    """
    bp = {
        'x': {'LET': [0,
                      0,
                      bt,
                      xt],
              'TET': [xt,
                      2*xt - bt,
                      1 + (0 - (1.5*kt*(xt - bt)**2 + yt))*Cot(r(wte)),
                      1],
              'LEC': [0,
                      bc*Cot(r(cle)),
                      xc - ((2*(bc - yc))/(3*kc))**0.5,
                      xc],
              'TEC': [xc,
                      xc + ((2*(bc - yc))/(3*kc))**0.5,
                      1 + (0 - bc)*Cot(r(cte)),
                      1]},
        'y': {'LET': [0,
                      (1.5*kt*(xt - bt)**2 + yt),
                      yt,
                      yt],
              'TET': [yt,
                      yt,
                      (1.5*kt*(xt - bt)**2 + yt),
                      0],
              'LEC': [0,
                      bc,
                      yc,
                      yc],
              'TEC': [yc,
                      yc,
                      bc,
                      0]}
    }
    c = {'x': [], 'yC': [], 'yT': []}
    dc = {'xC': [], 'xT': [], 'yC': [], 'yT': []}

    # Get x, y coordinates and slopes for LE and TE of T and C bezier curves
    n_points = 80
    for i in range(0, n_points + 1):
        x = (1 - m.cos(i*m.pi/n_points))/2
        c['x'].append(x)
        for j in ['C', 'T']:
            if j == 'C':
                xmax = xc
            elif j == 'T':
                xmax = xt
            if x <= xmax:
                loc = 'LE'
                bp_xi = bp['x'][loc + j]
            elif x > xmax:
                loc = 'TE'
                bp_xi = bp['x'][loc + j]

            u_x = [-bp_xi[0] + 3*bp_xi[1] - 3*bp_xi[2] + bp_xi[3],
                   3*bp_xi[0] - 6*bp_xi[1] + 3*bp_xi[2],
                   -3*bp_xi[0] + 3*bp_xi[1],
                   bp_xi[0] - x]
            u = [root.real for root in np.roots(u_x)
                 if root.imag == 0
                 and root >= 0
                 and root <= (1 + 1/(20*n_points))]
            u = u[0]
            c['y' + j].append(bp['y'][loc + j][0]*(1 - u)**3
                              + 3*bp['y'][loc + j][1]*(u)*(1 - u)**2
                              + 3*bp['y'][loc + j][2]*(u)**2*(1 - u)
                              + bp['y'][loc + j][3]*(u)**3)
            for k in ['x', 'y']:
                dc[k + j].append(bp[k][loc + j][0]*(-3*(1 - u)**2)
                                 + 3*bp[k][loc + j][1]*(3*u**2 - 4*u + 1)
                                 + 3*bp[k][loc + j][2]*(-3*u**2 + 2*u)
                                 + bp[k][loc + j][3]*(3*u**2))

    theta = []
    xu = []
    yu = []
    xl = []
    yl = []
    lines = []
    for i in range(0, len(c['x'])):  # theta[0] is a divide by 0 error
        theta.append(d(m.atan(dc['yC'][i]/dc['xC'][i])))
        xu.append(c['x'][i] - c['yT'][i]*m.sin(r(theta[i])))
        yu.append(c['yC'][i] + c['yT'][i]*m.cos(r(theta[i])))
        xl.append(c['x'][i] + c['yT'][i]*m.sin(r(theta[i])))
        yl.append(c['yC'][i] - c['yT'][i]*m.cos(r(theta[i])))
        lines.append([(xu[i], yu[i]), (xl[i], yl[i])])

    # Organize coordinates so that xfoil can read them
    xcoord = list(reversed(xu))[:-1] + xl
    ycoord = list(reversed(yu))[:-1] + yl
    # Create airfoil file for xfoil
    with open(filename, 'w') as DataFile:
        for i in range(len(xcoord)):
            DataFile.write(f'     {xcoord[i]:.6f}    {ycoord[i]:.6f}\n')
    if plot:
        lc = mc.LineCollection(lines)
        fig, ax = pl.subplots()
        ax.add_collection(lc)
        plt.scatter(c['x'], c['yC'], color='red', marker='^')
        plt.scatter(c['x'], c['yT'], color='blue', marker='^')
        plt.scatter(xu, yu, color='green')
        plt.scatter(xl, yl, color='green')
        plt.axis('equal')
        plt.show()

    # Prepare input for aeropy.xfoil_module.prepare_xfoil
    # upper = {'x': xu, 'y': yu}
    # lower = {'x': xl, 'y': yl}
    # coord = xf.prepare_xfoil(upper, lower, chord=1, reposition=False)
    # xcoord = [coord[i][0] for i in range(len(coord))]
    # ycoord = [coord[i][1] for i in range(len(coord))]
    # xf.create_input(xcoord, ycoord, filename=filename,
    #                 different_x_upper_lower=True)
