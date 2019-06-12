"""Writes variables to files."""

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

import csv
import json
import math as m
import os
import pathlib as pth
from math import degrees as d
from math import radians as r

import matplotlib.collections as mc
import matplotlib.pyplot as plt
import numpy as np

import utils as ut

# import pylab as pl - used in line 137 only... trying plt instead


project = os.getcwd()


def create_coord_file(filename: str, plot: bool,
                      xc, yc, kc, bc,
                      xt, yt, kt, bt,
                      cle, cte, wte, rle=0):  # pylint: disable=W0613
    """Generate airfoil coordinate file."""
    bp = {
        'x': {'LET': [0,
                      0,
                      bt,
                      xt],
              'TET': [xt,
                      2*xt - bt,
                      1 + (0 - (1.5*kt*(xt - bt)**2 + yt))*ut.cot(r(wte)),
                      1],
              'LEC': [0,
                      bc*ut.cot(r(cle)),
                      xc - ((2*(bc - yc))/(3*kc))**0.5,
                      xc],
              'TEC': [xc,
                      xc + ((2*(bc - yc))/(3*kc))**0.5,
                      1 + (0 - bc)*ut.cot(r(cte)),
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

    # x, y coordinates and slopes for LE and TE of T and C bezier curves
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
                 and 0 <= root <= (1 + 1/(20*n_points))]
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
    os.chdir('Input')
    with open(filename, 'w') as datafile:
        for i, _ in enumerate(xcoord):
            datafile.write(f'     {xcoord[i]:.6f}    {ycoord[i]:.6f}\n')
    os.chdir('..')

    if plot:
        lc = mc.LineCollection(lines)
        _, ax = plt.subplots()
        ax.add_collection(lc)
        plt.scatter(c['x'], c['yC'], color='red', marker='^')
        plt.scatter(c['x'], c['yT'], color='blue', marker='^')
        plt.scatter(xu, yu, color='green')
        plt.scatter(xl, yl, color='green')
        plt.axis('equal')
        plt.show()


def write_sw_config(stage, blade, filename):
    """Write Solidworks equation file."""
    keys = ['Blades', 'ID', 'OD']
    stations = ['root', 'mean', 'tip']
    for s in stations:
        keys.extend(['Chord_' + s, 'Stagger_' + s,
                     'xt_' + s, 'yt_' + s, 'kt_' + s, 'bt_' + s,
                     'xc_' + s, 'yc_' + s, 'kc_' + s, 'bc_' + s,
                     'cle_' + s, 'cte_' + s, 'rle_' + s, 'wte_' + s])

    sw_v = [blade.z, stage.root.radius*2/0.0254, stage.tip.radius*2/0.0254,
            blade.root.chord/0.0254, blade.root.stagger,
            blade.root.xt, blade.root.yt, blade.root.kt, blade.root.bt,
            blade.root.xc, blade.root.yc, blade.root.kc, blade.root.bc,
            blade.root.cle, blade.root.cte, blade.root.rle, blade.root.wte,
            blade.mean.chord/0.0254, blade.mean.stagger,
            blade.mean.xt, blade.mean.yt, blade.mean.kt, blade.mean.bt,
            blade.mean.xc, blade.mean.yc, blade.mean.kc, blade.mean.bc,
            blade.mean.cle, blade.mean.cte, blade.mean.rle, blade.mean.wte,
            blade.tip.chord/0.0254, blade.tip.stagger,
            blade.tip.xt, blade.tip.yt, blade.tip.kt, blade.tip.bt,
            blade.tip.xc, blade.tip.yc, blade.tip.kc, blade.tip.bc,
            blade.tip.cle, blade.tip.cte, blade.tip.rle, blade.tip.wte]
    blade_dict = dict(zip(keys, np.round(sw_v, 4)))

    os.chdir(pth.Path('../Solidworks/Vars/'))
    with open(filename, 'w') as sw_config:
        output = []
        for key, value in blade_dict.items():
            line = f'"{key}"= {value}'
            output.append(line)
        output = '\n'.join(output)
        sw_config.write(output)
    os.chdir(pth.Path('../../python'))


class DumpVars():
    """Writes all instantiated variables to csv and json dump files."""

    def __init__(self):
        """Instantiate dump properties."""

    @staticmethod
    def dump_blade_csv(filename, rotor, stator):
        """Write blade variables in CSV format."""
        os.chdir('Output')
        with open(filename, 'w', newline='') as csv_dump:
            blade_vals = {'rotor': {'root': [],
                                    'mean': [],
                                    'tip': []},
                          'stator': {'root': [],
                                     'mean': [],
                                     'tip': []}}
            blade_keys = []
            blades = ['rotor', 'stator']
            stations = ['root', 'mean', 'tip']
            blade_dicts = {'rotor': {'root': rotor.root.__dict__,
                                     'mean': rotor.mean.__dict__,
                                     'tip': rotor.tip.__dict__},
                           'stator': {'root': stator.root.__dict__,
                                      'mean': stator.mean.__dict__,
                                      'tip': stator.tip.__dict__}}
            blade_polar_dicts = {'rotor': {'root': rotor.root.polar,
                                           'mean': rotor.mean.polar,
                                           'tip': rotor.tip.polar},
                                 'stator': {'root': stator.root.polar,
                                            'mean': stator.mean.polar,
                                            'tip': stator.tip.polar}}
            for key, val in rotor.root.__dict__.items():
                blade_keys.append(key)
            blade_keys.pop(-6)
            for key, val in rotor.root.polar.items():
                blade_keys.append(key)
            data_by_row = [blade_keys]
            for bld in blades:
                for stn in stations:
                    for key, val in blade_dicts[bld][stn].items():
                        blade_vals[bld][stn].append(val)
                    blade_vals[bld][stn].pop(-6)
                    for key, val in blade_polar_dicts[bld][stn].items():
                        blade_vals[bld][stn].append(val)
                    data_by_row.append(blade_vals[bld][stn])
            data_by_column = list(zip(*data_by_row))
            fieldnames = ['DATA',
                          'rotor_root', 'rotor_mean', 'rotor_tip',
                          'stator_root', 'stator_mean', 'stator_tip']
            writer = csv.writer(csv_dump, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL)
            writer.writerow(fieldnames)
            for field in data_by_column:
                writer.writerow(field)
        os.chdir('..')

    @staticmethod
    def dump_stage_csv(filename, stage):
        """Write stage variables in CSV format."""
        os.chdir('Output')
        with open(filename, 'w', newline='') as csv_dump:
            stage_vals = {'root': [], 'mean': [], 'tip': []}
            stage_keys = []
            stations = ['root', 'mean', 'tip']
            stage_dicts = {'root': stage.root.__dict__,
                           'mean': stage.mean.__dict__,
                           'tip': stage.tip.__dict__}

            for key, val in stage.root.__dict__.items():
                stage_keys.append(key)
            data_by_row = [stage_keys]
            for stn in stations:
                for key, val in stage_dicts[stn].items():
                    val = round(val, 4)
                    stage_vals[stn].append(val)
                data_by_row.append(stage_vals[stn])
            fieldnames = ['DATA', 'stage_root', 'stage_mean', 'stage_tip']
            writer = csv.writer(csv_dump, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL)
            writer.writerow(fieldnames)
            for field in list(zip(*data_by_row)):
                writer.writerow(field)
        os.chdir('..')

    @staticmethod
    def dump_json(filename, stage, igv, rotor, stator):
        """Write variables in JSON format."""
        os.chdir('Output')
        with open(filename, 'w') as json_dump:
            json.dump(json.loads(json.dumps([[stage.root.__dict__,
                                              stage.mean.__dict__,
                                              stage.tip.__dict__],
                                             [igv.root.__dict__,
                                              igv.mean.__dict__,
                                              igv.tip.__dict__],
                                             [rotor.root.__dict__,
                                              rotor.mean.__dict__,
                                              rotor.tip.__dict__],
                                             [stator.root.__dict__,
                                              stator.mean.__dict__,
                                              stator.tip.__dict__]]),
                                 parse_float=lambda x: round(float(x), 3)),
                      json_dump, indent=2)
        os.chdir('..')
