"""Determines compressor blade geometry."""
import math as m
import pathlib as pth
from configparser import ConfigParser
from math import degrees as d
from math import radians as r
from typing import ClassVar, Dict

import matplotlib.pyplot as plt

import utils as ut
import write_file as wf
from aeropy import xfoil_module as xf
from dataclasses import dataclass

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

################################
# #Class: Airfoil
# Gets and calculates variables needed to parameterize BP3333 airfoil
# #Attributes: (all linear dimensions normalized to chord length)

# z: number of blades (float)
# hub: thickness of blade row (float) (m)
# ar: blade aspect ratio (float)
# df: Leibein's Diffusion Factor (incompressible) (float)
# dh: DeHaller Number (float)

# sigma: solidity (chord/space) (float)
# space: blade pitch (float) (m)
# chord: blade chord (float) (m)
# deflection: flow turning angle (float) (deg)
# incidence: difference between flow and blade angle at inlet (float) (deg)
# deviation: difference between flow and blade angle at exit (float) (deg)
# blade_angle_i: inlet blade angle (float) (deg) (reference to axial)
# blade_angle_e: exit blade angle (float) (deg) (reference to axial)
# camber: blade turning angle (float) (deg) (calculated from camber slope)
# blade_angle_e2: exit blade angle (float) (deg) (reference to axial)
# camber2: blade turning angle (float) (deg) (calculated from deviation)
# stagger: stagger angle (float) (deg) (used to transform angles to sw coord)
# aoa: angle of attack of airfoil (float) (deg) (used in xfoil)

# xcr: stagger ratio (also very nearly xc) (1 = avle, 0 = avte) (float)
# ycr: max camber altitude ratio (1 = triangular camber, 0 = no camber) (float)
# xc: location of max camber (float)
# yc: max camber ratio (float)
# kc: curvature of camber crest (float)
# bc: bezier point from camber (float)
# xt: location of max thickness (float)
# yt: max thickness ratio (float)
# kt: curvature of thickness crest (float)
# bt: bezier point from thickness (float)
# cle: leading edge camber angle (for solidworks, from chord) (float) (deg)
# cte: trailing edge camber angle (for solidworks, from chord) (float) (deg)
# rle: leading edge radius (float) (NACA definition)
# wte: trailing edge wedge angle (float) (deg) (NACA definition)

# rho: density of air under standard conditions (float) (kg/m^3)
# mu: dynamic viscosity of air under standard conditions (float) (kg/(m*s))
# Re: airfoil reynolds number (float)
# _psi: loading coefficient from lift coefficients (float)
# psi_opt: optimum loading coefficient from fan theory (float)
# polar: airfoil coefficients solved for by xfoil (lift, drag, moment, etc.)

# dx: elemental thrust per blade (float) (N)
# dtau: elemental ring torque (float) (Nm)
# delta_p: pressure change (float) (Pa)
# delta_T: temperature change (float) (K)
# efficiency: blade efficiency (float)
################################


@dataclass
class Airfoil:
    """Set Bezier-PARSEC variables."""

    polar: Dict = 0
    cp: Dict = 0
    sv: Dict = 0

    z: float = 0
    hub: float = 0
    ar: float = 0
    df: float = 0
    dh: float = 0

    sigma: float = 0
    space: float = 0
    chord: float = 0
    deflection: float = 0
    incidence: float = 0
    deviation: float = 0
    blade_angle_i: float = 0
    blade_angle_e: float = 0
    camber: float = 0
    blade_angle_e2: float = 0
    camber2: float = 0
    stagger: float = 0
    aoa: float = 0

    xcr: float = 0
    ycr: float = 0
    xc: float = 0
    yc: float = 0
    kc: float = 0
    bc: float = 0
    xt: float = 0
    yt: float = 0
    kt: float = 0
    bt: float = 0
    cle: float = 0
    cte: float = 0
    rle: float = 0
    wte: float = 0

    rho: float = 1.2
    mu: float = 0.00018
    Re: float = 0
    _psi: float = 0
    psi_opt: float = 0

    dX: float = 0
    dTau: float = 0
    delta_p: float = 0
    delta_T: float = 0
    efficiency: float = 0

    def get_airfoil_config(self, filename: str, station: str):
        """Read .ini file."""
        file = pth.Path(f'Config/{filename}')
        config = ConfigParser()
        config.read(file)

        self.z = config.getfloat('blade', 'z')
        self.hub = config.getfloat('blade', 'hub')
        self.ar = config.getfloat('blade', 'ar')

        self.df = config.getfloat(station, 'DF')
        self.xcr = config.getfloat(station, 'xcr')
        self.ycr = config.getfloat(station, 'ycr')
        self.xt = config.getfloat(station, 'xt')
        self.yt = config.getfloat(station, 'yt')

    def calcairfoil(self, stage, blade: str, station: str,
                    plot_airfoil=False, plot_ploar=False,
                    plot_cp=False, plot_sv=False):
        """Calculate Bezier-PARSEC variables and airfoil properties."""
        avle, avte, rvle, rvte, v1, v2, radius = ut.get_vars(stage=stage,
                                                             blade=blade,
                                                             station=station)
        cx = v1*m.cos(r(avle))
        u = stage.rpm/60*2*m.pi*radius
        ct2 = cx*m.tan(r(rvte))
        betam = d(m.atan((m.tan(r(avle)) + m.tan(r(avte)))/2))

        self.dh = v2/v1
        self.space = (2*m.pi*radius)/self.z
        # self.chord = stage.span/self.ar
        # self.sigma = self.chord/self.space
        # self.df = ((1 - m.cos(r(avle))/m.cos(r(avte)))
        #            + ((m.cos(r(avle))*(m.tan(r(avle)) - m.tan(r(avte))))
        #               / (2*self.sigma)))
        self.sigma = ((m.cos(r(avle))*(m.tan(r(avle)) - m.tan(r(avte))))
                      / (2*(m.cos(r(avle))/m.cos(r(avte)) - 1 + self.df)))
        self.chord = self.sigma*self.space

        # Incidence, deviation, and blade angle calculation
        self.deflection = abs(avle - avte)
        self.incidence, _n_slope = ut.calcincidence(thickness=2*self.yt,
                                                    solidity=self.sigma,
                                                    relative_inlet_angle=rvle)
        self.deviation, _m_slope = ut.calcdeviation(thickness=2*self.yt,
                                                    solidity=self.sigma,
                                                    relative_inlet_angle=rvle)
        self.blade_angle_i = avle - self.incidence
        self.blade_angle_e = avte - self.deviation
        self.camber = abs(self.blade_angle_i - self.blade_angle_e)
        # Compare to:
        # self.camber2 = (self.deflection
        #                 - (self.incidence
        #                    - self.deviation))/(1 - m_slope + n_slope)
        # self.blade_angle_e2 = self.blade_angle_i + self.camber2

        # Assumes double circular arc airfoil to estimate stagger
        self.stagger = (self.blade_angle_i*self.xcr
                        + self.blade_angle_e*(1 - self.xcr))
        self.aoa = abs(self.stagger - avle)
        # Redefine blade angles to sw blade angles (stagger independent)
        self.cle = self.blade_angle_i - self.stagger
        self.cte = self.stagger - self.blade_angle_e

        # Assume NACA 4 definition of wedge angle and leading edge radius
        self.rle = 1.1019*(2*self.yt)**2
        self.wte = 2*d(m.atan(1.16925*(2*self.yt)))  # TE radius done in sw

        # Calculate natural xc,yc from cle and cte using ycr altitutde ratio
        self.xc = m.tan(r(self.cte))/(m.tan(r(self.cle)) + m.tan(r(self.cte)))
        # self.ycr = (self.xc + 1)/2
        self.yc = self.ycr*self.xc*abs(m.tan(r(self.cle)))

        # Assume NACA 4 curvature
        self.kc, self.kt = ut.curvatures(xc=self.xc, yc=self.yc,
                                         xt=self.xt, yt=self.yt)

        self.bc, self.bt = ut.beziers(xc=self.xc, yc=self.yc, kc=self.kc,
                                      xt=self.xt, yt=self.yt, kt=self.kt,
                                      cle=self.cle, cte=self.cte, rle=self.rle,
                                      blade=blade, station=station)
        airfoil_name = blade + '_' + station
        wf.create_coord_file(filename=airfoil_name, plot=plot_airfoil,
                             xc=self.xc, yc=self.yc, kc=self.kc, bc=self.bc,
                             xt=self.xt, yt=self.yt, kt=self.kt, bt=self.bt,
                             cle=self.cle, cte=self.cte,
                             wte=self.wte, rle=self.rle)

        self.Re = self.rho*v1*self.chord/self.mu
        self.polar = xf.find_coefficients(airfoil=airfoil_name, indir='Input',
                                          outdir='Output', alpha=self.aoa,
                                          Reynolds=self.Re, iteration=500,
                                          echo=False, delete=False, NACA=False,
                                          PANE=True)
        self.cp = xf.find_pressure_coefficients(airfoil=airfoil_name,
                                                alpha=self.aoa, indir='Input',
                                                outdir='Output',
                                                Reynolds=self.Re,
                                                iteration=500, echo=False,
                                                NACA=False, chord=1.,
                                                PANE=True, delete=False)
        self.sv = {'x': list(), 'y': list(), 'v': list()}
        for i in range(len(self.cp['Cp'])):
            self.sv['x'].append(self.cp['x'][i])
            # self.sv['y'].append(self.cp['y'][i])
            self.sv['v'].append(m.sqrt(1 - self.cp['Cp'][i]))

        if plot_cp:
            plt.plot(self.cp['x'], self.cp['Cp'])
            plt.show()

        if plot_sv:
            plt.plot(self.sv['x'], self.sv['v'])
            plt.show()

        if plot_ploar:
            # Gather multiple angles of attack for airfoil and get polars
            alphas = list(range(-30, 30))
            polars = xf.find_coefficients(airfoil=airfoil_name,
                                          alpha=alphas,
                                          Reynolds=self.Re,
                                          iteration=500, delete=True,
                                          echo=False,
                                          NACA=False, PANE=True)
            # Plot airfoil
            plt.plot(polars['alpha'], polars['CL'])
            plt.show()

        try:
            cl = self.polar['CL']
            cd = self.polar['CD']
            gamma = m.atan(cd/cl)
            self._psi = (((cx/u)/2)*ut.sec(r(betam))*self.sigma
                         * (cl + cd*m.tan(r(betam))))
            self.psi_opt = (((cx/u)/m.sqrt(2))*self.sigma*(cl + cd))

            self.dX = ((self.rho*cx**2*self.chord*cl
                        / (2*m.cos(r(betam))**2))
                       * m.sin(r(betam) - gamma)/m.cos(gamma))
            self.dTau = ((self.rho*cx**2*self.chord*cl
                          * self.z*radius/(2*m.cos(r(betam))**2))
                         * m.cos(r(betam) - gamma)/m.cos(gamma))
            self.delta_p = cl*((self.rho*cx**2*self.chord
                                / (2*self.space*m.cos(r(betam))**2))
                               * m.sin(r(betam) - gamma)/m.cos(gamma))
            self.delta_T = (cl/1.005)*((u*cx*self.chord
                                        / (2*self.space*m.cos(r(betam))**2))
                                       * m.cos(r(betam) - gamma)/m.cos(gamma))
            self.efficiency = (cx/u)*m.tan(r(betam) - gamma) + ct2/(2*u)
        except TypeError:
            print('Design point did not converge in xfoil.')
            # pass

        print(f'{station} {blade} Performance')
        print(f'    DF: {self.df:.2f} (<=0.6)')
        print(f'    DH: {self.dh:.2f} (>=0.72)')
        print(f'    i : {self.incidence:.2f} deg')
        print(f'    d : {self.deviation:.2f} deg')
        print('\n')


@dataclass
class Blade:
    """Set blade properties."""

    root: ClassVar = Airfoil()
    mean: ClassVar = Airfoil()
    tip: ClassVar = Airfoil()

    z: float = 0
    hub: float = 0

    def get_blade_config(self, filename: str):
        """Read .ini files."""
        self.root.get_airfoil_config(filename, 'root')
        self.mean.get_airfoil_config(filename, 'mean')
        self.tip.get_airfoil_config(filename, 'tip')

        file = pth.Path(f'Config/{filename}')
        config = ConfigParser()
        config.read(file)

        self.z = config.getfloat('blade', 'z')
        self.hub = config.getfloat('blade', 'hub')

    def calcblade(self, stage, blade: str,
                  plot_airfoil=False, plot_ploar=False,
                  plot_cp=False, plot_sv=False):
        """Calculate blade properties."""
        self.root.calcairfoil(stage, blade, 'root',
                              plot_airfoil=plot_airfoil, plot_ploar=plot_ploar,
                              plot_cp=plot_cp, plot_sv=plot_sv)
        self.mean.calcairfoil(stage, blade, 'mean',
                              plot_airfoil=plot_airfoil, plot_ploar=plot_ploar,
                              plot_cp=plot_cp, plot_sv=plot_sv)
        self.tip.calcairfoil(stage, blade, 'tip',
                             plot_airfoil=plot_airfoil, plot_ploar=plot_ploar,
                             plot_cp=plot_cp, plot_sv=plot_sv)
