"""Determines compressor blade geometry."""
import math as m
import pathlib as pth
from configparser import ConfigParser
from math import degrees as d
from math import radians as r

import matplotlib.pyplot as plt

import utils as u
import write_file as wf
from aeropy import xfoil_module as xf

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


class Airfoil():
    """Set Bezier-PARSEC variables."""

    def __init__(self):
        """Instantiate airfoil properties."""
        self.z = 0
        self.hub = 0
        self.ar = 0
        self.df = 0
        self.dh = 0

        self.sigma = 0
        self.space = 0
        self.chord = 0
        self.deflection = 0
        self.incidence = 0
        self.deviation = 0
        self.blade_angle_i = 0
        self.blade_angle_e = 0
        self.camber = 0
        self.blade_angle_e2 = 0
        self.camber2 = 0
        self.stagger = 0
        self.aoa = 0

        self.xcr = 0
        self.ycr = 0
        self.xc = 0
        self.yc = 0
        self.kc = 0
        self.bc = 0
        self.xt = 0
        self.yt = 0
        self.kt = 0
        self.bt = 0
        self.cle = 0
        self.cte = 0
        self.rle = 0
        self.wte = 0

        self.rho = 1.2
        self.mu = 0.00018
        self.Re = 0
        self._psi = 0
        self.psi_opt = 0
        self.polar = {}

        self.dX = 0
        self.dTau = 0
        self.delta_p = 0
        self.delta_T = 0
        self.efficiency = 0

    def GetAirfoilConfig(self, filename: str, station: str):
        """Read .ini file."""
        cfp = ConfigParser()
        file = pth.Path(f'Config/{filename}')
        cfp.read(file)

        self.z = cfp.getfloat('blade', 'z')
        self.hub = cfp.getfloat('blade', 'hub')
        self.ar = cfp.getfloat('blade', 'ar')

        self.df = cfp.getfloat(station, 'DF')
        self.xcr = cfp.getfloat(station, 'xcr')
        self.ycr = cfp.getfloat(station, 'ycr')
        self.xt = cfp.getfloat(station, 'xt')
        self.yt = cfp.getfloat(station, 'yt')

    def CalcAirfoil(self, stage, blade: str, station: str,
                    plotAirfoil=False, plotPolar=False,
                    plotCp=False, plotSv=False):
        """Calculate Bezier-PARSEC variables and airfoil properties."""
        if station == 'root':
            FlowVars = u.GetRootFlowVars(stage=stage,
                                         blade=blade,
                                         station=station)
        elif station == 'mean':
            FlowVars = u.GetMeanFlowVars(stage=stage,
                                         blade=blade,
                                         station=station)
        elif station == 'tip':
            FlowVars = u.GetTipFlowVars(stage=stage,
                                        blade=blade,
                                        station=station)
        avle, avte, rvle, rvte, v1, v2, radius = FlowVars
        cx = v1*m.cos(r(avle))
        U = stage.rpm/60*2*m.pi*radius
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
        self.incidence, n_slope = u.Incidence(thickness=2*self.yt,
                                              solidity=self.sigma,
                                              relative_inlet_angle=rvle)
        self.deviation, m_slope = u.Deviation(thickness=2*self.yt,
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
        self.kc, self.kt = u.Curvatures(xc=self.xc, yc=self.yc,
                                        xt=self.xt, yt=self.yt)

        self.bc, self.bt = u.Beziers(xc=self.xc, yc=self.yc, kc=self.kc,
                                     xt=self.xt, yt=self.yt, kt=self.kt,
                                     cle=self.cle, cte=self.cte, rle=self.rle,
                                     blade=blade, station=station)
        airfoil_name = blade + '_' + station
        wf.Airfoil.CreateFile(filename=airfoil_name, plot=plotAirfoil,
                              xc=self.xc, yc=self.yc, kc=self.kc, bc=self.bc,
                              xt=self.xt, yt=self.yt, kt=self.kt, bt=self.bt,
                              cle=self.cle, cte=self.cte,
                              rle=self.rle, wte=self.wte)

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

        if plotCp:
            plt.plot(self.cp['x'], self.cp['Cp'])
            plt.show()

        if plotSv:
            plt.plot(self.sv['x'], self.sv['v'])
            plt.show()

        if plotPolar:
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
            self._psi = (((cx/U)/2)*u.Sec(r(betam))*self.sigma
                         * (cl + cd*m.tan(r(betam))))
            self.psi_opt = (((cx/U)/m.sqrt(2))*self.sigma*(cl + cd))

            self.dX = ((self.rho*cx**2*self.chord*cl
                        / (2*m.cos(r(betam))**2))
                       * m.sin(r(betam) - gamma)/m.cos(gamma))
            self.dTau = ((self.rho*cx**2*self.chord*cl
                          * self.z*radius/(2*m.cos(r(betam))**2))
                         * m.cos(r(betam) - gamma)/m.cos(gamma))
            self.delta_p = cl*((self.rho*cx**2*self.chord
                                / (2*self.space*m.cos(r(betam))**2))
                               * m.sin(r(betam) - gamma)/m.cos(gamma))
            self.delta_T = (cl/1.005)*((U*cx*self.chord
                                        / (2*self.space*m.cos(r(betam))**2))
                                       * m.cos(r(betam) - gamma)/m.cos(gamma))
            self.efficiency = (cx/U)*m.tan(r(betam) - gamma) + ct2/(2*U)
        except TypeError:
            print('Design point did not converge in xfoil.')
            pass

        print(f'{station} {blade} Performance')
        print(f'    DF: {self.df:.2f} (<=0.6)')
        print(f'    DH: {self.dh:.2f} (>=0.72)')
        print(f'    i : {self.incidence:.2f} deg')
        print(f'    d : {self.deviation:.2f} deg')
        print('\n')


class Blade():
    """Set blade properties."""

    def __init__(self):
        """Instantiate blade properties."""
        self.root = Airfoil()
        self.mean = Airfoil()
        self.tip = Airfoil()

        self.z = 0

    def GetBladeConfig(self, filename: str):
        """Read .ini files."""
        self.root.GetAirfoilConfig(filename, 'root')
        self.mean.GetAirfoilConfig(filename, 'mean')
        self.tip.GetAirfoilConfig(filename, 'tip')

        cfp = ConfigParser()
        file = pth.Path(f'Config/{filename}')
        cfp.read(file)

        self.z = cfp.getfloat('blade', 'z')
        self.hub = cfp.getfloat('blade', 'hub')

    def CalcBlade(self, stage, blade: str,
                  plotAirfoil=False, plotPolar=False,
                  plotCp=False, plotSv=False):
        """Calculate blade properties."""
        self.root.CalcAirfoil(stage, blade, 'root',
                              plotAirfoil=plotAirfoil, plotPolar=plotPolar,
                              plotCp=plotCp, plotSv=plotSv)
        self.mean.CalcAirfoil(stage, blade, 'mean',
                              plotAirfoil=plotAirfoil, plotPolar=plotPolar,
                              plotCp=plotCp, plotSv=plotSv)
        self.tip.CalcAirfoil(stage, blade, 'tip',
                             plotAirfoil=plotAirfoil, plotPolar=plotPolar,
                             plotCp=plotCp, plotSv=plotSv)
