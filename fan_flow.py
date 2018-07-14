"""Determines fan flow and geometric properties."""
import os
import math as m
from configparser import ConfigParser
from math import radians as r
from math import degrees as d


################################
# #Class: ThermoFlow
# Holds stage flow characteristics


# ---Flow Attributes:
# r: reaction coefficient
# phi: loading coefficient
# psi: flow coefficient
# rpm: shaft speed (rpm)
# radius: station in blade (in)

# CL: airfoil lift coefficient
# CD: airfoil drag coefficient

# u: blade speed (in/s)
# cx: axial air velocity (in/s)
# beta1: relative air inlet angle (deg)
# beta2: relative air exit angle (deg)
# betam: atan(0.5*(tan(beta1) + tan(beta2))
# alpha1: absolute air inlet angle (deg)
# alpha2: absolute air exit angle (deg)
# alpham: atan(0.5*(tan(alpha1) + tan(alpha2))
# w1: relative air inlet speed (in/s)
# w2: relative air exit speed (in/s)
# c1: absolute air inlet speed (in/s)
# c2: absolute air exit speed (in/s)
# c1: absolute air inlet speed (in/s)
# c2: absolute air exit speed (in/s)
# wt1: tangential relative air inlet speed (in/s)
# wt2: tangential relative air outlet speed (in/s)


# ---Thermodynamic Attributes
# M: axial mach number
# T1: inlet static temperature (K)
# T2: exit static temperature (K)
# T01: inlet stagnation temperature (K)
# T02: exit stagnation temperature (K)
# p1: inlet static pressure (bar)
# p2: exit static pressure
# p01: inlet stagnation pressure
# p02: exit stagnation pressure

# rho: density (air) (kg/m^3)
# cp: constant pressure specific heat (air) (J/kg*K)
# cv: constant volume specific heat (air) (J/kg*K)
# gamma: specific heat ratio (air)
# R: ideal gas constant (air) (J/kg*K)

# dT0: change in stagnation temperature
# PR: pressure ratio
# Cp: static pressure rise coefficient
# W: ideal work required
# tau: ideal torque required
################################
class FanFlow():
    """Set flow properties."""

    def __init__(self):
        """Instantiate flow properties."""
        self.r = 0
        self.phi = 0
        self.psi = 0
        self.rpm = 0
        self.radius = 0

        self.u = 0
        self.beta1 = 0
        self.beta2 = 0
        self.betam = 0
        self.alpha1 = 0
        self.alpha2 = 0
        self.alpham = 0
        self.cx = 0
        self.w1 = 0
        self.w2 = 0
        self.c1 = 0
        self.c2 = 0
        self.wt1 = 0
        self.wt2 = 0
        self.ct1 = 0
        self.ct2 = 0

        self.T1 = 0
        self.P1 = 0

    def GetFanFlowConfig(self, filename: str, station: str):
        """Read .ini file."""
        cfp = ConfigParser()

        os.chdir(r'.\Config')
        cfp.read(filename)  # Fan_Stage.ini

        self.r = cfp.getfloat('flow', 'r')
        self.phi = cfp.getfloat('flow', 'phi')
        self.rpm = cfp.getfloat('flow', 'rpm')
        self.alpha1 = cfp.getfloat('flow', 'alpha1')
        # self.P1 = cfp.getfloat('thermo', 'p1')
        # self.T1 = cfp.getfloat('thermo', 'T1')
        self.span = (cfp.getfloat('tip', 'radius')
                     - cfp.getfloat('root', 'radius'))
        if station == 'mean':
            self.radius = m.sqrt((cfp.getfloat('root', 'radius')**2
                                  + cfp.getfloat('tip', 'radius')**2)/2)
        else:
            self.radius = cfp.getfloat(station, 'radius')
        os.chdir('..')

    def CalcFanFlow(self, iphi=None):
        """Calculate flow angles."""
        # Flow
        if iphi is not None:
            self.phi = iphi
        self.psi = 2*(1 - self.r - self.phi*m.tan(r(self.alpha1)))
        self.u = self.rpm/60*2*m.pi*self.radius
        self.beta1 = d(m.atan((1/self.phi) - m.tan(r(self.alpha1))))
        self.beta2 = d(m.atan((2*self.r-1)/self.phi + m.tan(r(self.alpha1))))
        self.betam = d(m.atan((m.tan(r(self.beta1))
                               + m.tan(r(self.beta2)))/2))
        self.alpha2 = d(m.atan((1/self.phi) - m.tan(r(self.beta2))))
        self.alpham = d(m.atan((m.tan(r(self.alpha1))
                                + m.tan(r(self.alpha2)))/2))
        self.cx = self.phi*self.u
        self.w1 = self.cx/m.cos(r(self.beta1))
        self.w2 = self.cx/m.cos(r(self.beta2))
        self.c1 = self.cx/m.cos(r(self.alpha1))
        self.c2 = self.cx/m.cos(r(self.alpha2))
        self.wt1 = self.w1*m.sin(r(self.beta1))
        self.wt2 = self.w1*m.sin(r(self.beta2))
        self.ct1 = self.c1*m.sin(r(self.alpha1))
        self.ct2 = self.c2*m.sin(r(self.alpha2))


################################
# #Class: Stage
# Holds stage characteristics
# #Attributes:
# root: root (Class: ThermoFlow)
# mean: mean (Class: ThermoFlow)
# tip: tip (Class: ThermoFlow)

# rho: air density (kg/m^3)
# mu: air dynamic viscosity (N*s/m^2) (kg/(m*s))
# cp: kJ/(Kg*K)
# cv: kJ/(Kg*K)
# gamma: cp/cv
# R: cp-cv kJ/(Kg*K)

# m: mass flow rate(m^3/s)
# dT0: change in temperature (K)
# PR: pressure ratio
# Cp: coefficient of pressure
# P: power required (kW)
# tau: torque required (Nm)
################################
class Stage():
    """Set stage properties."""

    def __init__(self):
        """Instantiate stage properties."""
        self.root = FanFlow()
        self.mean = FanFlow()
        self.tip = FanFlow()

    def GetStageConfig(self, filename: str):
        """Read .inp files."""
        self.root.GetFanFlowConfig(filename=filename, station='root')
        self.mean.GetFanFlowConfig(filename=filename, station='mean')
        self.tip.GetFanFlowConfig(filename=filename, station='tip')

    def CalcStage(self):
        """Calculate stage properties (assumes cosntant dcx/dr)."""
        self.mean.CalcFanFlow()

        self.root.CalcFanFlow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        rootPhi = self.root.phi
        while (abs(self.mean.cx - self.root.cx) > 1e-3):
            if (self.root.cx < self.mean.cx):
                lphi = rootPhi
            else:
                hphi = rootPhi
            rootPhi = (hphi + lphi)/2
            self.root.CalcFanFlow(iphi=rootPhi)

        self.tip.CalcFanFlow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        tipPhi = self.tip.phi
        while (abs(self.mean.cx - self.tip.cx) > 1e-3):
            if (self.tip.cx < self.mean.cx):
                lphi = tipPhi
            else:
                hphi = tipPhi
            tipPhi = (hphi + lphi)/2
            self.tip.CalcFanFlow(iphi=tipPhi)
